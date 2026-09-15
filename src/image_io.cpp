#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdint>

using namespace Rcpp;

#pragma pack(push, 1)
struct TIFFTag {
  uint16_t tag;
  uint16_t type;
  uint32_t count;
  uint32_t value_offset;
};
#pragma pack(pop)

static bool write_tiff_cpp_impl(const std::string& filename, int w, int h, int channels, const unsigned char* data) {
  FILE* f = fopen(filename.c_str(), "wb");
  if (!f) return false;

  uint32_t spp = (channels >= 3) ? 3 : 1;
  uint32_t total_bytes = static_cast<uint32_t>(w) * h * spp;
  uint32_t ifd_offset = 8 + total_bytes;

  // Header (Little Endian "II")
  uint8_t header[8] = {
    'I', 'I', 42, 0,
    static_cast<uint8_t>(ifd_offset & 0xff),
    static_cast<uint8_t>((ifd_offset >> 8) & 0xff),
    static_cast<uint8_t>((ifd_offset >> 16) & 0xff),
    static_cast<uint8_t>((ifd_offset >> 24) & 0xff)
  };
  fwrite(header, 1, 8, f);

  // Write uncompressed pixel payload
  fwrite(data, 1, total_bytes, f);

  // IFD (9 tags)
  uint16_t num_tags = 9;
  fwrite(&num_tags, 1, 2, f);

  TIFFTag tags[9] = {
    {256, 4, 1, static_cast<uint32_t>(w)},           // ImageWidth
    {257, 4, 1, static_cast<uint32_t>(h)},           // ImageLength
    {258, 3, 1, 8},                                  // BitsPerSample
    {259, 3, 1, 1},                                  // Compression (1 = None)
    {262, 3, 1, static_cast<uint32_t>(spp == 1 ? 1 : 2)}, // Photometric (1 = Gray, 2 = RGB)
    {273, 4, 1, 8},                                  // StripOffsets
    {277, 3, 1, spp},                                // SamplesPerPixel
    {278, 4, 1, static_cast<uint32_t>(h)},           // RowsPerStrip
    {279, 4, 1, total_bytes}                         // StripByteCounts
  };
  fwrite(tags, sizeof(TIFFTag), 9, f);

  uint32_t next_ifd = 0;
  fwrite(&next_ifd, 1, 4, f);

  fclose(f);
  return true;
}

static bool read_tiff_cpp_impl(const std::string& filename, int& w, int& h, int& channels, std::vector<unsigned char>& out_data) {
  FILE* f = fopen(filename.c_str(), "rb");
  if (!f) return false;

  uint8_t hdr[8];
  if (fread(hdr, 1, 8, f) != 8) { fclose(f); return false; }

  bool little_endian = (hdr[0] == 'I' && hdr[1] == 'I');
  uint16_t version = little_endian ? (hdr[2] | (hdr[3] << 8)) : (hdr[3] | (hdr[2] << 8));
  if (version != 42) { fclose(f); return false; }

  uint32_t ifd_offset = little_endian ?
    (hdr[4] | (hdr[5] << 8) | (hdr[6] << 16) | (hdr[7] << 24)) :
    (hdr[7] | (hdr[6] << 8) | (hdr[5] << 16) | (hdr[4] << 24));

  fseek(f, ifd_offset, SEEK_SET);
  uint16_t num_tags = 0;
  if (fread(&num_tags, 1, 2, f) != 2) { fclose(f); return false; }
  if (!little_endian) num_tags = (num_tags >> 8) | (num_tags << 8);

  uint32_t width = 0, height = 0, spp = 1, strip_offset = 0, compression = 1;

  for (uint16_t i = 0; i < num_tags; ++i) {
    TIFFTag tag;
    if (fread(&tag, 1, 12, f) != 12) break;
    if (!little_endian) {
      tag.tag = (tag.tag >> 8) | (tag.tag << 8);
      tag.type = (tag.type >> 8) | (tag.type << 8);
      tag.count = ((tag.count >> 24) & 0xff) | ((tag.count >> 8) & 0xff00) |
                  ((tag.count << 8) & 0xff0000) | ((tag.count << 24) & 0xff000000);
      tag.value_offset = ((tag.value_offset >> 24) & 0xff) | ((tag.value_offset >> 8) & 0xff00) |
                         ((tag.value_offset << 8) & 0xff0000) | ((tag.value_offset << 24) & 0xff000000);
    }

    if (tag.tag == 256) width = tag.value_offset;          // ImageWidth
    else if (tag.tag == 257) height = tag.value_offset;     // ImageLength
    else if (tag.tag == 259) compression = tag.value_offset;// Compression
    else if (tag.tag == 273) strip_offset = tag.value_offset;// StripOffsets
    else if (tag.tag == 277) spp = tag.value_offset;        // SamplesPerPixel
  }

  if (width == 0 || height == 0 || strip_offset == 0 || compression != 1) {
    fclose(f);
    return false;
  }

  w = static_cast<int>(width);
  h = static_cast<int>(height);
  channels = (spp >= 3) ? 3 : 1;

  fseek(f, strip_offset, SEEK_SET);
  std::size_t total_bytes = static_cast<std::size_t>(w) * h * spp;
  out_data.resize(total_bytes);
  size_t read_bytes = fread(out_data.data(), 1, total_bytes, f);
  fclose(f);

  return read_bytes == total_bytes;
}

// [[Rcpp::export]]
SEXP read_image_cpp(std::string filename) {
  int w = 0, h = 0, channels = 0;
  unsigned char *data = nullptr;
  std::vector<unsigned char> tiff_buf;

  std::size_t dot_pos = filename.find_last_of(".");
  std::string ext = (dot_pos != std::string::npos) ? filename.substr(dot_pos + 1) : "";
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  bool is_tiff = (ext == "tif" || ext == "tiff");
  bool tiff_read_success = false;

  if (is_tiff) {
    tiff_read_success = read_tiff_cpp_impl(filename, w, h, channels, tiff_buf);
    if (tiff_read_success) {
      data = tiff_buf.data();
    }
  }

  if (!data) {
    data = stbi_load(filename.c_str(), &w, &h, &channels, 0);
  }

  if (!data) {
    std::string err_msg = stbi_failure_reason() ? stbi_failure_reason() : "Unknown STB/TIFF error";
    Rcpp::stop("Failed to read image '" + filename + "': " + err_msg);
  }

  int n_channels = (channels >= 3) ? 3 : 1;
  std::size_t total_pixels = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
  std::size_t total_size = total_pixels * n_channels;

  RawVector out(total_size);
  unsigned char* p_out = reinterpret_cast<unsigned char*>(out.begin());

  if (n_channels == 3) {
    #pragma omp parallel for collapse(2) if(total_pixels > 50000)
    for (int c = 0; c < 3; ++c) {
      for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
          std::size_t r_idx = x + static_cast<std::size_t>(w) * (y + static_cast<std::size_t>(h) * c);
          std::size_t stb_idx = (static_cast<std::size_t>(y) * w + x) * channels + c;
          p_out[r_idx] = data[stb_idx];
        }
      }
    }
  } else {
    #pragma omp parallel for collapse(2) if(total_pixels > 50000)
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        std::size_t r_idx = x + static_cast<std::size_t>(w) * y;
        std::size_t stb_idx = static_cast<std::size_t>(y) * w + x;
        p_out[r_idx] = data[stb_idx];
      }
    }
  }

  if (!tiff_read_success) {
    stbi_image_free(data);
  }

  if (n_channels == 3) {
    out.attr("dim") = IntegerVector::create(w, h, 3);
  } else {
    out.attr("dim") = IntegerVector::create(w, h);
  }
  out.attr("class") = CharacterVector::create("image", "array");
  out.attr("colormode") = (n_channels == 3) ? "Color" : "Grayscale";

  return out;
}

// [[Rcpp::export]]
bool write_image_cpp(SEXP img_sexp, std::string filename, int quality = 95) {
  RawVector raw_vec;
  if (TYPEOF(img_sexp) == RAWSXP) {
    raw_vec = img_sexp;
  } else {
    NumericVector num_vec(img_sexp);
    std::size_t n = num_vec.size();
    raw_vec = RawVector(n);
    const double* p_in = num_vec.begin();
    unsigned char* p_out = reinterpret_cast<unsigned char*>(raw_vec.begin());

    #pragma omp parallel for schedule(static) if(n > 50000)
    for (std::size_t i = 0; i < n; ++i) {
      double val = p_in[i];
      if (val <= 0.0) {
        p_out[i] = 0;
      } else if (val >= 1.0 || val >= 255.0) {
        p_out[i] = 255;
      } else {
        p_out[i] = static_cast<unsigned char>(val <= 1.0 ? val * 255.0 + 0.5 : val);
      }
    }
    raw_vec.attr("dim") = num_vec.attr("dim");
  }

  IntegerVector dims = raw_vec.attr("dim");
  if (dims.size() < 2) {
    Rcpp::stop("Image provided to write_image_cpp must have 2D or 3D dimensions.");
  }

  int w = dims[0];
  int h = dims[1];
  int channels = (dims.size() >= 3) ? dims[2] : 1;
  int out_channels = (channels >= 3) ? 3 : 1;

  std::size_t total_pixels = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
  std::vector<unsigned char> interleaved(total_pixels * out_channels);
  const unsigned char* p_in = reinterpret_cast<const unsigned char*>(raw_vec.begin());

  if (out_channels == 3) {
    #pragma omp parallel for collapse(2) if(total_pixels > 50000)
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        std::size_t stb_idx = (static_cast<std::size_t>(y) * w + x) * 3;
        for (int c = 0; c < 3; ++c) {
          std::size_t r_idx = x + static_cast<std::size_t>(w) * (y + static_cast<std::size_t>(h) * c);
          interleaved[stb_idx + c] = p_in[r_idx];
        }
      }
    }
  } else {
    #pragma omp parallel for collapse(2) if(total_pixels > 50000)
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        std::size_t stb_idx = static_cast<std::size_t>(y) * w + x;
        std::size_t r_idx = x + static_cast<std::size_t>(w) * y;
        interleaved[stb_idx] = p_in[r_idx];
      }
    }
  }

  std::size_t dot_pos = filename.find_last_of(".");
  std::string ext = (dot_pos != std::string::npos) ? filename.substr(dot_pos + 1) : "jpg";
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  int ret = 0;
  if (ext == "jpg" || ext == "jpeg") {
    ret = stbi_write_jpg(filename.c_str(), w, h, out_channels, interleaved.data(), quality);
  } else if (ext == "png") {
    ret = stbi_write_png(filename.c_str(), w, h, out_channels, interleaved.data(), w * out_channels);
  } else if (ext == "bmp") {
    ret = stbi_write_bmp(filename.c_str(), w, h, out_channels, interleaved.data());
  } else if (ext == "tif" || ext == "tiff") {
    ret = write_tiff_cpp_impl(filename, w, h, out_channels, interleaved.data()) ? 1 : 0;
  } else {
    ret = stbi_write_png(filename.c_str(), w, h, out_channels, interleaved.data(), w * out_channels);
  }

  return ret != 0;
}
