#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "font8x8_basic.h"
#include "font3x5.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

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

// [[Rcpp::export]]
Rcpp::NumericVector preprocess_yolo_file_cpp(
    std::string filename,
    int target_size = 640
) {
  int orig_w = 0, orig_h = 0, orig_ch = 0;
  // stbi_load with req_comp = 3 forces RGB 3-channel interleaved output (R,G,B, R,G,B, ...)
  unsigned char* data = stbi_load(filename.c_str(), &orig_w, &orig_h, &orig_ch, 3);
  if (!data) {
    std::string err = stbi_failure_reason() ? stbi_failure_reason() : "Failed to load image";
    Rcpp::stop("Cannot open image file '%s': %s", filename.c_str(), err.c_str());
  }

  float gain = std::min((float)target_size / (float)orig_w, (float)target_size / (float)orig_h);
  int new_w = std::max(1, (int)std::round(orig_w * gain));
  int new_h = std::max(1, (int)std::round(orig_h * gain));
  float pad_x = ((float)target_size - (float)new_w) / 2.0f;
  float pad_y = ((float)target_size - (float)new_h) / 2.0f;
  int x1 = (int)std::floor(pad_x);
  int y1 = (int)std::floor(pad_y);

  size_t plane_size = (size_t)target_size * target_size;
  size_t total_floats = 3 * plane_size;
  Rcpp::NumericVector out_tensor(total_floats);

  // Fill canvas with YOLO letterbox color: 114 / 255
  double fill_val = 114.0 / 255.0;
  std::fill(out_tensor.begin(), out_tensor.end(), fill_val);
  double* out_ptr = out_tensor.begin();
  struct HorizWeight {
    int x0;
    int x1_idx;
    float dx0;
    float dx1;
    size_t off00;
    size_t off10;
  };
  std::vector<HorizWeight> hw(new_w);
  for (int x = 0; x < new_w; ++x) {
    float src_x = ((float)x + 0.5f) / gain - 0.5f;
    if (src_x < 0.0f) src_x = 0.0f;
    if (src_x > (float)(orig_w - 1)) src_x = (float)(orig_w - 1);
    int x0 = (int)std::floor(src_x);
    int x1_idx = std::min(x0 + 1, orig_w - 1);
    float dx = src_x - (float)x0;
    hw[x].x0 = x0;
    hw[x].x1_idx = x1_idx;
    hw[x].dx0 = 1.0f - dx;
    hw[x].dx1 = dx;
    hw[x].off00 = (size_t)x0 * 3;
    hw[x].off10 = (size_t)x1_idx * 3;
  }

  size_t row_stride = (size_t)orig_w * 3;

  #pragma omp parallel for schedule(static) if(new_h > 30)
  for (int y = 0; y < new_h; ++y) {
    float src_y = ((float)y + 0.5f) / gain - 0.5f;
    if (src_y < 0.0f) src_y = 0.0f;
    if (src_y > (float)(orig_h - 1)) src_y = (float)(orig_h - 1);
    int y0 = (int)std::floor(src_y);
    int y1_idx = std::min(y0 + 1, orig_h - 1);
    float dy1 = src_y - (float)y0;
    float dy0 = 1.0f - dy1;

    int dst_y = y1 + y;
    size_t dst_base = (size_t)dst_y * target_size + x1;
    const unsigned char* row0 = data + (size_t)y0 * row_stride;
    const unsigned char* row1 = data + (size_t)y1_idx * row_stride;

    for (int x = 0; x < new_w; ++x) {
      const HorizWeight& w = hw[x];
      float w00 = w.dx0 * dy0;
      float w10 = w.dx1 * dy0;
      float w01 = w.dx0 * dy1;
      float w11 = w.dx1 * dy1;

      size_t dst_idx = dst_base + x;

      for (int c = 0; c < 3; ++c) {
        float val = (w00 * row0[w.off00 + c] +
                     w10 * row0[w.off10 + c] +
                     w01 * row1[w.off00 + c] +
                     w11 * row1[w.off10 + c]) / 255.0f;
        out_ptr[(size_t)c * plane_size + dst_idx] = (double)val;
      }
    }
  }

  stbi_image_free(data);

  out_tensor.attr("gain") = gain;
  out_tensor.attr("pad_x") = pad_x;
  out_tensor.attr("pad_y") = pad_y;
  out_tensor.attr("orig_w") = orig_w;
  out_tensor.attr("orig_h") = orig_h;

  return out_tensor;
}

static inline void set_pixel_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, uint8_t r, uint8_t g, uint8_t b) {
  if (x >= 0 && x < orig_w && y >= 0 && y < orig_h) {
    size_t idx = (size_t)x * 3 + (size_t)y * ((size_t)3 * orig_w);
    p[idx]     = r;
    p[idx + 1] = g;
    p[idx + 2] = b;
  }
}

static inline void blend_pixel_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, uint8_t r, uint8_t g, uint8_t b, float alpha) {
  if (x >= 0 && x < orig_w && y >= 0 && y < orig_h) {
    size_t idx = (size_t)x * 3 + (size_t)y * ((size_t)3 * orig_w);
    int a = (int)std::round(alpha * 256.0f);
    int inv_a = 256 - a;
    p[idx]     = (uint8_t)(((int)p[idx]     * inv_a + (int)r * a) >> 8);
    p[idx + 1] = (uint8_t)(((int)p[idx + 1] * inv_a + (int)g * a) >> 8);
    p[idx + 2] = (uint8_t)(((int)p[idx + 2] * inv_a + (int)b * a) >> 8);
  }
}

static void fill_rect_rgb(unsigned char* p, int orig_w, int orig_h, int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b) {
  x1 = std::max(0, std::min(x1, orig_w - 1));
  x2 = std::max(0, std::min(x2, orig_w - 1));
  y1 = std::max(0, std::min(y1, orig_h - 1));
  y2 = std::max(0, std::min(y2, orig_h - 1));
  if (x1 > x2) std::swap(x1, x2);
  if (y1 > y2) std::swap(y1, y2);

  for (int y = y1; y <= y2; ++y) {
    size_t row_start = (size_t)y * ((size_t)3 * orig_w);
    for (int x = x1; x <= x2; ++x) {
      size_t idx = (size_t)x * 3 + row_start;
      p[idx]     = r;
      p[idx + 1] = g;
      p[idx + 2] = b;
    }
  }
}

static void fill_rect_alpha_rgb(unsigned char* p, int orig_w, int orig_h, int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b, float alpha = 0.3f) {
  x1 = std::max(0, std::min(x1, orig_w - 1));
  x2 = std::max(0, std::min(x2, orig_w - 1));
  y1 = std::max(0, std::min(y1, orig_h - 1));
  y2 = std::max(0, std::min(y2, orig_h - 1));
  if (x1 > x2) std::swap(x1, x2);
  if (y1 > y2) std::swap(y1, y2);

  int a = (int)std::round(std::max(0.0f, std::min(1.0f, alpha)) * 256.0f);
  int inv_a = 256 - a;

  for (int y = y1; y <= y2; ++y) {
    size_t row_start = (size_t)y * ((size_t)3 * orig_w);
    for (int x = x1; x <= x2; ++x) {
      size_t idx = (size_t)x * 3 + row_start;
      p[idx]     = (uint8_t)(((int)p[idx]     * inv_a + (int)r * a) >> 8);
      p[idx + 1] = (uint8_t)(((int)p[idx + 1] * inv_a + (int)g * a) >> 8);
      p[idx + 2] = (uint8_t)(((int)p[idx + 2] * inv_a + (int)b * a) >> 8);
    }
  }
}

static void draw_box_rgb(unsigned char* p, int orig_w, int orig_h, int x1, int y1, int x2, int y2, int lwd, uint8_t r, uint8_t g, uint8_t b) {
  for (int t = 0; t < lwd; ++t) {
    fill_rect_rgb(p, orig_w, orig_h, x1, y1 + t, x2, y1 + t, r, g, b);
    fill_rect_rgb(p, orig_w, orig_h, x1, y2 - t, x2, y2 - t, r, g, b);
    fill_rect_rgb(p, orig_w, orig_h, x1 + t, y1, x1 + t, y2, r, g, b);
    fill_rect_rgb(p, orig_w, orig_h, x2 - t, y1, x2 - t, y2, r, g, b);
  }
}

static void draw_char_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, char ch, int scale, uint8_t r, uint8_t g, uint8_t b) {
  uint8_t c = (uint8_t)ch;
  if (c > 127) c = '?';
  const unsigned char* glyph = font8x8_basic[c];

  for (int row = 0; row < 8; ++row) {
    unsigned char byte = glyph[row];
    for (int col = 0; col < 8; ++col) {
      if ((byte >> col) & 1) {
        if (scale == 1) {
          set_pixel_rgb(p, orig_w, orig_h, x + col, y + row, r, g, b);
        } else {
          fill_rect_rgb(p, orig_w, orig_h, x + col * scale, y + row * scale, x + (col + 1) * scale - 1, y + (row + 1) * scale - 1, r, g, b);
        }
      }
    }
  }
}

static void draw_text_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, const std::string& text, int scale, uint8_t r, uint8_t g, uint8_t b) {
  int cur_x = x;
  for (size_t i = 0; i < text.size(); ++i) {
    draw_char_rgb(p, orig_w, orig_h, cur_x, y, text[i], scale, r, g, b);
    cur_x += 8 * scale;
  }
}

static void draw_char_micro_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, char ch, uint8_t r, uint8_t g, uint8_t b) {
  uint8_t c = (uint8_t)ch;
  if (c > 127) c = '?';
  const unsigned char* glyph = font3x5_basic[c];
  for (int row = 0; row < 5; ++row) {
    unsigned char byte = glyph[row];
    for (int col = 0; col < 3; ++col) {
      if ((byte >> (2 - col)) & 1) {
        set_pixel_rgb(p, orig_w, orig_h, x + col, y + row, r, g, b);
      }
    }
  }
}

static void draw_text_micro_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, const std::string& text, uint8_t r, uint8_t g, uint8_t b) {
  int cur_x = x;
  for (size_t i = 0; i < text.size(); ++i) {
    draw_char_micro_rgb(p, orig_w, orig_h, cur_x, y, text[i], r, g, b);
    cur_x += 4; // 3px glyph + 1px spacing
  }
}

static inline void parse_color_rgb(const std::string& hex, uint8_t& r, uint8_t& g, uint8_t& b) {
  if (hex.size() >= 7 && hex[0] == '#') {
    unsigned int rv = 0, gv = 0, bv = 0;
    if (std::sscanf(hex.c_str() + 1, "%02x%02x%02x", &rv, &gv, &bv) == 3) {
      r = (uint8_t)rv;
      g = (uint8_t)gv;
      b = (uint8_t)bv;
      return;
    }
  }
  // Default green #00CC66
  r = 0; g = 204; b = 102;
}

static void draw_line_rgb(unsigned char* p, int orig_w, int orig_h, int x0, int y0, int x1, int y1, int lwd, uint8_t r, uint8_t g, uint8_t b) {
  int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy, e2;

  while (true) {
    if (lwd <= 1) {
      set_pixel_rgb(p, orig_w, orig_h, x0, y0, r, g, b);
    } else {
      int half = lwd / 2;
      fill_rect_rgb(p, orig_w, orig_h, x0 - half, y0 - half, x0 + half, y0 + half, r, g, b);
    }
    if (x0 == x1 && y0 == y1) break;
    e2 = 2 * err;
    if (e2 >= dy) { err += dy; x0 += sx; }
    if (e2 <= dx) { err += dx; y0 += sy; }
  }
}

static void draw_circle_rgb(unsigned char* p, int orig_w, int orig_h, int cx, int cy, int radius, uint8_t r, uint8_t g, uint8_t b) {
  int r2 = radius * radius;
  for (int dy = -radius; dy <= radius; ++dy) {
    int y = cy + dy;
    if (y < 0 || y >= orig_h) continue;
    for (int dx = -radius; dx <= radius; ++dx) {
      int x = cx + dx;
      if (x < 0 || x >= orig_w) continue;
      if (dx * dx + dy * dy <= r2) {
        set_pixel_rgb(p, orig_w, orig_h, x, y, r, g, b);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// MODERN PROPORTIONAL TYPOGRAPHY & HUD TELEMETRY
// ---------------------------------------------------------------------------

static inline int get_char_width_modern(char ch) {
  uint8_t c = (uint8_t)ch;
  if (c > 127) return 6;
  if (c == ' ') return 4;
  if (c == ':' || c == '.' || c == ',' || c == '\'' || c == '!') return 3;
  if (c == '|' || c == '(' || c == ')' || c == '[' || c == ']') return 4;
  if (c == '1' || c == 'I' || c == 'i' || c == 'l') return 4;
  const unsigned char* glyph = font8x8_basic[c];
  int min_col = 7, max_col = 0;
  bool has_pixel = false;
  for (int row = 0; row < 8; ++row) {
    unsigned char byte = glyph[row];
    for (int col = 0; col < 8; ++col) {
      if ((byte >> col) & 1) {
        has_pixel = true;
        if (col < min_col) min_col = col;
        if (col > max_col) max_col = col;
      }
    }
  }
  if (!has_pixel) return 4;
  return std::max(3, max_col - min_col + 1);
}

static inline int get_text_width_modern(const std::string& text, int scale = 1) {
  int w = 0;
  for (size_t i = 0; i < text.size(); ++i) {
    w += (get_char_width_modern(text[i]) + 1) * scale;
  }
  return w;
}

static int draw_char_modern_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, char ch, int scale, uint8_t r, uint8_t g, uint8_t b, bool shadow = true) {
  uint8_t c = (uint8_t)ch;
  if (c > 127) c = '?';
  if (c == ' ') {
    return 4 * scale;
  }
  const unsigned char* glyph = font8x8_basic[c];

  int min_col = 0, max_col = 7;
  bool found = false;
  for (int col = 0; col < 8; ++col) {
    for (int row = 0; row < 8; ++row) {
      if ((glyph[row] >> col) & 1) {
        if (!found) { min_col = col; found = true; }
        max_col = col;
      }
    }
  }
  if (!found) return 4 * scale;
  int char_w = max_col - min_col + 1;

  // 1. Soft drop shadow (1px offset)
  if (shadow) {
    for (int row = 0; row < 8; ++row) {
      unsigned char byte = glyph[row];
      for (int col = min_col; col <= max_col; ++col) {
        if ((byte >> col) & 1) {
          int sx = x + (col - min_col) * scale + 1;
          int sy = y + row * scale + 1;
          if (scale == 1) {
            blend_pixel_rgb(p, orig_w, orig_h, sx, sy, 0, 0, 0, 0.65f);
          } else {
            fill_rect_alpha_rgb(p, orig_w, orig_h, sx, sy, sx + scale - 1, sy + scale - 1, 0, 0, 0, 0.65f);
          }
        }
      }
    }
  }

  // 2. Crisp foreground glyph
  for (int row = 0; row < 8; ++row) {
    unsigned char byte = glyph[row];
    for (int col = min_col; col <= max_col; ++col) {
      if ((byte >> col) & 1) {
        int px = x + (col - min_col) * scale;
        int py = y + row * scale;
        if (scale == 1) {
          set_pixel_rgb(p, orig_w, orig_h, px, py, r, g, b);
        } else {
          fill_rect_rgb(p, orig_w, orig_h, px, py, px + scale - 1, py + scale - 1, r, g, b);
        }
      }
    }
  }

  return (char_w + 1) * scale;
}

static void draw_text_modern_rgb(unsigned char* p, int orig_w, int orig_h, int x, int y, const std::string& text, int scale, uint8_t r, uint8_t g, uint8_t b, bool shadow = true) {
  int cur_x = x;
  for (size_t i = 0; i < text.size(); ++i) {
    cur_x += draw_char_modern_rgb(p, orig_w, orig_h, cur_x, y, text[i], scale, r, g, b, shadow);
  }
}

struct HudTelemetryItem {
  std::string label;
  std::string value;
  uint8_t val_r, val_g, val_b;
};

static inline std::string trim_str(const std::string& s) {
  size_t start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  size_t end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

static void draw_hud_modern(
    unsigned char* p,
    int orig_w,
    int orig_h,
    const std::string& hud_text,
    const std::string& hud_pos,
    const std::string& hud_layout,
    double font_scale = 1.0
) {
  if (hud_text.empty()) return;

  // 1. Split tokens by '|' or '\n'
  std::vector<HudTelemetryItem> items;
  std::string cur_tok = "";
  for (size_t i = 0; i <= hud_text.size(); ++i) {
    char ch = (i < hud_text.size()) ? hud_text[i] : '|';
    if (ch == '|' || ch == '\n') {
      std::string t = trim_str(cur_tok);
      cur_tok = "";
      if (t.empty()) continue;

      HudTelemetryItem item;
      std::string t_lower = t;
      for (size_t k = 0; k < t_lower.size(); ++k) t_lower[k] = (char)::tolower(t_lower[k]);

      if (t_lower.rfind("frame", 0) == 0) {
        item.label = "FRAME";
        size_t sp = t.find(' ');
        item.value = (sp != std::string::npos) ? trim_str(t.substr(sp + 1)) : t;
        item.val_r = 255; item.val_g = 255; item.val_b = 255; // Crisp White
      } else if (t_lower.find("fps") != std::string::npos) {
        item.label = "SPEED";
        item.value = t;
        item.val_r = 0; item.val_g = 255; item.val_b = 136; // Neon Emerald
      } else if (t_lower.find("in zone") != std::string::npos || t_lower.find("zone:") != std::string::npos) {
        item.label = "IN ZONE";
        size_t cp = t.find(':');
        item.value = (cp != std::string::npos) ? trim_str(t.substr(cp + 1)) : t;
        item.val_r = 255; item.val_g = 187; item.val_b = 0; // Glowing Amber
      } else if (t_lower.find("count") != std::string::npos) {
        item.label = "COUNTED";
        size_t cp = t.find(':');
        item.value = (cp != std::string::npos) ? trim_str(t.substr(cp + 1)) : t;
        item.val_r = 0; item.val_g = 229; item.val_b = 255; // Electric Cyan
      } else if (t_lower.find("det") != std::string::npos) {
        item.label = "OBJECTS";
        size_t dp = t.find(" Det");
        if (dp == std::string::npos) dp = t.find(" det");
        item.value = (dp != std::string::npos) ? trim_str(t.substr(0, dp)) : t;
        item.val_r = 255; item.val_g = 255; item.val_b = 255;
      } else if (t.find(':') != std::string::npos) {
        size_t cp = t.find(':');
        item.label = trim_str(t.substr(0, cp));
        item.value = trim_str(t.substr(cp + 1));
        item.val_r = 255; item.val_g = 255; item.val_b = 255;
      } else {
        item.label = t;
        item.value = "";
        item.val_r = 255; item.val_g = 255; item.val_b = 255;
      }
      items.push_back(item);
    } else {
      cur_tok += ch;
    }
  }

  if (items.empty()) return;

  int scale = (font_scale >= 1.5) ? (int)std::round(font_scale) : 1;
  int card_w = 0, card_h = 0;
  bool is_vertical = (hud_layout != "horizontal" && hud_layout != "horiz");

  if (is_vertical) {
    int max_label_w = 0, max_val_w = 0;
    for (size_t k = 0; k < items.size(); ++k) {
      int lw = get_text_width_modern(items[k].label, scale);
      int vw = get_text_width_modern(items[k].value, scale);
      if (lw > max_label_w) max_label_w = lw;
      if (vw > max_val_w) max_val_w = vw;
    }
    card_w = std::max(130 * scale, max_label_w + max_val_w + 30 * scale);
    card_h = (26 + (int)items.size() * 16 + 6) * scale;
  } else {
    // Horizontal layout
    int content_w = 26 * scale; // dot + padding
    for (size_t k = 0; k < items.size(); ++k) {
      int lw = get_text_width_modern(items[k].label + ":", scale);
      int vw = get_text_width_modern(items[k].value, scale);
      content_w += lw + 4 * scale + vw + 14 * scale;
    }
    card_w = content_w;
    card_h = 24 * scale;
  }

  // 2. Position calculation
  int margin = 14;
  std::string pos_str = hud_pos;
  for (size_t k = 0; k < pos_str.size(); ++k) pos_str[k] = (char)::tolower(pos_str[k]);

  int x1 = 0, y1 = 0, x2 = 0, y2 = 0;
  if (pos_str == "top-left" || pos_str == "topleft") {
    x1 = margin;
    y1 = margin;
    x2 = x1 + card_w;
    y2 = y1 + card_h;
  } else if (pos_str == "bottom-right" || pos_str == "bottomright") {
    x2 = orig_w - margin;
    x1 = x2 - card_w;
    y2 = orig_h - margin;
    y1 = y2 - card_h;
  } else if (pos_str == "bottom-left" || pos_str == "bottomleft") {
    x1 = margin;
    x2 = x1 + card_w;
    y2 = orig_h - margin;
    y1 = y2 - card_h;
  } else if (pos_str == "top") {
    x1 = (orig_w - card_w) / 2;
    x2 = x1 + card_w;
    y1 = margin;
    y2 = y1 + card_h;
  } else if (pos_str == "bottom") {
    x1 = (orig_w - card_w) / 2;
    x2 = x1 + card_w;
    y2 = orig_h - margin;
    y1 = y2 - card_h;
  } else {
    // Default to "top-right"
    x2 = orig_w - margin;
    x1 = x2 - card_w;
    y1 = margin;
    y2 = y1 + card_h;
  }

  // Clamp within image bounds
  x1 = std::max(0, x1); y1 = std::max(0, y1);
  x2 = std::min(orig_w - 1, x2); y2 = std::min(orig_h - 1, y2);

  // 3. Render Card Backdrop
  // Dark Obsidian Slate with 82% alpha
  fill_rect_alpha_rgb(p, orig_w, orig_h, x1, y1, x2, y2, 12, 18, 28, 0.82f);
  // 1px subtle tech border
  draw_box_rgb(p, orig_w, orig_h, x1, y1, x2, y2, 1, 35, 60, 85);
  // Glowing cyan top accent line
  fill_rect_rgb(p, orig_w, orig_h, x1, y1, x2, y1 + std::max(2, scale), 0, 229, 255);

  // 4. Render Content
  if (is_vertical) {
    // Header with glowing status dot
    draw_circle_rgb(p, orig_w, orig_h, x1 + 10 * scale, y1 + 12 * scale, 2 * scale, 0, 255, 136);
    draw_text_modern_rgb(p, orig_w, orig_h, x1 + 16 * scale, y1 + 8 * scale, "TELEMETRY", scale, 0, 229, 255);
    // Subtle horizontal divider line
    draw_line_rgb(p, orig_w, orig_h, x1 + 8 * scale, y1 + 21 * scale, x2 - 8 * scale, y1 + 21 * scale, 1, 35, 55, 78);

    // Rows
    for (size_t k = 0; k < items.size(); ++k) {
      int cur_y = y1 + (26 + (int)k * 16) * scale;
      // Label in cool slate
      draw_text_modern_rgb(p, orig_w, orig_h, x1 + 10 * scale, cur_y, items[k].label, scale, 148, 176, 196);
      // Value right-aligned in its high-contrast color
      if (!items[k].value.empty()) {
        int vw = get_text_width_modern(items[k].value, scale);
        draw_text_modern_rgb(p, orig_w, orig_h, x2 - 10 * scale - vw, cur_y, items[k].value, scale, items[k].val_r, items[k].val_g, items[k].val_b);
      }
    }
  } else {
    // Horizontal layout
    draw_circle_rgb(p, orig_w, orig_h, x1 + 10 * scale, y1 + 12 * scale, 2 * scale, 0, 255, 136);
    int cur_x = x1 + 18 * scale;
    for (size_t k = 0; k < items.size(); ++k) {
      if (k > 0) {
        draw_line_rgb(p, orig_w, orig_h, cur_x - 6 * scale, y1 + 5 * scale, cur_x - 6 * scale, y2 - 5 * scale, 1, 40, 65, 90);
      }
      std::string lbl_str = items[k].label + ":";
      draw_text_modern_rgb(p, orig_w, orig_h, cur_x, y1 + 8 * scale, lbl_str, scale, 148, 176, 196);
      cur_x += get_text_width_modern(lbl_str, scale) + 4 * scale;
      draw_text_modern_rgb(p, orig_w, orig_h, cur_x, y1 + 8 * scale, items[k].value, scale, items[k].val_r, items[k].val_g, items[k].val_b);
      cur_x += get_text_width_modern(items[k].value, scale) + 12 * scale;
    }
  }
}

// 16 COCO anatomical skeleton pairs (0-indexed)
// 0: nose, 1: left_eye, 2: right_eye, 3: left_ear, 4: right_ear,
// 5: left_shoulder, 6: right_shoulder, 7: left_elbow, 8: right_elbow,
// 9: left_wrist, 10: right_wrist, 11: left_hip, 12: right_hip,
// 13: left_knee, 14: right_knee, 15: left_ankle, 16: right_ankle
static const int SKELETON_PAIRS[16][2] = {
  {0, 1}, {0, 2}, {1, 3}, {2, 4},     // 0..3: facial (red)
  {5, 6}, {5, 11}, {6, 12}, {11, 12}, // 4..7: torso (orange)
  {5, 7}, {7, 9},                     // 8..9: left arm (green)
  {6, 8}, {8, 10},                    // 10..11: right arm (blue)
  {11, 13}, {13, 15},                 // 12..13: left leg (purple)
  {12, 14}, {14, 16}                  // 14..15: right leg (magenta)
};

static const uint8_t LIMB_COLORS[16][3] = {
  {255, 75, 75},   {255, 75, 75},   {255, 75, 75},   {255, 75, 75},   // facial (red)
  {255, 165, 0},   {255, 165, 0},   {255, 165, 0},   {255, 165, 0},   // torso (orange)
  {0, 204, 102},   {0, 204, 102},                                     // left arm (green)
  {0, 153, 255},   {0, 153, 255},                                     // right arm (blue)
  {153, 51, 255},  {153, 51, 255},                                    // left leg (purple)
  {255, 0, 255},   {255, 0, 255}                                      // right leg (magenta)
};

// [[Rcpp::export]]
bool draw_yolo_detections_bgr_cpp(
    Rcpp::RawVector bm,
    Rcpp::NumericVector xmin,
    Rcpp::NumericVector ymin,
    Rcpp::NumericVector xmax,
    Rcpp::NumericVector ymax,
    Rcpp::CharacterVector labels,
    Rcpp::NumericVector scores,
    Rcpp::CharacterVector colors,
    int lwd = 2,
    std::string hud_text = "",
    double font_scale = 1.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> keypoints = R_NilValue,
    double kpt_threshold = 0.3,
    int kpt_radius = 4,
    bool draw_skeleton = true,
    bool draw_boxes = true,
    Rcpp::Nullable<Rcpp::IntegerVector> track_ids = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> flash = R_NilValue,
    Rcpp::Nullable<Rcpp::List> history_x = R_NilValue,
    Rcpp::Nullable<Rcpp::List> history_y = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> roi = R_NilValue,
    std::string roi_label = "",
    Rcpp::Nullable<Rcpp::NumericVector> count_line = R_NilValue,
    std::string count_line_label = "",
    Rcpp::Nullable<Rcpp::LogicalVector> counted = R_NilValue,
    bool hide_outside_roi = true,
    bool hide_counted = true,
    bool flash_counted = true,
    bool show_text = true,
    bool show_conf = true,
    bool show_class = true,
    bool show_id = true,
    Rcpp::Nullable<Rcpp::IntegerMatrix> mask_labels = R_NilValue,
    double mask_alpha = 0.4,
    int mask_offset_x = 0,
    int mask_offset_y = 0,
    bool draw_masks = true,
    Rcpp::Nullable<Rcpp::IntegerVector> mask_ids = R_NilValue,
    std::string hud_pos = "top-right",
    std::string hud_layout = "vertical"
) {
  SEXP dim_attr = Rf_getAttrib(bm, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 3) return false;
  int* dims = INTEGER(dim_attr);
  int orig_w = dims[1];
  int orig_h = dims[2];
  unsigned char* p = RAW(bm);

  // 1. Draw ROI box if specified
  bool has_roi_bounds = false;
  int rx1 = 0, ry1 = 0, rx2 = 0, ry2 = 0;
  if (roi.isNotNull()) {
    Rcpp::NumericVector rbox(roi.get());
    if (rbox.size() >= 4) {
      has_roi_bounds = true;
      rx1 = (int)std::round(rbox[0]);
      ry1 = (int)std::round(rbox[1]);
      rx2 = (int)std::round(rbox[2]);
      ry2 = (int)std::round(rbox[3]);
      if (rx1 > rx2) std::swap(rx1, rx2);
      if (ry1 > ry2) std::swap(ry1, ry2);

      // Semi-transparent salmon fill (RGB: 250, 128, 114, alpha: 0.3)
      fill_rect_alpha_rgb(p, orig_w, orig_h, rx1, ry1, rx2, ry2, 250, 128, 114, 0.3f);
      draw_box_rgb(p, orig_w, orig_h, rx1, ry1, rx2, ry2, 1, 250, 128, 114);

      std::string rlabel = !roi_label.empty() ? roi_label : "ROI / MONITORED ZONE";
      int rw = (int)rlabel.size() * 8;
      int r_y = std::max(0, ry1 - 16);
      if (ry1 < 16) {
        fill_rect_rgb(p, orig_w, orig_h, rx1, ry1, rx1 + rw + 8, ry1 + 16, 25, 25, 25);
        draw_box_rgb(p, orig_w, orig_h, rx1, ry1, rx1 + rw + 8, ry1 + 16, 1, 250, 128, 114);
        draw_text_rgb(p, orig_w, orig_h, rx1 + 4, ry1 + 4, rlabel, 1, 250, 128, 114);
      } else {
        fill_rect_rgb(p, orig_w, orig_h, rx1, r_y, rx1 + rw + 8, ry1, 25, 25, 25);
        draw_box_rgb(p, orig_w, orig_h, rx1, r_y, rx1 + rw + 8, ry1, 1, 250, 128, 114);
        draw_text_rgb(p, orig_w, orig_h, rx1 + 4, r_y + 4, rlabel, 1, 250, 128, 114);
      }
    }
  }

  // 2. Draw virtual counting line if specified
  if (count_line.isNotNull()) {
    Rcpp::NumericVector cline(count_line.get());
    if (cline.size() >= 4) {
      int lx1 = (int)std::round(cline[0]);
      int ly1 = (int)std::round(cline[1]);
      int lx2 = (int)std::round(cline[2]);
      int ly2 = (int)std::round(cline[3]);
      draw_line_rgb(p, orig_w, orig_h, lx1, ly1, lx2, ly2, 3, 0, 229, 255);
      draw_circle_rgb(p, orig_w, orig_h, lx1, ly1, 5, 0, 229, 255);
      draw_circle_rgb(p, orig_w, orig_h, lx2, ly2, 5, 0, 229, 255);

      if (!count_line_label.empty()) {
        std::string clabel = count_line_label;
        int lw = (int)clabel.size() * 8;

        // Position label at the border/extremity instead of the center
        int tag_x1 = 0, tag_y1 = 0, tag_x2 = 0, tag_y2 = 0;
        bool is_vertical = std::abs(ly2 - ly1) >= std::abs(lx2 - lx1);

        if (is_vertical) {
          int top_y = std::min(ly1, ly2);
          tag_y1 = (top_y < 18) ? 0 : (top_y - 18);
          tag_y2 = tag_y1 + 16;
          int lx = (ly1 <= ly2) ? lx1 : lx2;
          if (lx - lw - 8 >= 0) {
            tag_x2 = lx;
            tag_x1 = tag_x2 - lw - 8;
          } else {
            tag_x1 = lx;
            tag_x2 = std::min(orig_w - 1, tag_x1 + lw + 8);
          }
        } else {
          int left_x = std::min(lx1, lx2);
          tag_x1 = std::max(0, left_x);
          tag_x2 = std::min(orig_w - 1, tag_x1 + lw + 8);
          int ly = (lx1 <= lx2) ? ly1 : ly2;
          tag_y1 = (ly < 18) ? ly : (ly - 18);
          tag_y2 = tag_y1 + 16;
        }

        fill_rect_rgb(p, orig_w, orig_h, tag_x1, tag_y1, tag_x2, tag_y2, 25, 25, 25);
        draw_box_rgb(p, orig_w, orig_h, tag_x1, tag_y1, tag_x2, tag_y2, 1, 0, 229, 255);
        draw_text_rgb(p, orig_w, orig_h, tag_x1 + 4, tag_y1 + 4, clabel, 1, 0, 229, 255);
      }
    }
  }

  int n_dets = xmin.size();

  bool has_kpts_input = false;
  Rcpp::NumericMatrix kpt_mat;
  if (keypoints.isNotNull()) {
    kpt_mat = Rcpp::as<Rcpp::NumericMatrix>(keypoints);
    if (kpt_mat.nrow() == n_dets && kpt_mat.ncol() >= 51) {
      has_kpts_input = true;
    }
  }

  bool has_tracks = track_ids.isNotNull();
  Rcpp::IntegerVector t_ids;
  if (has_tracks) t_ids = Rcpp::as<Rcpp::IntegerVector>(track_ids);

  bool has_flash = flash.isNotNull();
  Rcpp::IntegerVector flash_vec;
  if (has_flash) flash_vec = Rcpp::as<Rcpp::IntegerVector>(flash);

  bool has_history = (history_x.isNotNull() && history_y.isNotNull());
  Rcpp::List hx_list, hy_list;
  if (has_history) {
    hx_list = Rcpp::as<Rcpp::List>(history_x);
    hy_list = Rcpp::as<Rcpp::List>(history_y);
  }

  bool has_counted = counted.isNotNull();
  Rcpp::LogicalVector counted_vec;
  if (has_counted) counted_vec = Rcpp::as<Rcpp::LogicalVector>(counted);

  bool has_masks_input = draw_masks && mask_labels.isNotNull();
  Rcpp::IntegerMatrix ml;
  const int* m_ptr = NULL;
  int mw = 0, mh = 0;
  int alpha_val = (int)std::round(std::max(0.0, std::min(1.0, mask_alpha)) * 256.0);
  int inv_alpha = 256 - alpha_val;
  if (has_masks_input) {
    ml = Rcpp::as<Rcpp::IntegerMatrix>(mask_labels);
    mw = ml.nrow();
    mh = ml.ncol();
    m_ptr = INTEGER(ml);
  }

  bool has_mids = mask_ids.isNotNull();
  Rcpp::IntegerVector mids_vec;
  if (has_mids) mids_vec = Rcpp::as<Rcpp::IntegerVector>(mask_ids);

  // Draw semi-transparent filled instance segmentation masks (no contours)
  if (has_masks_input && alpha_val > 0 && m_ptr != NULL) {
    for (int i = 0; i < n_dets; ++i) {
      int x1 = (int)std::round(xmin[i]);
      int y1 = (int)std::round(ymin[i]);
      int x2 = (int)std::round(xmax[i]);
      int y2 = (int)std::round(ymax[i]);

      if (hide_outside_roi && has_roi_bounds) {
        int cx = (x1 + x2) / 2;
        int cy = (y1 + y2) / 2;
        if (cx < rx1 || cx > rx2 || cy < ry1 || cy > ry2) continue;
      }
      if (hide_counted && has_counted && i < counted_vec.size() && counted_vec[i]) {
        int cur_flash = (has_flash && i < flash_vec.size()) ? flash_vec[i] : 0;
        if (!flash_counted || cur_flash <= 0) continue;
      }

      uint8_t mr = 0, mg = 204, mb = 102;
      if (colors.size() > 0) {
        std::string c_str = Rcpp::as<std::string>(colors[i % colors.size()]);
        parse_color_rgb(c_str, mr, mg, mb);
      }

      int target_lbl = (has_mids && i < mids_vec.size()) ? mids_vec[i] : (i + 1);

      int bx1 = std::max(0, std::min(orig_w - 1, x1));
      int bx2 = std::max(0, std::min(orig_w - 1, x2));
      int by1 = std::max(0, std::min(orig_h - 1, y1));
      int by2 = std::max(0, std::min(orig_h - 1, y2));

      for (int y = by1; y <= by2; ++y) {
        int my = y - mask_offset_y;
        if (my < 0 || my >= mh) continue;
        size_t row_start = (size_t)y * ((size_t)3 * orig_w);
        size_t m_col_start = (size_t)my * mw;

        for (int x = bx1; x <= bx2; ++x) {
          int mx = x - mask_offset_x;
          if (mx < 0 || mx >= mw) continue;

          if (m_ptr[m_col_start + mx] == target_lbl) {
            size_t idx = (size_t)x * 3 + row_start;
            p[idx]     = (uint8_t)(((int)p[idx]     * inv_alpha + (int)mr * alpha_val) >> 8);
            p[idx + 1] = (uint8_t)(((int)p[idx + 1] * inv_alpha + (int)mg * alpha_val) >> 8);
            p[idx + 2] = (uint8_t)(((int)p[idx + 2] * inv_alpha + (int)mb * alpha_val) >> 8);
          }
        }
      }
    }
  }

  for (int i = 0; i < n_dets; ++i) {
    int x1 = (int)std::round(xmin[i]);
    int y1 = (int)std::round(ymin[i]);
    int x2 = (int)std::round(xmax[i]);
    int y2 = (int)std::round(ymax[i]);

    // 1. Suppress if outside ROI
    if (hide_outside_roi && has_roi_bounds) {
      int cx = (x1 + x2) / 2;
      int cy = (y1 + y2) / 2;
      if (cx < rx1 || cx > rx2 || cy < ry1 || cy > ry2) {
        continue;
      }
    }

    // 2. Suppress if already crossed virtual counting line
    int cur_flash = (has_flash && i < flash_vec.size()) ? flash_vec[i] : 0;
    if (hide_counted && has_counted && i < counted_vec.size() && counted_vec[i]) {
      if (!flash_counted || cur_flash <= 0) {
        continue;
      }
    }

    uint8_t r = 0, g = 204, b = 102;
    if (colors.size() > 0) {
      std::string c_str = Rcpp::as<std::string>(colors[i % colors.size()]);
      parse_color_rgb(c_str, r, g, b);
    }

    // Flash bright neon green if just crossed counting line
    int cur_lwd = lwd;
    if (has_flash && i < flash_vec.size() && flash_vec[i] > 0) {
      r = 0; g = 255; b = 128;
      cur_lwd = lwd + 1;
    }

    // Draw motion trail behind object
    if (has_history && i < hx_list.size() && i < hy_list.size()) {
      Rcpp::NumericVector hx = hx_list[i];
      Rcpp::NumericVector hy = hy_list[i];
      int h_len = hx.size();
      for (int h = 0; h < h_len - 1; ++h) {
        int px0 = (int)std::round(hx[h]);
        int py0 = (int)std::round(hy[h]);
        int px1 = (int)std::round(hx[h + 1]);
        int py1 = (int)std::round(hy[h + 1]);
        draw_line_rgb(p, orig_w, orig_h, px0, py0, px1, py1, std::max(1, lwd - 1), r, g, b);
      }
      if (h_len > 0) {
        int cx = (int)std::round(hx[h_len - 1]);
        int cy = (int)std::round(hy[h_len - 1]);
        draw_circle_rgb(p, orig_w, orig_h, cx, cy, 3, r, g, b);
      }
    }

    if (draw_boxes) {
      // 1. Draw bounding box border
      draw_box_rgb(p, orig_w, orig_h, x1, y1, x2, y2, cur_lwd, r, g, b);

      bool draw_label = show_text && (font_scale > 0.05);
      bool use_micro = (font_scale < 0.75);
      int scale_int = std::max(1, (int)std::round(font_scale));

      // 2. Format label string
      std::string lbl = "";
      if (draw_label) {
        if (show_id && has_tracks && i < t_ids.size()) {
          lbl += "#" + std::to_string(t_ids[i]);
        }
        if (show_class && i < labels.size()) {
          std::string cl = Rcpp::as<std::string>(labels[i]);
          if (!cl.empty()) {
            if (!lbl.empty()) lbl += " ";
            lbl += cl;
          }
        }
        if (show_conf && i < scores.size() && !NumericVector::is_na(scores[i])) {
          char buf[32];
          std::snprintf(buf, sizeof(buf), " %.2f", scores[i]);
          lbl += buf;
        }
      }

      if (!lbl.empty()) {
        int char_w = use_micro ? 4 : (8 * scale_int);
        int char_h = use_micro ? 5 : (8 * scale_int);
        int tag_pad = use_micro ? 1 : (2 * scale_int);
        int text_w = (int)lbl.size() * char_w;
        int text_h = char_h;
        int tag_w = text_w + 2 * tag_pad;
        int tag_h = text_h + 2 * tag_pad;

        int tag_y1 = y1 - tag_h;
        int tag_y2 = y1;
        if (tag_y1 < 0) {
          tag_y1 = y1;
          tag_y2 = y1 + tag_h;
        }
        int tag_x1 = x1;
        int tag_x2 = x1 + tag_w;
        if (tag_x2 >= orig_w) {
          tag_x2 = orig_w - 1;
          tag_x1 = std::max(0, tag_x2 - tag_w);
        }

        // Draw tag background
        fill_rect_rgb(p, orig_w, orig_h, tag_x1, tag_y1, tag_x2, tag_y2, r, g, b);

        // Contrast text color
        uint8_t tr = 255, tg = 255, tb = 255;
        double luma = 0.299 * r + 0.587 * g + 0.114 * b;
        if (luma > 180.0) {
          tr = 0; tg = 0; tb = 0;
        }

        // Draw text
        if (use_micro) {
          draw_text_micro_rgb(p, orig_w, orig_h, tag_x1 + tag_pad, tag_y1 + tag_pad, lbl, tr, tg, tb);
        } else {
          draw_text_rgb(p, orig_w, orig_h, tag_x1 + tag_pad, tag_y1 + tag_pad, lbl, scale_int, tr, tg, tb);
        }
      }
    }

    // 3. Draw pose keypoints & skeleton limbs if present
    if (has_kpts_input && draw_skeleton) {
      int limb_lwd = std::max(1, lwd);
      // a) Draw skeleton limbs
      for (int pair_idx = 0; pair_idx < 16; ++pair_idx) {
        int p1 = SKELETON_PAIRS[pair_idx][0];
        int p2 = SKELETON_PAIRS[pair_idx][1];
        double conf1 = kpt_mat(i, p1 * 3 + 2);
        double conf2 = kpt_mat(i, p2 * 3 + 2);
        if (conf1 >= kpt_threshold && conf2 >= kpt_threshold) {
          int px1 = (int)std::round(kpt_mat(i, p1 * 3 + 0));
          int py1 = (int)std::round(kpt_mat(i, p1 * 3 + 1));
          int px2 = (int)std::round(kpt_mat(i, p2 * 3 + 0));
          int py2 = (int)std::round(kpt_mat(i, p2 * 3 + 1));
          uint8_t lr = LIMB_COLORS[pair_idx][0];
          uint8_t lg = LIMB_COLORS[pair_idx][1];
          uint8_t lb = LIMB_COLORS[pair_idx][2];
          draw_line_rgb(p, orig_w, orig_h, px1, py1, px2, py2, limb_lwd, lr, lg, lb);
        }
      }

      // b) Draw keypoint circles
      int rad = std::max(2, kpt_radius);
      for (int k = 0; k < 17; ++k) {
        double conf = kpt_mat(i, k * 3 + 2);
        if (conf >= kpt_threshold) {
          int kx = (int)std::round(kpt_mat(i, k * 3 + 0));
          int ky = (int)std::round(kpt_mat(i, k * 3 + 1));
          draw_circle_rgb(p, orig_w, orig_h, kx, ky, rad + 1, 255, 255, 255);
          draw_circle_rgb(p, orig_w, orig_h, kx, ky, rad - 1, 0, 229, 255);
        }
      }
    }
  }

  // 4. Draw Modern Telemetry HUD badge if requested
  if (!hud_text.empty()) {
    draw_hud_modern(p, orig_w, orig_h, hud_text, hud_pos, hud_layout, font_scale);
  }

  return true;
}

// [[Rcpp::export]]
Rcpp::RawVector crop_bgr_cpp(Rcpp::RawVector bm, int x1, int y1, int x2, int y2) {
  SEXP dim_attr = Rf_getAttrib(bm, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 3) return bm;
  int* dims = INTEGER(dim_attr);
  int orig_w = dims[1];
  int orig_h = dims[2];

  x1 = std::max(0, std::min(orig_w - 1, x1));
  y1 = std::max(0, std::min(orig_h - 1, y1));
  x2 = std::max(x1 + 1, std::min(orig_w, x2));
  y2 = std::max(y1 + 1, std::min(orig_h, y2));

  int cw = x2 - x1;
  int ch = y2 - y1;
  const unsigned char* src = RAW(bm);

  Rcpp::RawVector out_crop(3 * cw * ch);
  unsigned char* dst = RAW(out_crop);

  #pragma omp parallel for if(cw * ch > 20000)
  for (int cy = 0; cy < ch; ++cy) {
    int sy = y1 + cy;
    for (int cx = 0; cx < cw; ++cx) {
      int sx = x1 + cx;
      size_t src_idx = (size_t)sx * 3 + (size_t)sy * ((size_t)3 * orig_w);
      size_t dst_idx = (size_t)cx * 3 + (size_t)cy * ((size_t)3 * cw);
      dst[dst_idx + 0] = src[src_idx + 0];
      dst[dst_idx + 1] = src[src_idx + 1];
      dst[dst_idx + 2] = src[src_idx + 2];
    }
  }

  Rcpp::IntegerVector new_dims(3);
  new_dims[0] = 3;
  new_dims[1] = cw;
  new_dims[2] = ch;
  out_crop.attr("dim") = new_dims;

  return out_crop;
}

// [[Rcpp::export]]
bool save_frame_bgr_cpp(Rcpp::RawVector bm, std::string filename, int quality = 90) {
  SEXP dim_attr = Rf_getAttrib(bm, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 3) return false;
  int* dims = INTEGER(dim_attr);
  int orig_w = dims[1];
  int orig_h = dims[2];
  const unsigned char* p = RAW(bm);

  if (filename.size() >= 4) {
    std::string ext = filename.substr(filename.size() - 4);
    for (char &c : ext) c = (char)::tolower(c);
    if (ext == ".bmp") {
      int ret = stbi_write_bmp(filename.c_str(), orig_w, orig_h, 3, p);
      return ret != 0;
    }
  }

  int ret = stbi_write_jpg(filename.c_str(), orig_w, orig_h, 3, p, quality);
  return ret != 0;
}

// [[Rcpp::export]]
bool swap_rgb_bgr_inplace_cpp(Rcpp::RawVector bm) {
  SEXP dim_attr = Rf_getAttrib(bm, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 3) return false;
  int* dims = INTEGER(dim_attr);
  int orig_w = dims[1];
  int orig_h = dims[2];
  unsigned char* p = RAW(bm);
  size_t total_px = (size_t)orig_w * orig_h;

  #pragma omp parallel for if(total_px > 50000)
  for (size_t i = 0; i < total_px; ++i) {
    size_t idx = i * 3;
    unsigned char tmp = p[idx];
    p[idx] = p[idx + 2];
    p[idx + 2] = tmp;
  }
  return true;
}

// [[Rcpp::export]]
Rcpp::RawVector load_frame_bgr_cpp(std::string filename) {
  int orig_w = 0, orig_h = 0, orig_ch = 0;
  unsigned char* data = stbi_load(filename.c_str(), &orig_w, &orig_h, &orig_ch, 3);
  if (!data) return Rcpp::RawVector(0);

  size_t total_bytes = (size_t)orig_w * orig_h * 3;
  Rcpp::RawVector bm(total_bytes);
  std::memcpy(RAW(bm), data, total_bytes);

  stbi_image_free(data);
  bm.attr("dim") = Rcpp::IntegerVector::create(3, orig_w, orig_h);
  return bm;
}
