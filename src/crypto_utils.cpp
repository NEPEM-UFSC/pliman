#include <Rcpp.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cctype>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
static std::wstring utf8_to_wstring(const std::string& str) {
    if (str.empty()) return std::wstring();
    int size_needed = MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), NULL, 0);
    std::wstring wstrTo(size_needed, 0);
    MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), &wstrTo[0], size_needed);
    return wstrTo;
}
#endif

// =============================================================================
// SHA-256 (FIPS 180-4) Implementation - Zero External Dependencies
// =============================================================================

static const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static inline uint32_t rotr(uint32_t x, uint32_t n) {
    return (x >> n) | (x << (32 - n));
}

static inline uint32_t choose(uint32_t e, uint32_t f, uint32_t g) {
    return (e & f) ^ (~e & g);
}

static inline uint32_t majority(uint32_t a, uint32_t b, uint32_t c) {
    return (a & b) ^ (a & c) ^ (b & c);
}

static inline uint32_t sig0(uint32_t x) {
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);
}

static inline uint32_t sig1(uint32_t x) {
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10);
}

static inline uint32_t ep0(uint32_t x) {
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22);
}

static inline uint32_t ep1(uint32_t x) {
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25);
}

struct SHA256_CTX {
    uint8_t data[64];
    uint32_t datalen;
    uint64_t bitlen;
    uint32_t state[8];

    void init() {
        datalen = 0;
        bitlen = 0;
        state[0] = 0x6a09e667;
        state[1] = 0xbb67ae85;
        state[2] = 0x3c6ef372;
        state[3] = 0xa54ff53a;
        state[4] = 0x510e527f;
        state[5] = 0x9b05688c;
        state[6] = 0x1f83d9ab;
        state[7] = 0x5be0cd19;
    }

    void transform() {
        uint32_t m[64];
        for (int i = 0, j = 0; i < 16; ++i, j += 4) {
            m[i] = ((uint32_t)data[j] << 24) |
                   ((uint32_t)data[j + 1] << 16) |
                   ((uint32_t)data[j + 2] << 8) |
                   ((uint32_t)data[j + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            m[i] = sig1(m[i - 2]) + m[i - 7] + sig0(m[i - 15]) + m[i - 16];
        }

        uint32_t a = state[0];
        uint32_t b = state[1];
        uint32_t c = state[2];
        uint32_t d = state[3];
        uint32_t e = state[4];
        uint32_t f = state[5];
        uint32_t g = state[6];
        uint32_t h = state[7];

        for (int i = 0; i < 64; ++i) {
            uint32_t t1 = h + ep1(e) + choose(e, f, g) + K[i] + m[i];
            uint32_t t2 = ep0(a) + majority(a, b, c);
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }

        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
        state[4] += e;
        state[5] += f;
        state[6] += g;
        state[7] += h;
    }

    void update(const uint8_t* in, size_t len) {
        for (size_t i = 0; i < len; ++i) {
            data[datalen] = in[i];
            datalen++;
            if (datalen == 64) {
                transform();
                bitlen += 512;
                datalen = 0;
            }
        }
    }

    void final(uint8_t hash[32]) {
        uint32_t i = datalen;

        if (datalen < 56) {
            data[i++] = 0x80;
            while (i < 56) {
                data[i++] = 0x00;
            }
        } else {
            data[i++] = 0x80;
            while (i < 64) {
                data[i++] = 0x00;
            }
            transform();
            std::memset(data, 0, 56);
        }

        bitlen += (uint64_t)datalen * 8;
        data[63] = (uint8_t)(bitlen);
        data[62] = (uint8_t)(bitlen >> 8);
        data[61] = (uint8_t)(bitlen >> 16);
        data[60] = (uint8_t)(bitlen >> 24);
        data[59] = (uint8_t)(bitlen >> 32);
        data[58] = (uint8_t)(bitlen >> 40);
        data[57] = (uint8_t)(bitlen >> 48);
        data[56] = (uint8_t)(bitlen >> 56);
        transform();

        for (i = 0; i < 4; ++i) {
            hash[i]      = (uint8_t)((state[0] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 4]  = (uint8_t)((state[1] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 8]  = (uint8_t)((state[2] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 12] = (uint8_t)((state[3] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 16] = (uint8_t)((state[4] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 20] = (uint8_t)((state[5] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 24] = (uint8_t)((state[6] >> (24 - i * 8)) & 0x000000ff);
            hash[i + 28] = (uint8_t)((state[7] >> (24 - i * 8)) & 0x000000ff);
        }
    }
};

static std::string to_hex_string(const uint8_t* hash, size_t len = 32) {
    static const char hex_digits[] = "0123456789abcdef";
    std::string out;
    out.reserve(len * 2);
    for (size_t i = 0; i < len; ++i) {
        out.push_back(hex_digits[(hash[i] >> 4) & 0x0F]);
        out.push_back(hex_digits[hash[i] & 0x0F]);
    }
    return out;
}

// =============================================================================
// Base64 (RFC 4648) Implementation - Zero External Dependencies
// =============================================================================

static const char b64_chars[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";

static std::string base64_encode_bytes(const uint8_t* bytes_to_encode, size_t in_len) {
    std::string ret;
    ret.reserve(((in_len + 2) / 3) * 4);
    int i = 0;
    int j = 0;
    uint8_t char_array_3[3];
    uint8_t char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for (i = 0; i < 4; i++) {
                ret.push_back(b64_chars[char_array_4[i]]);
            }
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 3; j++) {
            char_array_3[j] = '\0';
        }

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);

        for (j = 0; j < i + 1; j++) {
            ret.push_back(b64_chars[char_array_4[j]]);
        }

        while (i++ < 3) {
            ret.push_back('=');
        }
    }

    return ret;
}

static inline bool is_base64(uint8_t c) {
    return (std::isalnum(c) || (c == '+') || (c == '/'));
}

static std::vector<uint8_t> base64_decode_string(const std::string& encoded_string) {
    size_t in_len = encoded_string.size();
    int i = 0;
    int j = 0;
    size_t in_ = 0;
    uint8_t char_array_4[4], char_array_3[3];
    std::vector<uint8_t> ret;
    ret.reserve((in_len * 3) / 4);

    while (in_len-- && (encoded_string[in_] != '=') && (is_base64((uint8_t)encoded_string[in_]) || std::isspace((uint8_t)encoded_string[in_]))) {
        if (std::isspace((uint8_t)encoded_string[in_])) {
            in_++;
            continue;
        }
        char_array_4[i++] = (uint8_t)encoded_string[in_]; in_++;
        if (i == 4) {
            for (i = 0; i < 4; i++) {
                const char* p = std::strchr(b64_chars, char_array_4[i]);
                char_array_4[i] = p ? (uint8_t)(p - b64_chars) : 0;
            }

            char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (i = 0; i < 3; i++) {
                ret.push_back(char_array_3[i]);
            }
            i = 0;
        }
    }

    if (i) {
        for (j = 0; j < i; j++) {
            const char* p = std::strchr(b64_chars, char_array_4[j]);
            char_array_4[j] = p ? (uint8_t)(p - b64_chars) : 0;
        }

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);

        for (j = 0; j < i - 1; j++) {
            ret.push_back(char_array_3[j]);
        }
    }

    return ret;
}

// =============================================================================
// Rcpp Exported Functions
// =============================================================================

//' Compute SHA-256 Checksum of a File
//'
//' Streams the file in 64 KB binary chunks and computes the cryptographic SHA-256
//' checksum with zero external dependencies.
//'
//' @param file_path Character string. Path to the file.
//' @return Character string containing the 64-character lowercase hexadecimal hash.
//' @export
// [[Rcpp::export]]
std::string rcpp_sha256_file(std::string file_path) {
#ifdef _WIN32
    FILE* fp = _wfopen(utf8_to_wstring(file_path).c_str(), L"rb");
#else
    FILE* fp = std::fopen(file_path.c_str(), "rb");
#endif
    if (!fp) {
        Rcpp::stop("Cannot open file for SHA-256 calculation: " + file_path);
    }
    SHA256_CTX ctx;
    ctx.init();
    std::vector<uint8_t> buffer(65536);
    size_t bytes_read = 0;
    while ((bytes_read = std::fread(buffer.data(), 1, buffer.size(), fp)) > 0) {
        ctx.update(buffer.data(), bytes_read);
    }
    std::fclose(fp);
    uint8_t hash[32];
    ctx.final(hash);
    return to_hex_string(hash, 32);
}

//' Compute SHA-256 Checksum of Raw Bytes or String
//'
//' @param bytes Raw vector of bytes to hash.
//' @return Character string containing the 64-character lowercase hexadecimal hash.
//' @export
// [[Rcpp::export]]
std::string rcpp_sha256(Rcpp::RawVector bytes) {
    SHA256_CTX ctx;
    ctx.init();
    if (bytes.size() > 0) {
        ctx.update(reinterpret_cast<const uint8_t*>(&bytes[0]), bytes.size());
    }
    uint8_t hash[32];
    ctx.final(hash);
    return to_hex_string(hash, 32);
}

//' Encode Raw Bytes to Base64
//'
//' @param bytes Raw vector of bytes.
//' @return Base64-encoded character string.
//' @export
// [[Rcpp::export]]
std::string rcpp_base64_encode(Rcpp::RawVector bytes) {
    if (bytes.size() == 0) return "";
    return base64_encode_bytes(reinterpret_cast<const uint8_t*>(&bytes[0]), bytes.size());
}

//' Encode File Directly to Base64
//'
//' Reads a binary or text file and returns its Base64-encoded string directly,
//' without requiring external packages like \code{base64enc}.
//'
//' @param file_path Character string. Path to the file.
//' @return Base64-encoded character string.
//' @export
// [[Rcpp::export]]
std::string rcpp_base64_encode_file(std::string file_path) {
#ifdef _WIN32
    FILE* fp = _wfopen(utf8_to_wstring(file_path).c_str(), L"rb");
#else
    FILE* fp = std::fopen(file_path.c_str(), "rb");
#endif
    if (!fp) {
        Rcpp::stop("Cannot open file for Base64 encoding: " + file_path);
    }
    std::fseek(fp, 0, SEEK_END);
    long size = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);

    if (size <= 0) {
        std::fclose(fp);
        return "";
    }

    std::vector<uint8_t> buffer((size_t)size);
    size_t read_bytes = std::fread(buffer.data(), 1, (size_t)size, fp);
    std::fclose(fp);
    return base64_encode_bytes(buffer.data(), read_bytes);
}

//' Decode Base64 String to Raw Bytes
//'
//' @param text Base64-encoded character string.
//' @return Raw vector of decoded bytes.
//' @export
// [[Rcpp::export]]
Rcpp::RawVector rcpp_base64_decode(std::string text) {
    std::vector<uint8_t> dec = base64_decode_string(text);
    Rcpp::RawVector out(dec.size());
    if (!dec.empty()) {
        std::memcpy(&out[0], dec.data(), dec.size());
    }
    return out;
}
