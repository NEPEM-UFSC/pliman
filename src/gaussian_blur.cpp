// gaussian_blur.cpp — Highly optimized separable Gaussian blur, brightness, and noise with OpenMP
//
// Optimizations:
//   1. Separable convolution: O(k) per pixel instead of O(k²)
//   2. OpenMP multi-threading across rows, columns, and channels
//   3. Line-buffered cache-contiguous memory access (eliminates strided cache misses)
//   4. Boundary reflection separated from interior loops (eliminates branches in inner loop)
//   5. SIMD vectorization friendly

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

static std::vector<double> make_gaussian_kernel(double sigma) {
  int r = (int)std::ceil(3.0 * sigma);
  if (r < 1) r = 1;
  int ksize = 2 * r + 1;
  std::vector<double> kernel(ksize);

  double sum = 0.0;
  double inv2s2 = 1.0 / (2.0 * sigma * sigma);
  for (int i = 0; i < ksize; ++i) {
    double d = (double)(i - r);
    kernel[i] = std::exp(-d * d * inv2s2);
    sum += kernel[i];
  }
  double inv_sum = 1.0 / sum;
  for (int i = 0; i < ksize; ++i) kernel[i] *= inv_sum;

  return kernel;
}

// Pass 1: horizontal convolution (along cols) for channel plane
// src is nrow x ncol (column-major)
// Writes to tmp of size nrow x ncol
static void blur_horizontal_fast(const double* src, double* tmp,
                                 int nrow, int ncol, int r,
                                 const double* k_data, int ksize) {
  #pragma omp parallel
  {
    std::vector<double> line(ncol + 2 * r);
    #pragma omp for schedule(static)
    for (int row = 0; row < nrow; ++row) {
      // 1. Fill line buffer with boundary reflection
      // Left boundary reflection
      for (int k = 0; k < r; ++k) {
        int jj = r - k;
        if (jj >= ncol) jj = ncol - 1;
        line[k] = src[row + (size_t)jj * nrow];
      }
      // Interior
      for (int col = 0; col < ncol; ++col) {
        line[r + col] = src[row + (size_t)col * nrow];
      }
      // Right boundary reflection
      for (int k = 0; k < r; ++k) {
        int jj = ncol - 2 - k;
        if (jj < 0) jj = 0;
        line[r + ncol + k] = src[row + (size_t)jj * nrow];
      }

      // 2. Convolve line with kernel (100% contiguous, no branches)
      for (int col = 0; col < ncol; ++col) {
        const double* l_ptr = &line[col];
        double acc = 0.0;
        #pragma omp simd
        for (int k = 0; k < ksize; ++k) {
          acc += l_ptr[k] * k_data[k];
        }
        tmp[row + (size_t)col * nrow] = acc;
      }
    }
  }
}

// Pass 2: vertical convolution (along rows) for channel plane
// src (tmp) is nrow x ncol (column-major)
// dst is nrow x ncol
static void blur_vertical_fast(const double* tmp, double* dst,
                               int nrow, int ncol, int r,
                               const double* k_data, int ksize) {
  #pragma omp parallel
  {
    std::vector<double> line(nrow + 2 * r);
    #pragma omp for schedule(static)
    for (int col = 0; col < ncol; ++col) {
      const double* c_src = tmp + (size_t)col * nrow;
      double* c_dst = dst + (size_t)col * nrow;

      // 1. Fill line buffer with boundary reflection
      // Top boundary reflection
      for (int k = 0; k < r; ++k) {
        int ii = r - k;
        if (ii >= nrow) ii = nrow - 1;
        line[k] = c_src[ii];
      }
      // Interior
      for (int row = 0; row < nrow; ++row) {
        line[r + row] = c_src[row];
      }
      // Bottom boundary reflection
      for (int k = 0; k < r; ++k) {
        int ii = nrow - 2 - k;
        if (ii < 0) ii = 0;
        line[r + nrow + k] = c_src[ii];
      }

      // 2. Convolve line with kernel (100% contiguous, no branches)
      for (int row = 0; row < nrow; ++row) {
        const double* l_ptr = &line[row];
        double acc = 0.0;
        #pragma omp simd
        for (int k = 0; k < ksize; ++k) {
          acc += l_ptr[k] * k_data[k];
        }
        c_dst[row] = acc;
      }
    }
  }
}

// [[Rcpp::export]]
SEXP gaussian_blur_cpp(SEXP img_sexp, double sigma) {
  if (sigma <= 0.0) {
    return Rcpp::clone(img_sexp);
  }

  std::vector<double> kernel = make_gaussian_kernel(sigma);
  int ksize = (int)kernel.size();
  int r = ksize / 2;
  const double* k_data = kernel.data();

  Rcpp::NumericVector img_vec(img_sexp);
  Rcpp::IntegerVector dims = img_vec.attr("dim");

  int nrow, ncol, nch;
  if (dims.size() == 2) {
    nrow = dims[0];
    ncol = dims[1];
    nch  = 1;
  } else if (dims.size() == 3) {
    nrow = dims[0];
    ncol = dims[1];
    nch  = dims[2];
  } else {
    Rcpp::stop("gaussian_blur_cpp: input must be a 2-D matrix or 3-D array.");
    return R_NilValue;
  }

  size_t plane_size = (size_t)nrow * ncol;
  Rcpp::NumericVector out = Rcpp::no_init(img_vec.size());

  for (int ch = 0; ch < nch; ++ch) {
    const double* src = &img_vec[ch * plane_size];
    double* dst = &out[ch * plane_size];
    std::vector<double> tmp(plane_size);

    blur_horizontal_fast(src, tmp.data(), nrow, ncol, r, k_data, ksize);
    blur_vertical_fast(tmp.data(), dst, nrow, ncol, r, k_data, ksize);
  }

  out.attr("dim") = Rcpp::clone(dims);
  return out;
}

// [[Rcpp::export]]
SEXP image_brightness_cpp(SEXP img_sexp, double factor, double offset) {
  Rcpp::NumericVector img(img_sexp);
  Rcpp::IntegerVector dims = img.attr("dim");
  size_t n = (size_t)img.size();
  Rcpp::NumericVector out = Rcpp::no_init(n);

  const double* src = img.begin();
  double* dst = out.begin();

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (size_t i = 0; i < n; ++i) {
    double v = src[i] * factor + offset;
    if (v < 0.0) v = 0.0;
    else if (v > 1.0) v = 1.0;
    dst[i] = v;
  }

  out.attr("dim") = Rcpp::clone(dims);
  return out;
}

// [[Rcpp::export]]
SEXP image_noise_cpp(SEXP img_sexp, int type, double mean, double sd, double p, int seed) {
  Rcpp::NumericVector img(img_sexp);
  Rcpp::IntegerVector dims = img.attr("dim");
  size_t n = (size_t)img.size();
  Rcpp::NumericVector out = Rcpp::no_init(n);

  const double* src = img.begin();
  double* dst = out.begin();

  // type: 1 = Gaussian, 2 = Salt & Pepper, 3 = Speckle
  if (type == 1) { // Gaussian
    #ifdef _OPENMP
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::mt19937 gen((seed == 0) ? 12345 + tid * 997 : seed + tid * 997);
      std::normal_distribution<double> dist(mean, sd);
      #pragma omp for schedule(static)
      for (size_t i = 0; i < n; ++i) {
        double v = src[i] + dist(gen);
        if (v < 0.0) v = 0.0;
        else if (v > 1.0) v = 1.0;
        dst[i] = v;
      }
    }
    #else
    std::mt19937 gen(seed == 0 ? 12345 : seed);
    std::normal_distribution<double> dist(mean, sd);
    for (size_t i = 0; i < n; ++i) {
      double v = src[i] + dist(gen);
      if (v < 0.0) v = 0.0;
      else if (v > 1.0) v = 1.0;
      dst[i] = v;
    }
    #endif
  } else if (type == 2) { // Salt & Pepper
    #ifdef _OPENMP
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::mt19937 gen((seed == 0) ? 12345 + tid * 997 : seed + tid * 997);
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      double half_p = p * 0.5;
      #pragma omp for schedule(static)
      for (size_t i = 0; i < n; ++i) {
        double r = dist(gen);
        if (r < half_p) {
          dst[i] = 0.0; // pepper
        } else if (r < p) {
          dst[i] = 1.0; // salt
        } else {
          dst[i] = src[i];
        }
      }
    }
    #else
    std::mt19937 gen(seed == 0 ? 12345 : seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double half_p = p * 0.5;
    for (size_t i = 0; i < n; ++i) {
      double r = dist(gen);
      if (r < half_p) {
        dst[i] = 0.0;
      } else if (r < p) {
        dst[i] = 1.0;
      } else {
        dst[i] = src[i];
      }
    }
    #endif
  } else { // Speckle (multiplicative Gaussian: I = I + I * N)
    #ifdef _OPENMP
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::mt19937 gen((seed == 0) ? 12345 + tid * 997 : seed + tid * 997);
      std::normal_distribution<double> dist(mean, sd);
      #pragma omp for schedule(static)
      for (size_t i = 0; i < n; ++i) {
        double v = src[i] + src[i] * dist(gen);
        if (v < 0.0) v = 0.0;
        else if (v > 1.0) v = 1.0;
        dst[i] = v;
      }
    }
    #else
    std::mt19937 gen(seed == 0 ? 12345 : seed);
    std::normal_distribution<double> dist(mean, sd);
    for (size_t i = 0; i < n; ++i) {
      double v = src[i] + src[i] * dist(gen);
      if (v < 0.0) v = 0.0;
      else if (v > 1.0) v = 1.0;
      dst[i] = v;
    }
    #endif
  }

  out.attr("dim") = Rcpp::clone(dims);
  return out;
}
