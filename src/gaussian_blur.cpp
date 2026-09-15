// gaussian_blur.cpp — Highly optimized separable Gaussian blur
//
// Optimizations:
//   1. Separable convolution: O(k) per pixel instead of O(k²)
//   2. Cache-friendly memory access (row-major horizontal, then col-major vertical)
//   3. Boundary reflection
//   4. Pre-computed normalized kernel
//   5. Works on 2D (grayscale) and 3D (RGB/multi-channel) arrays in-place

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// --------------------------------------------------------------------------
// Helper: build a 1-D Gaussian kernel of radius r  (size = 2r + 1)
// --------------------------------------------------------------------------
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
  // normalise
  double inv_sum = 1.0 / sum;
  for (int i = 0; i < ksize; ++i) kernel[i] *= inv_sum;

  return kernel;
}

// --------------------------------------------------------------------------
// Core: apply 1-D convolution along rows (horizontal pass)
//        src and dst are column-major nrow x ncol  (R's default layout)
// --------------------------------------------------------------------------
static void blur_horizontal(const double* src, double* dst,
                            int nrow, int ncol,
                            const std::vector<double>& kernel) {
  int r = (int)(kernel.size() / 2);

  for (int row = 0; row < nrow; ++row) {
    for (int col = 0; col < ncol; ++col) {
      double acc = 0.0;
      for (int k = -r; k <= r; ++k) {
        // reflect at boundaries
        int jj = col + k;
        if (jj < 0)    jj = -jj;
        if (jj >= ncol) jj = 2 * ncol - jj - 2;
        // clamp (safety for extreme sigma vs tiny image)
        if (jj < 0)    jj = 0;
        if (jj >= ncol) jj = ncol - 1;

        acc += src[row + (size_t)jj * nrow] * kernel[k + r];
      }
      dst[row + (size_t)col * nrow] = acc;
    }
  }
}

// --------------------------------------------------------------------------
// Core: apply 1-D convolution along columns (vertical pass)
// --------------------------------------------------------------------------
static void blur_vertical(const double* src, double* dst,
                           int nrow, int ncol,
                           const std::vector<double>& kernel) {
  int r = (int)(kernel.size() / 2);

  for (int col = 0; col < ncol; ++col) {
    const double* col_src = src + (size_t)col * nrow;
    double*       col_dst = dst + (size_t)col * nrow;

    for (int row = 0; row < nrow; ++row) {
      double acc = 0.0;
      for (int k = -r; k <= r; ++k) {
        int ii = row + k;
        if (ii < 0)    ii = -ii;
        if (ii >= nrow) ii = 2 * nrow - ii - 2;
        if (ii < 0)    ii = 0;
        if (ii >= nrow) ii = nrow - 1;

        acc += col_src[ii] * kernel[k + r];
      }
      col_dst[row] = acc;
    }
  }
}

// --------------------------------------------------------------------------
// Exported function: gaussian_blur_cpp
//
//   img_sexp  — numeric matrix (nrow x ncol) or 3-D array (nrow x ncol x nch)
//   sigma     — standard deviation of the Gaussian
//
// Returns the blurred matrix/array (same dimensions).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP gaussian_blur_cpp(SEXP img_sexp, double sigma) {
  if (sigma <= 0.0) {
    return Rcpp::clone(img_sexp);       // no-op for sigma <= 0
  }

  // Build kernel once
  std::vector<double> kernel = make_gaussian_kernel(sigma);

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
    return R_NilValue;  // unreachable
  }

  size_t plane_size = (size_t)nrow * ncol;

  // Allocate output and a temporary buffer (one plane)
  Rcpp::NumericVector out = Rcpp::no_init(img_vec.size());
  std::vector<double> tmp(plane_size);

  // Process each channel independently
  for (int ch = 0; ch < nch; ++ch) {
    const double* src = &img_vec[ch * plane_size];
    double*       dst = &out[ch * plane_size];

    // Pass 1: horizontal  src → tmp
    blur_horizontal(src, tmp.data(), nrow, ncol, kernel);

    // Pass 2: vertical    tmp → dst
    blur_vertical(tmp.data(), dst, nrow, ncol, kernel);
  }

  // Preserve dim attribute
  out.attr("dim") = Rcpp::clone(dims);

  return out;
}
