// image_alignment.cpp — Native C++ Sub-Pixel Image Alignment & Registration for pliman
// Multi-scale Gaussian Pyramid Enhanced Correlation Coefficient (ECC) Maximization
// Zero Python, Zero OpenCV, Zero External Dependencies.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Bilinear sampling helper for a single channel (double)
static inline double sample_bilinear_clamp(const double* img, int w, int h, double x, double y) {
  if (x <= 0.0) x = 0.0;
  else if (x >= w - 1.0) x = w - 1.0;
  if (y <= 0.0) y = 0.0;
  else if (y >= h - 1.0) y = h - 1.0;
  
  int x0 = (int)x;
  int y0 = (int)y;
  int x1 = (x0 + 1 < w) ? x0 + 1 : x0;
  int y1 = (y0 + 1 < h) ? y0 + 1 : y0;
  
  double fx = x - (double)x0;
  double fy = y - (double)y0;
  
  double v00 = img[x0 + y0 * w];
  double v10 = img[x1 + y0 * w];
  double v01 = img[x0 + y1 * w];
  double v11 = img[x1 + y1 * w];
  
  return (1.0 - fx) * ((1.0 - fy) * v00 + fy * v01) +
         fx * ((1.0 - fy) * v10 + fy * v11);
}

// 5x5 separable Gaussian smoothing [1, 4, 6, 4, 1] / 16
static void gaussian_blur_5x5(const double* src, double* dst, int w, int h) {
  std::vector<double> tmp(w * h);
  
  // Horizontal pass
  for (int y = 0; y < h; ++y) {
    const double* src_row = src + y * w;
    double* tmp_row = tmp.data() + y * w;
    for (int x = 0; x < w; ++x) {
      int xm2 = std::max(0, x - 2);
      int xm1 = std::max(0, x - 1);
      int xp1 = std::min(w - 1, x + 1);
      int xp2 = std::min(w - 1, x + 2);
      tmp_row[x] = (src_row[xm2] + 4.0 * src_row[xm1] + 6.0 * src_row[x] + 4.0 * src_row[xp1] + src_row[xp2]) * (1.0 / 16.0);
    }
  }
  
  // Vertical pass
  for (int x = 0; x < w; ++x) {
    for (int y = 0; y < h; ++y) {
      int ym2 = std::max(0, y - 2);
      int ym1 = std::max(0, y - 1);
      int yp1 = std::min(h - 1, y + 1);
      int yp2 = std::min(h - 1, y + 2);
      dst[x + y * w] = (tmp[x + ym2 * w] + 4.0 * tmp[x + ym1 * w] + 6.0 * tmp[x + y * w] + 4.0 * tmp[x + yp1 * w] + tmp[x + yp2 * w]) * (1.0 / 16.0);
    }
  }
}

// Downsample 2D image by 2x using 2x2 box average
static void downsample_2x(const double* src, int w, int h, double* dst) {
  int w2 = w / 2;
  int h2 = h / 2;
  for (int y = 0; y < h2; ++y) {
    for (int x = 0; x < w2; ++x) {
      int sx = 2 * x;
      int sy = 2 * y;
      dst[x + y * w2] = 0.25 * (src[sx + sy * w] + src[sx + 1 + sy * w] +
                               src[sx + (sy + 1) * w] + src[sx + 1 + (sy + 1) * w]);
    }
  }
}

// Central difference gradients
static void compute_gradients(const double* src, double* gx, double* gy, int w, int h) {
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int xm1 = std::max(0, x - 1);
      int xp1 = std::min(w - 1, x + 1);
      int ym1 = std::max(0, y - 1);
      int yp1 = std::min(h - 1, y + 1);
      gx[x + y * w] = 0.5 * (src[xp1 + y * w] - src[xm1 + y * w]);
      gy[x + y * w] = 0.5 * (src[x + yp1 * w] - src[x + ym1 * w]);
    }
  }
}

// Structure for Pyramid level
struct PyramidLevel {
  int w;
  int h;
  std::vector<double> img;
  std::vector<double> smooth;
  std::vector<double> gx;
  std::vector<double> gy;
};

// Build image pyramid
static std::vector<PyramidLevel> build_pyramid(const double* full_img, int w, int h, int num_levels) {
  std::vector<PyramidLevel> pyr(num_levels);
  
  pyr[0].w = w;
  pyr[0].h = h;
  pyr[0].img.assign(full_img, full_img + w * h);
  pyr[0].smooth.resize(w * h);
  pyr[0].gx.resize(w * h);
  pyr[0].gy.resize(w * h);
  gaussian_blur_5x5(pyr[0].img.data(), pyr[0].smooth.data(), w, h);
  compute_gradients(pyr[0].smooth.data(), pyr[0].gx.data(), pyr[0].gy.data(), w, h);
  
  for (int l = 1; l < num_levels; ++l) {
    pyr[l].w = pyr[l - 1].w / 2;
    pyr[l].h = pyr[l - 1].h / 2;
    int pw = pyr[l].w;
    int ph = pyr[l].h;
    pyr[l].img.resize(pw * ph);
    pyr[l].smooth.resize(pw * ph);
    pyr[l].gx.resize(pw * ph);
    pyr[l].gy.resize(pw * ph);
    
    downsample_2x(pyr[l - 1].smooth.data(), pyr[l - 1].w, pyr[l - 1].h, pyr[l].img.data());
    gaussian_blur_5x5(pyr[l].img.data(), pyr[l].smooth.data(), pw, ph);
    compute_gradients(pyr[l].smooth.data(), pyr[l].gx.data(), pyr[l].gy.data(), pw, ph);
  }
  
  return pyr;
}

// Single-scale ECC optimization
static bool ecc_single_scale(const PyramidLevel& ref_lvl,
                             const PyramidLevel& src_lvl,
                             arma::mat& W,
                             bool is_affine,
                             int max_iter,
                             double eps) {
  const int w = ref_lvl.w;
  const int h = ref_lvl.h;
  const int num_params = is_affine ? 6 : 2;
  
  const double* T_smooth = ref_lvl.smooth.data();
  const double* I_smooth = src_lvl.smooth.data();
  const double* Gx_src = src_lvl.gx.data();
  const double* Gy_src = src_lvl.gy.data();
  
  std::vector<double> Iw(w * h);
  std::vector<double> Gxw(w * h);
  std::vector<double> Gyw(w * h);
  std::vector<uint8_t> mask(w * h);
  
  double rho = -1.0;
  double last_rho = -1.0;
  
  for (int it = 0; it < max_iter; ++it) {
    double m00 = W(0, 0);
    double m01 = W(0, 1);
    double m02 = W(0, 2);
    double m10 = W(1, 0);
    double m11 = W(1, 1);
    double m12 = W(1, 2);
    
    int n_valid = 0;
    double sum_T = 0.0;
    double sum_Iw = 0.0;
    
    // 1. Warp image & gradients and compute means over valid overlap
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        double xw = m00 * x + m01 * y + m02;
        double yw = m10 * x + m11 * y + m12;
        
        int idx = x + y * w;
        if (xw >= 0.0 && xw <= w - 1.0 && yw >= 0.0 && yw <= h - 1.0) {
          mask[idx] = 1;
          double iw_val = sample_bilinear_clamp(I_smooth, w, h, xw, yw);
          double gxw_val = sample_bilinear_clamp(Gx_src, w, h, xw, yw);
          double gyw_val = sample_bilinear_clamp(Gy_src, w, h, xw, yw);
          Iw[idx] = iw_val;
          Gxw[idx] = gxw_val;
          Gyw[idx] = gyw_val;
          
          sum_T += T_smooth[idx];
          sum_Iw += iw_val;
          n_valid++;
        } else {
          mask[idx] = 0;
        }
      }
    }
    
    if (n_valid < 0.25 * w * h) {
      return false; // Insufficient overlap
    }
    
    double mean_T = sum_T / n_valid;
    double mean_Iw = sum_Iw / n_valid;
    
    double norm_T_sq = 0.0;
    double norm_Iw_sq = 0.0;
    double correlation = 0.0;
    
    // 2. Accumulate Hessian and projections
    arma::mat H(num_params, num_params, arma::fill::zeros);
    arma::vec v_i(num_params, arma::fill::zeros);
    arma::vec v_t(num_params, arma::fill::zeros);
    
    double j[6];
    
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        int idx = x + y * w;
        if (!mask[idx]) continue;
        
        double t_bar = T_smooth[idx] - mean_T;
        double iw_bar = Iw[idx] - mean_Iw;
        
        norm_T_sq += t_bar * t_bar;
        norm_Iw_sq += iw_bar * iw_bar;
        correlation += t_bar * iw_bar;
        
        double gx = Gxw[idx];
        double gy = Gyw[idx];
        
        if (is_affine) {
          j[0] = x * gx;
          j[1] = x * gy;
          j[2] = y * gx;
          j[3] = y * gy;
          j[4] = gx;
          j[5] = gy;
        } else {
          j[0] = gx;
          j[1] = gy;
        }
        
        for (int p = 0; p < num_params; ++p) {
          v_i(p) += j[p] * iw_bar;
          v_t(p) += j[p] * t_bar;
          for (int q = p; q < num_params; ++q) {
            H(p, q) += j[p] * j[q];
          }
        }
      }
    }
    
    // Symmetrize H
    for (int p = 0; p < num_params; ++p) {
      for (int q = 0; q < p; ++q) {
        H(p, q) = H(q, p);
      }
    }
    
    double norm_T = std::sqrt(norm_T_sq);
    double norm_Iw = std::sqrt(norm_Iw_sq);
    if (norm_T < 1e-10 || norm_Iw < 1e-10) return false;
    
    last_rho = rho;
    rho = correlation / (norm_T * norm_Iw);
    
    if (it > 0 && std::fabs(rho - last_rho) < eps) {
      break; // Converged
    }
    
    // Regularize H
    for (int p = 0; p < num_params; ++p) {
      H(p, p) += 1e-6 * H(p, p) + 1e-8;
    }
    
    arma::mat H_inv;
    if (!arma::inv(H_inv, H)) {
      if (!arma::pinv(H_inv, H)) return false;
    }
    
    arma::vec q_proj = H_inv * v_i;
    arma::vec p_proj = H_inv * v_t;
    
    double lambda_n = norm_Iw_sq - arma::dot(v_i, q_proj);
    double lambda_d = correlation - arma::dot(v_t, q_proj);
    
    if (lambda_d <= 0.0) {
      return false; // Correlation would decrease
    }
    
    double lambda = lambda_n / lambda_d;
    arma::vec delta_p = lambda * p_proj - q_proj;
    
    if (is_affine) {
      W(0, 0) += delta_p(0);
      W(1, 0) += delta_p(1);
      W(0, 1) += delta_p(2);
      W(1, 1) += delta_p(3);
      W(0, 2) += delta_p(4);
      W(1, 2) += delta_p(5);
    } else {
      W(0, 2) += delta_p(0);
      W(1, 2) += delta_p(1);
    }
    
    if (arma::norm(delta_p, 2) < eps) {
      break; // Step converged
    }
  }
  
  return true;
}

// Convert arbitrary R image / array to double grayscale vector [w * h]
static std::vector<double> extract_grayscale(SEXP img_sexp, int& out_w, int& out_h, int& out_c) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);
  if (dims.size() == 2) {
    out_w = dims[0];
    out_h = dims[1];
    out_c = 1;
  } else if (dims.size() >= 3) {
    out_w = dims[0];
    out_h = dims[1];
    out_c = dims[2];
  } else {
    stop("extract_grayscale: Input must have 2 or 3 dimensions.");
  }
  
  int total_pixels = out_w * out_h;
  std::vector<double> gray(total_pixels);
  
  if (TYPEOF(img_sexp) == RAWSXP) {
    Rcpp::RawVector raw_v(img_sexp);
    const Rbyte* p = raw_v.begin();
    if (out_c >= 3) {
      const Rbyte* pr = p;
      const Rbyte* pg = p + total_pixels;
      const Rbyte* pb = p + 2 * total_pixels;
      for (int i = 0; i < total_pixels; ++i) {
        gray[i] = (0.299 * pr[i] + 0.587 * pg[i] + 0.114 * pb[i]) / 255.0;
      }
    } else {
      for (int i = 0; i < total_pixels; ++i) {
        gray[i] = p[i] / 255.0;
      }
    }
  } else if (TYPEOF(img_sexp) == REALSXP) {
    Rcpp::NumericVector num_v(img_sexp);
    const double* p = num_v.begin();
    if (out_c >= 3) {
      const double* pr = p;
      const double* pg = p + total_pixels;
      const double* pb = p + 2 * total_pixels;
      for (int i = 0; i < total_pixels; ++i) {
        gray[i] = 0.299 * pr[i] + 0.587 * pg[i] + 0.114 * pb[i];
      }
    } else {
      for (int i = 0; i < total_pixels; ++i) {
        gray[i] = p[i];
      }
    }
  } else {
    stop("extract_grayscale: Unsupported image data type (must be raw or numeric).");
  }
  
  return gray;
}

// Invert 2x3 affine matrix
static bool invert_affine_2x3(const arma::mat& W, arma::mat& W_inv) {
  double a00 = W(0, 0);
  double a01 = W(0, 1);
  double a02 = W(0, 2);
  double a10 = W(1, 0);
  double a11 = W(1, 1);
  double a12 = W(1, 2);
  
  double det = a00 * a11 - a01 * a10;
  if (std::fabs(det) < 1e-12) {
    W_inv = arma::eye(2, 3);
    return false;
  }
  
  double inv_det = 1.0 / det;
  double inv_a00 = a11 * inv_det;
  double inv_a01 = -a01 * inv_det;
  double inv_a10 = -a10 * inv_det;
  double inv_a11 = a00 * inv_det;
  
  double inv_t0 = -(inv_a00 * a02 + inv_a01 * a12);
  double inv_t1 = -(inv_a10 * a02 + inv_a11 * a12);
  
  W_inv.set_size(2, 3);
  W_inv(0, 0) = inv_a00;
  W_inv(0, 1) = inv_a01;
  W_inv(0, 2) = inv_t0;
  W_inv(1, 0) = inv_a10;
  W_inv(1, 1) = inv_a11;
  W_inv(1, 2) = inv_t1;
  
  return true;
}

// Warp and crop a single image
static SEXP warp_and_crop_image(SEXP img_sexp,
                                const arma::mat& W,
                                int crop_x1, int crop_y1,
                                int out_w, int out_h) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);
  int in_w = dims[0];
  int in_h = dims[1];
  int nch = (dims.size() >= 3) ? dims[2] : 1;
  
  double m00 = W(0, 0);
  double m01 = W(0, 1);
  double m02 = W(0, 2);
  double m10 = W(1, 0);
  double m11 = W(1, 1);
  double m12 = W(1, 2);
  
  std::ptrdiff_t in_plane = (std::ptrdiff_t)in_w * in_h;
  std::ptrdiff_t out_plane = (std::ptrdiff_t)out_w * out_h;
  
  if (TYPEOF(img_sexp) == RAWSXP) {
    Rcpp::RawVector in_vec(img_sexp);
    const Rbyte* src = in_vec.begin();
    
    Rcpp::RawVector out_vec = Rcpp::no_init((std::size_t)out_w * out_h * nch);
    Rbyte* dst = out_vec.begin();
    
    #pragma omp parallel for collapse(2) schedule(static)
    for (int c = 0; c < nch; ++c) {
      for (int yo = 0; yo < out_h; ++yo) {
        int y_ref = crop_y1 + yo;
        const Rbyte* src_ch = src + c * in_plane;
        Rbyte* dst_row = dst + c * out_plane + yo * out_w;
        
        for (int xo = 0; xo < out_w; ++xo) {
          int x_ref = crop_x1 + xo;
          double xw = m00 * x_ref + m01 * y_ref + m02;
          double yw = m10 * x_ref + m11 * y_ref + m12;
          
          if (xw <= 0.0) xw = 0.0;
          else if (xw >= in_w - 1.0) xw = in_w - 1.0;
          if (yw <= 0.0) yw = 0.0;
          else if (yw >= in_h - 1.0) yw = in_h - 1.0;
          
          int x0 = (int)xw;
          int y0 = (int)yw;
          int x1 = (x0 + 1 < in_w) ? x0 + 1 : x0;
          int y1 = (y0 + 1 < in_h) ? y0 + 1 : y0;
          
          double fx = xw - (double)x0;
          double fy = yw - (double)y0;
          
          double v00 = src_ch[x0 + y0 * in_w];
          double v10 = src_ch[x1 + y0 * in_w];
          double v01 = src_ch[x0 + y1 * in_w];
          double v11 = src_ch[x1 + y1 * in_w];
          
          double val = (1.0 - fx) * ((1.0 - fy) * v00 + fy * v01) +
                        fx * ((1.0 - fy) * v10 + fy * v11);
          
          int ival = (int)std::round(val);
          if (ival < 0) ival = 0;
          else if (ival > 255) ival = 255;
          
          dst_row[xo] = (Rbyte)ival;
        }
      }
    }
    
    if (nch == 1) {
      out_vec.attr("dim") = IntegerVector::create(out_w, out_h);
      out_vec.attr("colormode") = "Grayscale";
    } else {
      out_vec.attr("dim") = IntegerVector::create(out_w, out_h, nch);
      out_vec.attr("colormode") = "Color";
    }
    out_vec.attr("class") = CharacterVector::create("image", "Image", "array");
    out_vec.attr("colspace") = "RGB";
    out_vec.attr("gamma") = 1.0;
    return out_vec;
  } else {
    Rcpp::NumericVector in_vec(img_sexp);
    const double* src = in_vec.begin();
    
    Rcpp::NumericVector out_vec = Rcpp::no_init((std::size_t)out_w * out_h * nch);
    double* dst = out_vec.begin();
    
    #pragma omp parallel for collapse(2) schedule(static)
    for (int c = 0; c < nch; ++c) {
      for (int yo = 0; yo < out_h; ++yo) {
        int y_ref = crop_y1 + yo;
        const double* src_ch = src + c * in_plane;
        double* dst_row = dst + c * out_plane + yo * out_w;
        
        for (int xo = 0; xo < out_w; ++xo) {
          int x_ref = crop_x1 + xo;
          double xw = m00 * x_ref + m01 * y_ref + m02;
          double yw = m10 * x_ref + m11 * y_ref + m12;
          
          if (xw <= 0.0) xw = 0.0;
          else if (xw >= in_w - 1.0) xw = in_w - 1.0;
          if (yw <= 0.0) yw = 0.0;
          else if (yw >= in_h - 1.0) yw = in_h - 1.0;
          
          int x0 = (int)xw;
          int y0 = (int)yw;
          int x1 = (x0 + 1 < in_w) ? x0 + 1 : x0;
          int y1 = (y0 + 1 < in_h) ? y0 + 1 : y0;
          
          double fx = xw - (double)x0;
          double fy = yw - (double)y0;
          
          double v00 = src_ch[x0 + y0 * in_w];
          double v10 = src_ch[x1 + y0 * in_w];
          double v01 = src_ch[x0 + y1 * in_w];
          double v11 = src_ch[x1 + y1 * in_w];
          
          double val = (1.0 - fx) * ((1.0 - fy) * v00 + fy * v01) +
                        fx * ((1.0 - fy) * v10 + fy * v11);
          
          if (val < 0.0) val = 0.0;
          else if (val > 1.0 && val < 255.0) val = std::min(1.0, val);
          
          dst_row[xo] = val;
        }
      }
    }
    
    if (nch == 1) {
      out_vec.attr("dim") = IntegerVector::create(out_w, out_h);
      out_vec.attr("colormode") = "Grayscale";
    } else {
      out_vec.attr("dim") = IntegerVector::create(out_w, out_h, nch);
      out_vec.attr("colormode") = "Color";
    }
    out_vec.attr("class") = CharacterVector::create("image", "Image", "array");
    out_vec.attr("colspace") = "RGB";
    out_vec.attr("gamma") = 1.0;
    return out_vec;
  }
}

//' Native C++ Sub-Pixel Image Stack Alignment (ECC Maximization)
//'
//' High-performance C++ implementation of Enhanced Correlation Coefficient (ECC)
//' maximization with multi-scale Gaussian pyramids for sub-pixel image alignment.
//'
//' @param images List of `Image` objects or arrays.
//' @param ref_idx 0-indexed reference image index.
//' @param method Character: `"ecc"` (affine) or `"translation"`.
//' @param crop Logical: automatically crop empty affine borders.
//' @param levels Integer: number of pyramid levels (default: 3).
//' @param max_iter Integer: max iterations per level (default: 50).
//' @param eps Numeric: convergence tolerance (default: 1e-4).
//' @param verbose Logical: print progress.
//'
//' @return A named list containing `aligned`, `matrices`, `statuses`, and `crop_box`.
//' @export
// [[Rcpp::export]]
Rcpp::List align_stack_cpp(Rcpp::List images,
                           int ref_idx = 0,
                           std::string method = "ecc",
                           bool crop = true,
                           int levels = 3,
                           int max_iter = 50,
                           double eps = 1e-4,
                           bool verbose = true) {
  int n = images.size();
  if (n == 0) {
    stop("align_stack_cpp: Empty image list.");
  }
  
  if (ref_idx < 0 || ref_idx >= n) {
    ref_idx = n / 2;
  }
  
  bool is_affine = (method == "ecc" || method == "affine");
  
  // Extract reference image grayscale
  int ref_w = 0, ref_h = 0, ref_c = 0;
  std::vector<double> ref_gray = extract_grayscale(images[ref_idx], ref_w, ref_h, ref_c);
  
  // Auto-tune levels if needed
  int min_dim = std::min(ref_w, ref_h);
  int effective_levels = levels;
  while (effective_levels > 1 && (min_dim >> (effective_levels - 1)) < 32) {
    effective_levels--;
  }
  if (effective_levels < 1) effective_levels = 1;
  
  std::vector<PyramidLevel> pyr_ref = build_pyramid(ref_gray.data(), ref_w, ref_h, effective_levels);
  
  std::vector<arma::mat> matrices(n);
  std::vector<std::string> statuses(n);
  
  for (int i = 0; i < n; ++i) {
    if (i == ref_idx) {
      matrices[i] = arma::eye(2, 3);
      statuses[i] = "reference";
      continue;
    }
    
    int src_w = 0, src_h = 0, src_c = 0;
    std::vector<double> src_gray = extract_grayscale(images[i], src_w, src_h, src_c);
    
    if (src_w != ref_w || src_h != ref_h) {
      matrices[i] = arma::eye(2, 3);
      statuses[i] = "dimension_mismatch";
      continue;
    }
    
    std::vector<PyramidLevel> pyr_src = build_pyramid(src_gray.data(), src_w, src_h, effective_levels);
    
    arma::mat W = arma::eye(2, 3);
    bool ok = true;
    
    for (int l = effective_levels - 1; l >= 0; --l) {
      if (l < effective_levels - 1) {
        W(0, 2) *= 2.0;
        W(1, 2) *= 2.0;
      }
      
      bool step_ok = ecc_single_scale(pyr_ref[l], pyr_src[l], W, is_affine, max_iter, eps);
      if (!step_ok && l == 0) {
        ok = false;
      }
    }
    
    if (ok) {
      matrices[i] = W;
      statuses[i] = "success";
    } else {
      matrices[i] = arma::eye(2, 3);
      statuses[i] = "fallback_identity";
    }
  }
  
  // Compute border-free crop coordinates across all aligned images
  int crop_x1 = 0;
  int crop_y1 = 0;
  int crop_x2 = ref_w - 1;
  int crop_y2 = ref_h - 1;
  
  if (crop) {
    double min_x = 0.0;
    double max_x = ref_w - 1.0;
    double min_y = 0.0;
    double max_y = ref_h - 1.0;
    
    for (int i = 0; i < n; ++i) {
      if (statuses[i] == "fallback_identity" || statuses[i] == "dimension_mismatch") continue;
      
      arma::mat W_inv;
      if (!invert_affine_2x3(matrices[i], W_inv)) continue;
      
      // Source corners: (0, 0), (w-1, 0), (w-1, h-1), (0, h-1)
      double c_src[4][2] = {
        {0.0, 0.0},
        {(double)(ref_w - 1), 0.0},
        {(double)(ref_w - 1), (double)(ref_h - 1)},
        {0.0, (double)(ref_h - 1)}
      };
      
      double c_ref[4][2];
      for (int k = 0; k < 4; ++k) {
        c_ref[k][0] = W_inv(0, 0) * c_src[k][0] + W_inv(0, 1) * c_src[k][1] + W_inv(0, 2);
        c_ref[k][1] = W_inv(1, 0) * c_src[k][0] + W_inv(1, 1) * c_src[k][1] + W_inv(1, 2);
      }
      
      double x_left = std::max(c_ref[0][0], c_ref[3][0]);
      double x_right = std::min(c_ref[1][0], c_ref[2][0]);
      double y_top = std::max(c_ref[0][1], c_ref[1][1]);
      double y_bottom = std::min(c_ref[2][1], c_ref[3][1]);
      
      min_x = std::max(min_x, x_left);
      max_x = std::min(max_x, x_right);
      min_y = std::max(min_y, y_top);
      max_y = std::min(max_y, y_bottom);
    }
    
    int cx1 = std::max(0, (int)std::ceil(min_x));
    int cy1 = std::max(0, (int)std::ceil(min_y));
    int cx2 = std::min(ref_w - 1, (int)std::floor(max_x));
    int cy2 = std::min(ref_h - 1, (int)std::floor(max_y));
    
    if (cx2 - cx1 >= 32 && cy2 - cy1 >= 32) {
      crop_x1 = cx1;
      crop_y1 = cy1;
      crop_x2 = cx2;
      crop_y2 = cy2;
    }
  }
  
  int out_w = crop_x2 - crop_x1 + 1;
  int out_h = crop_y2 - crop_y1 + 1;
  
  // Warp and crop all images
  Rcpp::List aligned_stack(n);
  for (int i = 0; i < n; ++i) {
    aligned_stack[i] = warp_and_crop_image(images[i], matrices[i], crop_x1, crop_y1, out_w, out_h);
  }
  
  // Package matrices into 3D array [2, 3, n]
  NumericVector mats_arr(Dimension(2, 3, n));
  for (int i = 0; i < n; ++i) {
    for (int r = 0; r < 2; ++r) {
      for (int c = 0; c < 3; ++c) {
        mats_arr[r + c * 2 + i * 6] = matrices[i](r, c);
      }
    }
  }
  
  // Return crop box as 1-indexed for R: c(xmin, ymin, xmax, ymax)
  IntegerVector crop_box = IntegerVector::create(crop_x1 + 1, crop_y1 + 1, crop_x2 + 1, crop_y2 + 1);
  
  return List::create(
    Named("aligned") = aligned_stack,
    Named("matrices") = mats_arr,
    Named("statuses") = wrap(statuses),
    Named("crop_box") = crop_box
  );
}

// Separable box blur for focus measure smoothing
static void box_blur_2d_slice(const double* src, double* dst, int w, int h, int r) {
  std::vector<double> tmp(w * h);
  double inv_len = 1.0 / (2 * r + 1);
  
  for (int y = 0; y < h; ++y) {
    const double* s_row = src + y * w;
    double* t_row = tmp.data() + y * w;
    double acc = 0.0;
    for (int x = -r; x <= r; ++x) {
      int cx = std::max(0, std::min(w - 1, x));
      acc += s_row[cx];
    }
    t_row[0] = acc * inv_len;
    for (int x = 1; x < w; ++x) {
      int prev_x = std::max(0, x - 1 - r);
      int next_x = std::min(w - 1, x + r);
      acc += s_row[next_x] - s_row[prev_x];
      t_row[x] = acc * inv_len;
    }
  }
  
  for (int x = 0; x < w; ++x) {
    double acc = 0.0;
    for (int y = -r; y <= r; ++y) {
      int cy = std::max(0, std::min(h - 1, y));
      acc += tmp[x + cy * w];
    }
    dst[x] = acc * inv_len;
    for (int y = 1; y < h; ++y) {
      int prev_y = std::max(0, y - 1 - r);
      int next_y = std::min(h - 1, y + r);
      acc += tmp[x + next_y * w] - tmp[x + prev_y * w];
      dst[x + y * w] = acc * inv_len;
    }
  }
}

//' Native C++ Extended Depth of Field (EDF) Focus Stacking via Modified Laplacian
//'
//' Ultra-fast analytical multi-focus image fusion using spatial Modified Laplacian (ML)
//' energy maps. Ideal for large image stacks (50-500 images) where deep learning
//' networks are computationally heavy.
//'
//' @param images List of aligned `Image` objects or arrays.
//' @param blend Character: `"soft"` (weighted softmax blending) or `"hard"` (maximum sharpness selection).
//' @param radius Integer: smoothing window radius for local sharpness aggregation (default: 3).
//' @param power Numeric: exponent weight for soft blending (default: 6.0).
//'
//' @return A fused `Image` object.
//' @export
// [[Rcpp::export]]
SEXP fuse_laplacian_cpp(Rcpp::List images,
                        std::string blend = "soft",
                        int radius = 3,
                        double power = 6.0) {
  int n = images.size();
  if (n == 0) stop("fuse_laplacian_cpp: Empty images list.");
  if (n == 1) return images[0];
  
  IntegerVector d0 = Rf_getAttrib(images[0], R_DimSymbol);
  int w = d0[0];
  int h = d0[1];
  int nch = (d0.size() >= 3) ? d0[2] : 1;
  int npix = w * h;
  int in_type = TYPEOF(images[0]);
  
  std::vector<std::vector<double>> focus_maps(n, std::vector<double>(npix));
  
  for (int k = 0; k < n; ++k) {
    std::vector<double> lum(npix);
    SEXP im_k = images[k];
    
    if (TYPEOF(im_k) == RAWSXP) {
      Rcpp::RawVector rv(im_k);
      const Rbyte* p = rv.begin();
      if (nch >= 3) {
        const Rbyte* pr = p;
        const Rbyte* pg = p + npix;
        const Rbyte* pb = p + 2 * npix;
        for (int i = 0; i < npix; ++i) {
          lum[i] = (0.299 * pr[i] + 0.587 * pg[i] + 0.114 * pb[i]) / 255.0;
        }
      } else {
        for (int i = 0; i < npix; ++i) lum[i] = p[i] / 255.0;
      }
    } else {
      Rcpp::NumericVector nv(im_k);
      const double* p = nv.begin();
      if (nch >= 3) {
        const double* pr = p;
        const double* pg = p + npix;
        const double* pb = p + 2 * npix;
        for (int i = 0; i < npix; ++i) {
          lum[i] = 0.299 * pr[i] + 0.587 * pg[i] + 0.114 * pb[i];
        }
      } else {
        for (int i = 0; i < npix; ++i) lum[i] = p[i];
      }
    }
    
    std::vector<double> ml(npix, 0.0);
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        int xm1 = std::max(0, x - 1);
        int xp1 = std::min(w - 1, x + 1);
        int ym1 = std::max(0, y - 1);
        int yp1 = std::min(h - 1, y + 1);
        
        double v = lum[x + y * w];
        double dxx = std::abs(2.0 * v - lum[xm1 + y * w] - lum[xp1 + y * w]);
        double dyy = std::abs(2.0 * v - lum[x + ym1 * w] - lum[x + yp1 * w]);
        ml[x + y * w] = (dxx + dyy) * (dxx + dyy);
      }
    }
    
    box_blur_2d_slice(ml.data(), focus_maps[k].data(), w, h, std::max(1, radius));
  }
  
  bool is_hard = (blend == "hard");
  
  if (in_type == RAWSXP) {
    std::vector<const Rbyte*> raw_ptrs(n);
    for (int k = 0; k < n; ++k) {
      Rcpp::RawVector rv(images[k]);
      raw_ptrs[k] = rv.begin();
    }
    
    Rcpp::RawVector out = Rcpp::no_init((std::size_t)w * h * nch);
    Rbyte* dst = out.begin();
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < npix; ++i) {
      if (is_hard) {
        int best_k = 0;
        double max_s = focus_maps[0][i];
        for (int k = 1; k < n; ++k) {
          if (focus_maps[k][i] > max_s) {
            max_s = focus_maps[k][i];
            best_k = k;
          }
        }
        for (int c = 0; c < nch; ++c) {
          dst[i + c * npix] = raw_ptrs[best_k][i + c * npix];
        }
      } else {
        double max_s = 1e-12;
        for (int k = 0; k < n; ++k) {
          if (focus_maps[k][i] > max_s) max_s = focus_maps[k][i];
        }
        
        double sum_w = 0.0;
        std::vector<double> weights(n);
        for (int k = 0; k < n; ++k) {
          double ratio = focus_maps[k][i] / max_s;
          double w_val = std::pow(ratio, power);
          weights[k] = w_val;
          sum_w += w_val;
        }
        
        double inv_sum = 1.0 / sum_w;
        for (int c = 0; c < nch; ++c) {
          double blended = 0.0;
          for (int k = 0; k < n; ++k) {
            blended += weights[k] * (double)raw_ptrs[k][i + c * npix];
          }
          int ival = (int)std::round(blended * inv_sum);
          if (ival < 0) ival = 0;
          else if (ival > 255) ival = 255;
          dst[i + c * npix] = (Rbyte)ival;
        }
      }
    }
    
    if (nch == 1) {
      out.attr("dim") = IntegerVector::create(w, h);
      out.attr("colormode") = "Grayscale";
    } else {
      out.attr("dim") = IntegerVector::create(w, h, nch);
      out.attr("colormode") = "Color";
    }
    out.attr("class") = CharacterVector::create("image", "Image", "array");
    out.attr("colspace") = "RGB";
    out.attr("gamma") = 1.0;
    return out;
  } else {
    std::vector<const double*> dbl_ptrs(n);
    for (int k = 0; k < n; ++k) {
      Rcpp::NumericVector nv(images[k]);
      dbl_ptrs[k] = nv.begin();
    }
    
    Rcpp::NumericVector out = Rcpp::no_init((std::size_t)w * h * nch);
    double* dst = out.begin();
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < npix; ++i) {
      if (is_hard) {
        int best_k = 0;
        double max_s = focus_maps[0][i];
        for (int k = 1; k < n; ++k) {
          if (focus_maps[k][i] > max_s) {
            max_s = focus_maps[k][i];
            best_k = k;
          }
        }
        for (int c = 0; c < nch; ++c) {
          dst[i + c * npix] = dbl_ptrs[best_k][i + c * npix];
        }
      } else {
        double max_s = 1e-12;
        for (int k = 0; k < n; ++k) {
          if (focus_maps[k][i] > max_s) max_s = focus_maps[k][i];
        }
        
        double sum_w = 0.0;
        std::vector<double> weights(n);
        for (int k = 0; k < n; ++k) {
          double ratio = focus_maps[k][i] / max_s;
          double w_val = std::pow(ratio, power);
          weights[k] = w_val;
          sum_w += w_val;
        }
        
        double inv_sum = 1.0 / sum_w;
        for (int c = 0; c < nch; ++c) {
          double blended = 0.0;
          for (int k = 0; k < n; ++k) {
            blended += weights[k] * dbl_ptrs[k][i + c * npix];
          }
          dst[i + c * npix] = blended * inv_sum;
        }
      }
    }
    
    if (nch == 1) {
      out.attr("dim") = IntegerVector::create(w, h);
      out.attr("colormode") = "Grayscale";
    } else {
      out.attr("dim") = IntegerVector::create(w, h, nch);
      out.attr("colormode") = "Color";
    }
    out.attr("class") = CharacterVector::create("image", "Image", "array");
    out.attr("colspace") = "RGB";
    out.attr("gamma") = 1.0;
    return out;
  }
}
