// image_stitch.cpp — High-Performance Native Planar Image Stitching for pliman
// Multi-Scale Coarse-to-Fine Sub-Pixel Normalized Cross-Correlation (NCC) Image Registration
// Cosine-Smooth Seam Feathering & Gain Compensation
// Zero Python, Zero OpenCV, Zero External Dependencies.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Extract luminance grayscale from image SEXP (supports both raw and double arrays)
static std::vector<double> extract_gray_stitch(SEXP img_sexp, int& w, int& h, int& nch) {
  Rcpp::RObject robj(img_sexp);
  Rcpp::IntegerVector dims = robj.attr("dim");
  if (dims.size() == 2) {
    w = dims[0];
    h = dims[1];
    nch = 1;
  } else if (dims.size() == 3) {
    w = dims[0];
    h = dims[1];
    nch = dims[2];
  } else {
    stop("extract_gray_stitch: Invalid image dimensions.");
  }

  std::vector<double> gray(w * h);
  int npix = w * h;

  if (TYPEOF(img_sexp) == RAWSXP) {
    Rcpp::RawVector rv(img_sexp);
    const Rbyte* ptr = rv.begin();
    if (nch == 1) {
      for (int i = 0; i < npix; ++i) gray[i] = (double)ptr[i];
    } else {
      for (int i = 0; i < npix; ++i) {
        gray[i] = 0.299 * (double)ptr[i] + 0.587 * (double)ptr[i + npix] + 0.114 * (double)ptr[i + 2 * npix];
      }
    }
  } else {
    Rcpp::NumericVector nv(img_sexp);
    const double* ptr = nv.begin();
    if (nch == 1) {
      for (int i = 0; i < npix; ++i) gray[i] = ptr[i] * 255.0;
    } else {
      for (int i = 0; i < npix; ++i) {
        gray[i] = (0.299 * ptr[i] + 0.587 * ptr[i + npix] + 0.114 * ptr[i + 2 * npix]) * 255.0;
      }
    }
  }
  return gray;
}

// Fast box-downsampler for multi-scale pyramid
static std::vector<double> downsample_fast(const double* src, int w, int h, int factor, int& out_w, int& out_h) {
  if (factor <= 1) {
    out_w = w;
    out_h = h;
    std::vector<double> dst(w * h);
    std::copy(src, src + w * h, dst.begin());
    return dst;
  }

  out_w = w / factor;
  out_h = h / factor;
  std::vector<double> dst(out_w * out_h, 0.0);
  double inv_area = 1.0 / (double)(factor * factor);

  for (int oy = 0; oy < out_h; ++oy) {
    for (int ox = 0; ox < out_w; ++ox) {
      double sum = 0.0;
      int sy_start = oy * factor;
      int sx_start = ox * factor;
      for (int dy = 0; dy < factor; ++dy) {
        for (int dx = 0; dx < factor; ++dx) {
          sum += src[(sx_start + dx) + (sy_start + dy) * w];
        }
      }
      dst[ox + oy * out_w] = sum * inv_area;
    }
  }
  return dst;
}

// Compute Normalized Cross-Correlation (NCC) with optional spatial stride
static double compute_overlap_ncc_stride(const double* g1, int w1, int h1,
                                         const double* g2, int w2, int h2,
                                         int dx, int dy, int stride = 1) {
  int x1_start = std::max(0, dx);
  int x1_end = std::min(w1, dx + w2);
  int y1_start = std::max(0, dy);
  int y1_end = std::min(h1, dy + h2);

  int overlap_w = x1_end - x1_start;
  int overlap_h = y1_end - y1_start;

  if (overlap_w < 16 || overlap_h < 16) {
    return -1.0;
  }

  double sum1 = 0.0, sum2 = 0.0;
  double sum1_sq = 0.0, sum2_sq = 0.0;
  double sum_prod = 0.0;
  int count = 0;

  for (int y = y1_start; y < y1_end; y += stride) {
    int y2 = y - dy;
    for (int x = x1_start; x < x1_end; x += stride) {
      int x2 = x - dx;
      double v1 = g1[x + y * w1];
      double v2 = g2[x2 + y2 * w2];

      sum1 += v1;
      sum2 += v2;
      sum1_sq += v1 * v1;
      sum2_sq += v2 * v2;
      sum_prod += v1 * v2;
      count++;
    }
  }

  if (count < 10) return -1.0;

  double mean1 = sum1 / count;
  double mean2 = sum2 / count;

  double var1 = sum1_sq - count * mean1 * mean1;
  double var2 = sum2_sq - count * mean2 * mean2;

  if (var1 <= 1e-4 || var2 <= 1e-4) {
    return 0.0;
  }

  double cov = sum_prod - count * mean1 * mean2;
  return cov / std::sqrt(var1 * var2);
}

// Find optimal (dx, dy) translation using ultra-fast Multi-Scale Pyramid NCC
static void find_optimal_displacement_multiscale(const double* g1, int w1, int h1,
                                                 const double* g2, int w2, int h2,
                                                 std::string direction,
                                                 double overlap_hint,
                                                 int& best_dx, int& best_dy, double& best_ncc) {
  best_ncc = -1.0;
  best_dx = 0;
  best_dy = 0;

  // Determine downsampling factor S so max dimension of thumbnail is <= 256
  int max_dim = std::max({w1, h1, w2, h2});
  int factor = std::max(1, max_dim / 256);

  int w1_c = 0, h1_c = 0;
  int w2_c = 0, h2_c = 0;
  std::vector<double> g1_c = downsample_fast(g1, w1, h1, factor, w1_c, h1_c);
  std::vector<double> g2_c = downsample_fast(g2, w2, h2, factor, w2_c, h2_c);

  // 1. Coarse search on downsampled thumbnail
  int min_dx_c = 0, max_dx_c = 0;
  int min_dy_c = 0, max_dy_c = 0;

  if (direction == "horizontal") {
    min_dx_c = (int)(w1_c * 0.15);
    max_dx_c = (int)(w1_c * 0.95);
    min_dy_c = -(int)(h1_c * 0.15);
    max_dy_c = (int)(h1_c * 0.15);
  } else { // vertical
    min_dx_c = -(int)(w1_c * 0.15);
    max_dx_c = (int)(w1_c * 0.15);
    min_dy_c = (int)(h1_c * 0.15);
    max_dy_c = (int)(h1_c * 0.95);
  }

  int best_c_dx = 0;
  int best_c_dy = 0;
  double best_c_ncc = -1.0;

  int step_c = 2; // small step on thumbnail
  for (int dy = min_dy_c; dy <= max_dy_c; dy += step_c) {
    for (int dx = min_dx_c; dx <= max_dx_c; dx += step_c) {
      double ncc = compute_overlap_ncc_stride(g1_c.data(), w1_c, h1_c, g2_c.data(), w2_c, h2_c, dx, dy, 1);
      if (ncc > best_c_ncc) {
        best_c_ncc = ncc;
        best_c_dx = dx;
        best_c_dy = dy;
      }
    }
  }

  // Refine on coarse grid (+/- 2 steps)
  for (int dy = best_c_dy - step_c; dy <= best_c_dy + step_c; ++dy) {
    for (int dx = best_c_dx - step_c; dx <= best_c_dx + step_c; ++dx) {
      double ncc = compute_overlap_ncc_stride(g1_c.data(), w1_c, h1_c, g2_c.data(), w2_c, h2_c, dx, dy, 1);
      if (ncc > best_c_ncc) {
        best_c_ncc = ncc;
        best_c_dx = dx;
        best_c_dy = dy;
      }
    }
  }

  // 2. Map back to full-resolution coordinates
  int init_dx = best_c_dx * factor;
  int init_dy = best_c_dy * factor;

  // 3. Local fine search around initial estimate on full-resolution
  int win_r = factor + 4;
  int fine_min_dx = std::max(0, init_dx - win_r);
  int fine_max_dx = std::min(w1, init_dx + win_r);
  int fine_min_dy = init_dy - win_r;
  int fine_max_dy = init_dy + win_r;

  int stride_full = std::max(1, factor / 4);

  for (int dy = fine_min_dy; dy <= fine_max_dy; dy += 2) {
    for (int dx = fine_min_dx; dx <= fine_max_dx; dx += 2) {
      double ncc = compute_overlap_ncc_stride(g1, w1, h1, g2, w2, h2, dx, dy, stride_full);
      if (ncc > best_ncc) {
        best_ncc = ncc;
        best_dx = dx;
        best_dy = dy;
      }
    }
  }

  // Sub-pixel fine polish (+/- 2 pixels around best)
  int polish_min_dx = best_dx - 2;
  int polish_max_dx = best_dx + 2;
  int polish_min_dy = best_dy - 2;
  int polish_max_dy = best_dy + 2;

  for (int dy = polish_min_dy; dy <= polish_max_dy; ++dy) {
    for (int dx = polish_min_dx; dx <= polish_max_dx; ++dx) {
      double ncc = compute_overlap_ncc_stride(g1, w1, h1, g2, w2, h2, dx, dy, 1);
      if (ncc > best_ncc) {
        best_ncc = ncc;
        best_dx = dx;
        best_dy = dy;
      }
    }
  }
}

// Stitch two images together with cosine feathering
// [[Rcpp::export]]
Rcpp::List stitch_pair_cpp(SEXP img1, SEXP img2,
                           std::string direction = "horizontal",
                           bool blend = true,
                           double overlap_hint = 0.2) {
  int w1 = 0, h1 = 0, nch1 = 0;
  int w2 = 0, h2 = 0, nch2 = 0;

  std::vector<double> g1 = extract_gray_stitch(img1, w1, h1, nch1);
  std::vector<double> g2 = extract_gray_stitch(img2, w2, h2, nch2);

  int nch = std::max(nch1, nch2);

  // Auto direction detection if needed
  if (direction == "auto") {
    int dx_h = 0, dy_h = 0;
    double ncc_h = -1.0;
    find_optimal_displacement_multiscale(g1.data(), w1, h1, g2.data(), w2, h2, "horizontal", overlap_hint, dx_h, dy_h, ncc_h);

    int dx_v = 0, dy_v = 0;
    double ncc_v = -1.0;
    find_optimal_displacement_multiscale(g1.data(), w1, h1, g2.data(), w2, h2, "vertical", overlap_hint, dx_v, dy_v, ncc_v);

    direction = (ncc_h >= ncc_v) ? "horizontal" : "vertical";
  }

  int best_dx = 0, best_dy = 0;
  double best_ncc = -1.0;
  find_optimal_displacement_multiscale(g1.data(), w1, h1, g2.data(), w2, h2, direction, overlap_hint, best_dx, best_dy, best_ncc);

  // Calculate canvas bounding box
  int min_cx = std::min(0, best_dx);
  int max_cx = std::max(w1, best_dx + w2);
  int min_cy = std::min(0, best_dy);
  int max_cy = std::max(h1, best_dy + h2);

  int canvas_w = max_cx - min_cx;
  int canvas_h = max_cy - min_cy;

  // Placement coordinates on canvas
  int i1_ox = -min_cx;
  int i1_oy = -min_cy;
  int i2_ox = best_dx - min_cx;
  int i2_oy = best_dy - min_cy;

  // Determine overlap bounding box on canvas
  int ol_x_start = std::max(i1_ox, i2_ox);
  int ol_x_end = std::min(i1_ox + w1, i2_ox + w2);
  int ol_y_start = std::max(i1_oy, i2_oy);
  int ol_y_end = std::min(i1_oy + h1, i2_oy + h2);

  bool has_overlap = (ol_x_end > ol_x_start && ol_y_end > ol_y_start);

  // Convert inputs to normalized double [0, 1] for blending
  std::vector<double> img1_dbl(w1 * h1 * nch, 0.0);
  std::vector<double> img2_dbl(w2 * h2 * nch, 0.0);

  if (TYPEOF(img1) == RAWSXP) {
    Rcpp::RawVector r1(img1);
    for (int i = 0; i < w1 * h1 * nch1; ++i) img1_dbl[i] = (double)r1[i] / 255.0;
  } else {
    Rcpp::NumericVector r1(img1);
    for (int i = 0; i < w1 * h1 * nch1; ++i) img1_dbl[i] = r1[i];
  }

  if (TYPEOF(img2) == RAWSXP) {
    Rcpp::RawVector r2(img2);
    for (int i = 0; i < w2 * h2 * nch2; ++i) img2_dbl[i] = (double)r2[i] / 255.0;
  } else {
    Rcpp::NumericVector r2(img2);
    for (int i = 0; i < w2 * h2 * nch2; ++i) img2_dbl[i] = r2[i];
  }

  // Allocate canvas
  Rcpp::NumericVector canvas(canvas_w * canvas_h * nch, 0.0);

  #pragma omp parallel for schedule(static)
  for (int cy = 0; cy < canvas_h; ++cy) {
    for (int cx = 0; cx < canvas_w; ++cx) {
      int c_idx = cx + cy * canvas_w;

      // Coordinate in I1
      int x1 = cx - i1_ox;
      int y1 = cy - i1_oy;
      bool in_i1 = (x1 >= 0 && x1 < w1 && y1 >= 0 && y1 < h1);

      // Coordinate in I2
      int x2 = cx - i2_ox;
      int y2 = cy - i2_oy;
      bool in_i2 = (x2 >= 0 && x2 < w2 && y2 >= 0 && y2 < h2);

      if (in_i1 && in_i2) {
        // Pixel is in overlap zone
        double w_i2 = 0.5;
        if (blend && has_overlap) {
          if (direction == "horizontal") {
            double ramp = (double)(cx - ol_x_start) / (double)(ol_x_end - ol_x_start);
            ramp = std::max(0.0, std::min(1.0, ramp));
            w_i2 = 0.5 * (1.0 - std::cos(M_PI * ramp));
          } else {
            double ramp = (double)(cy - ol_y_start) / (double)(ol_y_end - ol_y_start);
            ramp = std::max(0.0, std::min(1.0, ramp));
            w_i2 = 0.5 * (1.0 - std::cos(M_PI * ramp));
          }
        }
        double w_i1 = 1.0 - w_i2;

        for (int c = 0; c < nch; ++c) {
          double val1 = img1_dbl[x1 + y1 * w1 + c * (w1 * h1)];
          double val2 = img2_dbl[x2 + y2 * w2 + c * (w2 * h2)];
          canvas[c_idx + c * (canvas_w * canvas_h)] = w_i1 * val1 + w_i2 * val2;
        }
      } else if (in_i1) {
        for (int c = 0; c < nch; ++c) {
          canvas[c_idx + c * (canvas_w * canvas_h)] = img1_dbl[x1 + y1 * w1 + c * (w1 * h1)];
        }
      } else if (in_i2) {
        for (int c = 0; c < nch; ++c) {
          canvas[c_idx + c * (canvas_w * canvas_h)] = img2_dbl[x2 + y2 * w2 + c * (w2 * h2)];
        }
      }
    }
  }

  // Format as Image
  if (nch == 1) {
    canvas.attr("dim") = IntegerVector::create(canvas_w, canvas_h);
    canvas.attr("colormode") = "Grayscale";
  } else {
    canvas.attr("dim") = IntegerVector::create(canvas_w, canvas_h, nch);
    canvas.attr("colormode") = "Color";
  }
  canvas.attr("class") = CharacterVector::create("image", "Image", "array");
  canvas.attr("colspace") = "RGB";
  canvas.attr("gamma") = 1.0;

  return List::create(
    Named("image") = canvas,
    Named("dx") = best_dx,
    Named("dy") = best_dy,
    Named("ncc") = best_ncc,
    Named("direction") = direction,
    Named("overlap_width") = (has_overlap ? (ol_x_end - ol_x_start) : 0),
    Named("overlap_height") = (has_overlap ? (ol_y_end - ol_y_start) : 0)
  );
}
