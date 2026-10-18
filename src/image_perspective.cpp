// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#include <type_traits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// Invert a 3x3 matrix. Returns false if singular.
inline bool invert3x3(const double m[3][3], double inv[3][3]) {
  double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  if (std::fabs(det) < 1e-12) return false;

  double invDet = 1.0 / det;

  inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet;
  inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;

  inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
  inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * invDet;

  inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet;
  inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * invDet;
  inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * invDet;

  return true;
}

// Multiply two 3x3 matrices: C = A * B
inline void matMul3x3(const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      C[r][c] = A[r][0] * B[0][c] + A[r][1] * B[1][c] + A[r][2] * B[2][c];
    }
  }
}

template <int RTYPE>
SEXP do_warp_perspective(SEXP img_sexp, int W, int H, int nch,
                         int oW, int oH, const double Minv[3][3],
                         NumericVector bg_color) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  VecType img(img_sexp);
  const T* src = img.begin();

  std::vector<double> bg(nch, 1.0);
  for (int c = 0; c < std::min((int)bg_color.size(), nch); ++c)
    bg[c] = bg_color[c];

  const std::ptrdiff_t src_plane = (std::ptrdiff_t)W * H;
  const std::ptrdiff_t dst_plane = (std::ptrdiff_t)oW * oH;

  VecType out = Rcpp::no_init((std::size_t)oW * oH * nch);

  // Pre-fill with background color
  for (int ch = 0; ch < nch; ++ch) {
    T* dp = out.begin() + (std::ptrdiff_t)ch * dst_plane;
    T bg_val;
    if constexpr (std::is_same_v<T, Rbyte>) {
      bg_val = (T)std::max(0.0, std::min(255.0, std::round(bg[ch] * 255.0)));
    } else {
      bg_val = (T)bg[ch];
    }
    std::fill(dp, dp + dst_plane, bg_val);
  }

  const double m00 = Minv[0][0], m01 = Minv[0][1], m02 = Minv[0][2];
  const double m10 = Minv[1][0], m11 = Minv[1][1], m12 = Minv[1][2];
  const double m20 = Minv[2][0], m21 = Minv[2][1], m22 = Minv[2][2];

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const T* sp_ch = src + (std::ptrdiff_t)ch * src_plane;
    T*       dp_ch = out.begin() + (std::ptrdiff_t)ch * dst_plane;

    for (int oy = 0; oy < oH; ++oy) {
      double sx_h = m01 * oy + m02;
      double sy_h = m11 * oy + m12;
      double sw_h = m21 * oy + m22;

      T* row_dst = dp_ch + (std::ptrdiff_t)oy * oW;

      for (int ox = 0; ox < oW; ++ox, sx_h += m00, sy_h += m10, sw_h += m20) {
        if (sw_h <= 1e-9) continue; // Behind camera / vanishing horizon

        const double inv_w = 1.0 / sw_h;
        const double sx = sx_h * inv_w;
        const double sy = sy_h * inv_w;

        if (sx < 0.0 || sx > W - 1.0 || sy < 0.0 || sy > H - 1.0) {
          continue;
        }

        const int x0 = (int)sx;
        const int y0 = (int)sy;
        const int x1 = (x0 + 1 < W) ? (x0 + 1) : x0;
        const int y1 = (y0 + 1 < H) ? (y0 + 1) : y0;
        const double tx = sx - x0;
        const double ty = sy - y0;

        const double v00 = sp_ch[x0 + (std::ptrdiff_t)y0 * W];
        const double v10 = sp_ch[x1 + (std::ptrdiff_t)y0 * W];
        const double v01 = sp_ch[x0 + (std::ptrdiff_t)y1 * W];
        const double v11 = sp_ch[x1 + (std::ptrdiff_t)y1 * W];

        const double acc = (1.0 - tx) * (1.0 - ty) * v00 +
                           tx * (1.0 - ty) * v10 +
                           (1.0 - tx) * ty * v01 +
                           tx * ty * v11;

        if constexpr (std::is_same_v<T, Rbyte>) {
          row_dst[ox] = (T)std::max(0.0, std::min(255.0, std::round(acc)));
        } else {
          row_dst[ox] = (T)acc;
        }
      }
    }
  }

  out.attr("dim") = IntegerVector::create(oW, oH, nch);
  return out;
}

// [[Rcpp::export]]
List perspective_matrix_cpp(int W, int H,
                            double theta_x, double theta_y, double theta_z,
                            double scale = 1.0, double focal_fac = 1.0,
                            bool keep_size = true) {
  const double cx_src = (W - 1) * 0.5;
  const double cy_src = (H - 1) * 0.5;

  const double ax = theta_x * M_PI / 180.0;
  const double ay = theta_y * M_PI / 180.0;
  const double az = theta_z * M_PI / 180.0;

  // Rx: pitch (around X axis)
  const double Rx[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, std::cos(ax), -std::sin(ax)},
    {0.0, std::sin(ax),  std::cos(ax)}
  };

  // Ry: yaw (around Y axis)
  const double Ry[3][3] = {
    { std::cos(ay), 0.0, std::sin(ay)},
    { 0.0,         1.0, 0.0},
    {-std::sin(ay), 0.0, std::cos(ay)}
  };

  // Rz: roll (around Z axis)
  const double Rz[3][3] = {
    {std::cos(az), -std::sin(az), 0.0},
    {std::sin(az),  std::cos(az), 0.0},
    {0.0,          0.0,          1.0}
  };

  double Ryx[3][3];
  matMul3x3(Ry, Rx, Ryx);
  double R[3][3];
  matMul3x3(Rz, Ryx, R);

  // Focal distance d: default is max(W, H) * focal_fac
  double max_dim = (double)std::max(W, H);
  double d = (focal_fac <= 0.0 ? max_dim : max_dim * focal_fac);

  // Projection matrix A centered at origin
  double A[3][3] = {
    {scale * d * R[0][0], scale * d * R[0][1], 0.0},
    {scale * d * R[1][0], scale * d * R[1][1], 0.0},
    {R[2][0],             R[2][1],             d}
  };

  double Tin[3][3] = {
    {1.0, 0.0, -cx_src},
    {0.0, 1.0, -cy_src},
    {0.0, 0.0, 1.0}
  };

  double M0[3][3];
  matMul3x3(A, Tin, M0);

  // Corner projection to find bounding box if keep_size == false
  const double corners[4][2] = {
    {0.0, 0.0},
    {(double)(W - 1), 0.0},
    {(double)(W - 1), (double)(H - 1)},
    {0.0, (double)(H - 1)}
  };

  double proj_x[4], proj_y[4];
  for (int k = 0; k < 4; ++k) {
    double x = corners[k][0], y = corners[k][1];
    double u1 = M0[0][0] * x + M0[0][1] * y + M0[0][2];
    double u2 = M0[1][0] * x + M0[1][1] * y + M0[1][2];
    double u3 = M0[2][0] * x + M0[2][1] * y + M0[2][2];
    if (std::fabs(u3) < 1e-9) u3 = 1e-9;
    proj_x[k] = u1 / u3;
    proj_y[k] = u2 / u3;
  }

  int oW = W, oH = H;
  double Tout[3][3];

  if (keep_size) {
    oW = W;
    oH = H;
    const double cx_dst = (oW - 1) * 0.5;
    const double cy_dst = (oH - 1) * 0.5;
    Tout[0][0] = 1.0; Tout[0][1] = 0.0; Tout[0][2] = cx_dst;
    Tout[1][0] = 0.0; Tout[1][1] = 1.0; Tout[1][2] = cy_dst;
    Tout[2][0] = 0.0; Tout[2][1] = 0.0; Tout[2][2] = 1.0;
  } else {
    double min_x = proj_x[0], max_x = proj_x[0];
    double min_y = proj_y[0], max_y = proj_y[0];
    for (int k = 1; k < 4; ++k) {
      if (proj_x[k] < min_x) min_x = proj_x[k];
      if (proj_x[k] > max_x) max_x = proj_x[k];
      if (proj_y[k] < min_y) min_y = proj_y[k];
      if (proj_y[k] > max_y) max_y = proj_y[k];
    }
    oW = (int)std::ceil(max_x - min_x + 1.0);
    oH = (int)std::ceil(max_y - min_y + 1.0);
    Tout[0][0] = 1.0; Tout[0][1] = 0.0; Tout[0][2] = -min_x;
    Tout[1][0] = 0.0; Tout[1][1] = 1.0; Tout[1][2] = -min_y;
    Tout[2][0] = 0.0; Tout[2][1] = 0.0; Tout[2][2] = 1.0;
  }

  double M[3][3];
  matMul3x3(Tout, M0, M);

  // Normalize so M[2][2] == 1.0 if possible
  if (std::fabs(M[2][2]) > 1e-12) {
    double s = 1.0 / M[2][2];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        M[r][c] *= s;
      }
    }
  }

  double Minv[3][3];
  bool ok = invert3x3(M, Minv);
  if (!ok) {
    Rcpp::stop("Singular perspective matrix (cannot invert).");
  }

  NumericMatrix M_rcpp(3, 3);
  NumericMatrix Minv_rcpp(3, 3);
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      M_rcpp(r, c) = M[r][c];
      Minv_rcpp(r, c) = Minv[r][c];
    }
  }

  return List::create(
    Named("M") = M_rcpp,
    Named("Minv") = Minv_rcpp,
    Named("oW") = oW,
    Named("oH") = oH
  );
}

// [[Rcpp::export]]
SEXP image_warp_perspective_cpp(SEXP img_sexp,
                                NumericMatrix Minv_mat,
                                int oW, int oH,
                                NumericVector bg_color) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);
  int W, H, nch;
  if (dims.size() == 2) {
    W = dims[0]; H = dims[1]; nch = 1;
  } else if (dims.size() == 3) {
    W = dims[0]; H = dims[1]; nch = dims[2];
  } else {
    Rcpp::stop("image_warp_perspective_cpp: input must be 2D matrix or 3D array.");
    return R_NilValue;
  }

  double Minv[3][3];
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      Minv[r][c] = Minv_mat(r, c);
    }
  }

  if (TYPEOF(img_sexp) == RAWSXP) {
    return do_warp_perspective<RAWSXP>(img_sexp, W, H, nch, oW, oH, Minv, bg_color);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    return do_warp_perspective<REALSXP>(img_sexp, W, H, nch, oW, oH, Minv, bg_color);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    return do_warp_perspective<INTSXP>(img_sexp, W, H, nch, oW, oH, Minv, bg_color);
  } else {
    Rcpp::stop("image_warp_perspective_cpp: unsupported input type.");
    return R_NilValue;
  }
}

// [[Rcpp::export]]
SEXP image_perspective_cpp(SEXP img_sexp,
                           double theta_x, double theta_y, double theta_z,
                           double scale, double focal_fac,
                           NumericVector bg_color, bool keep_size) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);
  int W, H;
  if (dims.size() == 2) {
    W = dims[0]; H = dims[1];
  } else if (dims.size() == 3) {
    W = dims[0]; H = dims[1];
  } else {
    Rcpp::stop("image_perspective_cpp: input must be 2D matrix or 3D array.");
    return R_NilValue;
  }

  List pmat = perspective_matrix_cpp(W, H, theta_x, theta_y, theta_z, scale, focal_fac, keep_size);
  NumericMatrix Minv = pmat["Minv"];
  NumericMatrix M = pmat["M"];
  int oW = pmat["oW"];
  int oH = pmat["oH"];

  SEXP out = image_warp_perspective_cpp(img_sexp, Minv, oW, oH, bg_color);
  Rf_setAttrib(out, Rf_install("homography"), M);
  return out;
}
