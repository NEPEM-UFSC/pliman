// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <type_traits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

template <int RTYPE>
SEXP rotate90(SEXP img_sexp, int W, int H, int nch) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  VecType img(img_sexp);
  const T* src = img.begin();
  
  VecType out = Rcpp::no_init((std::size_t)W * H * nch);
  const int oW = H, oH = W;
  const std::ptrdiff_t sp = (std::ptrdiff_t)W * H;
  const std::ptrdiff_t dp = (std::ptrdiff_t)oW * oH;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const T* sc = src + (std::ptrdiff_t)ch * sp;
    T*       dc = out.begin() + (std::ptrdiff_t)ch * dp;
    for (int oy = 0; oy < oH; ++oy) {
      for (int ox = 0; ox < oW; ++ox) {
        dc[ox + (std::ptrdiff_t)oy * oW] = sc[oy + (std::ptrdiff_t)(H - 1 - ox) * W];
      }
    }
  }
  out.attr("dim") = IntegerVector::create(oW, oH, nch);
  return out;
}

template <int RTYPE>
SEXP rotate180(SEXP img_sexp, int W, int H, int nch) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  VecType img(img_sexp);
  const T* src = img.begin();
  
  VecType out = Rcpp::no_init((std::size_t)W * H * nch);
  const std::ptrdiff_t sp = (std::ptrdiff_t)W * H;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const T* sc = src  + (std::ptrdiff_t)ch * sp;
    T*       dc = out.begin() + (std::ptrdiff_t)ch * sp;
    for (std::ptrdiff_t i = 0; i < sp; ++i)
      dc[i] = sc[sp - 1 - i];
  }
  out.attr("dim") = IntegerVector::create(W, H, nch);
  return out;
}

template <int RTYPE>
SEXP rotate270(SEXP img_sexp, int W, int H, int nch) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  VecType img(img_sexp);
  const T* src = img.begin();
  
  VecType out = Rcpp::no_init((std::size_t)W * H * nch);
  const int oW = H, oH = W;
  const std::ptrdiff_t sp = (std::ptrdiff_t)W * H;
  const std::ptrdiff_t dp = (std::ptrdiff_t)oW * oH;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const T* sc = src + (std::ptrdiff_t)ch * sp;
    T*       dc = out.begin() + (std::ptrdiff_t)ch * dp;
    for (int oy = 0; oy < oH; ++oy) {
      for (int ox = 0; ox < oW; ++ox) {
        dc[ox + (std::ptrdiff_t)oy * oW] = sc[(W - 1 - oy) + (std::ptrdiff_t)ox * W];
      }
    }
  }
  out.attr("dim") = IntegerVector::create(oW, oH, nch);
  return out;
}

template <int RTYPE>
SEXP do_rotate(SEXP img_sexp, int W, int H, int nch, double angle_deg, NumericVector bg_color) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  
  VecType img(img_sexp);
  const T* src = img.begin();

  std::vector<double> bg(nch, 1.0);
  for (int c = 0; c < std::min((int)bg_color.size(), nch); ++c)
    bg[c] = bg_color[c];

  const double rad = -angle_deg * M_PI / 180.0;
  const double cosA = std::cos(rad);
  const double sinA = std::sin(rad);

  const double cx_src = (W - 1) * 0.5;
  const double cy_src = (H - 1) * 0.5;

  const double hw = (W - 1) * 0.5, hh = (H - 1) * 0.5;
  const double corners_x[4] = { hw,  hw, -hw, -hw};
  const double corners_y[4] = { hh, -hh,  hh, -hh};
  double max_ox = 0, max_oy = 0;
  for (int k = 0; k < 4; ++k) {
    double rx = std::fabs(cosA * corners_x[k] - sinA * corners_y[k]);
    double ry = std::fabs(sinA * corners_x[k] + cosA * corners_y[k]);
    if (rx > max_ox) max_ox = rx;
    if (ry > max_oy) max_oy = ry;
  }
  const int oW = (int)std::ceil(max_ox * 2.0 + 1.0);
  const int oH = (int)std::ceil(max_oy * 2.0 + 1.0);

  const double cx_dst = (oW - 1) * 0.5;
  const double cy_dst = (oH - 1) * 0.5;

  const std::ptrdiff_t src_plane = (std::ptrdiff_t)W * H;
  const std::ptrdiff_t dst_plane = (std::ptrdiff_t)oW * oH;

  VecType out = Rcpp::no_init((std::size_t)oW * oH * nch);

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

  const double iCos =  cosA;
  const double iSin = -sinA;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const T* sp_ch = src + (std::ptrdiff_t)ch * src_plane;
    T*       dp_ch = out.begin() + (std::ptrdiff_t)ch * dst_plane;

    for (int oy = 0; oy < oH; ++oy) {
      const double dy = oy - cy_dst;
      double sx0 = iCos * (-cx_dst) + iSin * dy + cx_src;
      double sy0 = -iSin * (-cx_dst) + iCos * dy + cy_src;

      T* row_dst = dp_ch + (std::ptrdiff_t)oy * oW;

      for (int ox = 0; ox < oW; ++ox, sx0 += iCos, sy0 += (-iSin)) {
        const double sx = sx0;
        const double sy = sy0;

        if (sx < 0.0 || sx > W - 1.0 || sy < 0.0 || sy > H - 1.0) {
          continue;
        }

        const int x0 = (int)sx;
        const int y0 = (int)sy;
        const int x1 = x0 + 1 < W ? x0 + 1 : x0;
        const int y1 = y0 + 1 < H ? y0 + 1 : y0;
        const double tx = sx - x0;
        const double ty = sy - y0;

        const double v00 = sp_ch[x0 + (std::ptrdiff_t)y0 * W];
        const double v10 = sp_ch[x1 + (std::ptrdiff_t)y0 * W];
        const double v01 = sp_ch[x0 + (std::ptrdiff_t)y1 * W];
        const double v11 = sp_ch[x1 + (std::ptrdiff_t)y1 * W];

        double acc = (1.0 - ty) * ((1.0 - tx) * v00 + tx * v10)
                   +        ty  * ((1.0 - tx) * v01 + tx * v11);
                   
        if constexpr (std::is_same_v<T, Rbyte>) {
           if (acc < 0.0) acc = 0.0;
           else if (acc > 255.0) acc = 255.0;
           row_dst[ox] = (T)std::round(acc);
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
SEXP image_rotate_cpp(SEXP img_sexp,
                      double angle_deg,
                      NumericVector bg_color) {

  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);

  int W, H, nch;
  if (dims.size() == 2) {
    W = dims[0]; H = dims[1]; nch = 1;
  } else if (dims.size() == 3) {
    W = dims[0]; H = dims[1]; nch = dims[2];
  } else {
    Rcpp::stop("image_rotate_cpp: input must be a 2-D matrix or 3-D array.");
    return R_NilValue;
  }

  double norm = std::fmod(angle_deg, 360.0);
  if (norm < 0) norm += 360.0;

  if (std::fabs(norm - 0.0)   < 1e-9 || std::fabs(norm - 360.0) < 1e-9) {
    return Rcpp::clone(img_sexp);
  }

  if (TYPEOF(img_sexp) == RAWSXP) {
    if (std::fabs(norm - 90.0)  < 1e-9) { return rotate90<RAWSXP>(img_sexp, W, H, nch);  }
    if (std::fabs(norm - 180.0) < 1e-9) { return rotate180<RAWSXP>(img_sexp, W, H, nch); }
    if (std::fabs(norm - 270.0) < 1e-9) { return rotate270<RAWSXP>(img_sexp, W, H, nch); }
    return do_rotate<RAWSXP>(img_sexp, W, H, nch, angle_deg, bg_color);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    if (std::fabs(norm - 90.0)  < 1e-9) { return rotate90<REALSXP>(img_sexp, W, H, nch);  }
    if (std::fabs(norm - 180.0) < 1e-9) { return rotate180<REALSXP>(img_sexp, W, H, nch); }
    if (std::fabs(norm - 270.0) < 1e-9) { return rotate270<REALSXP>(img_sexp, W, H, nch); }
    return do_rotate<REALSXP>(img_sexp, W, H, nch, angle_deg, bg_color);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    if (std::fabs(norm - 90.0)  < 1e-9) { return rotate90<INTSXP>(img_sexp, W, H, nch);  }
    if (std::fabs(norm - 180.0) < 1e-9) { return rotate180<INTSXP>(img_sexp, W, H, nch); }
    if (std::fabs(norm - 270.0) < 1e-9) { return rotate270<INTSXP>(img_sexp, W, H, nch); }
    return do_rotate<INTSXP>(img_sexp, W, H, nch, angle_deg, bg_color);
  } else {
    Rcpp::stop("image_rotate_cpp: unsupported input type.");
    return R_NilValue;
  }
}
