// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

// =============================================================================
// FIELD OCR HELPERS (image_ocr_dl)
//
// Image convention (pliman / EBImage): the input array is W x H x C stored
// column-major, so pixel (x, y, c) lives at x + y * W + c * W * H.
// Output tensors are NCHW (row-major per plane): c * (H * W) + y * W + x.
// =============================================================================

namespace {

struct ImgView {
  int w = 0, h = 0, nch = 1;
  const uint8_t* raw = nullptr;
  const double* dbl = nullptr;
  const int* itg = nullptr;
  double scale = 1.0;  // multiplier that maps stored values to [0, 255]

  inline double at(int x, int y, int c) const {
    size_t idx = (size_t)x + (size_t)y * w + (size_t)c * w * h;
    if (raw) return (double)raw[idx];
    if (dbl) return dbl[idx] * scale;
    return (double)itg[idx] * scale;
  }

  // Returns the RGB value (0-255) of pixel (x, y); alpha is blended over white,
  // grayscale is replicated.
  inline void rgb(int x, int y, double* out) const {
    if (nch >= 3) {
      out[0] = at(x, y, 0); out[1] = at(x, y, 1); out[2] = at(x, y, 2);
      if (nch >= 4) {
        double a = at(x, y, 3) / 255.0;
        for (int k = 0; k < 3; ++k) out[k] = out[k] * a + 255.0 * (1.0 - a);
      }
    } else {
      double g = at(x, y, 0);
      if (nch == 2) {
        double a = at(x, y, 1) / 255.0;
        g = g * a + 255.0 * (1.0 - a);
      }
      out[0] = out[1] = out[2] = g;
    }
  }

  // Bilinear sample at continuous coordinate (fx, fy); pixel centres at integers.
  inline void sample(double fx, double fy, double* out) const {
    if (fx < 0) fx = 0; if (fx > w - 1) fx = w - 1;
    if (fy < 0) fy = 0; if (fy > h - 1) fy = h - 1;
    int x0 = (int)std::floor(fx), y0 = (int)std::floor(fy);
    int x1 = std::min(x0 + 1, w - 1), y1 = std::min(y0 + 1, h - 1);
    double ax = fx - x0, ay = fy - y0;
    double p00[3], p10[3], p01[3], p11[3];
    rgb(x0, y0, p00); rgb(x1, y0, p10); rgb(x0, y1, p01); rgb(x1, y1, p11);
    for (int k = 0; k < 3; ++k) {
      out[k] = (1 - ax) * (1 - ay) * p00[k] + ax * (1 - ay) * p10[k] +
               (1 - ax) * ay * p01[k] + ax * ay * p11[k];
    }
  }
};

ImgView make_view(SEXP img, double scale) {
  ImgView v;
  SEXP dims = Rf_getAttrib(img, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_length(dims) < 2) Rcpp::stop("ocr: image must be a 2D or 3D array.");
  v.w = INTEGER(dims)[0];
  v.h = INTEGER(dims)[1];
  v.nch = Rf_length(dims) >= 3 ? INTEGER(dims)[2] : 1;
  v.scale = scale;
  switch (TYPEOF(img)) {
  case RAWSXP:  v.raw = RAW(img); break;
  case REALSXP: v.dbl = REAL(img); break;
  case INTSXP:
  case LGLSXP:  v.itg = INTEGER(img); break;
  default: Rcpp::stop("ocr: unsupported image storage mode.");
  }
  return v;
}

} // namespace

// -----------------------------------------------------------------------------
// ocr_warp_crops_cpp
//
// Samples N quadrilateral regions of an image into a single NCHW float tensor.
// `quads` is an N x 8 matrix with the corners TL, TR, BR, BL given as
// (x1, y1, x2, y2, x3, y3, x4, y4) in 0-based continuous pixel coordinates
// (pixel centres at integers; the outer edge of the image is -0.5 / W - 0.5).
// Each crop i is rendered at out_h x out_w[i] (bilinear quad mapping, which is
// exact for rectangles/parallelograms) and left-aligned inside a canvas of
// width pad_w; the remainder is left at `pad_value` (after normalisation).
// Downscaling uses k x k supersampling (area-like filtering) so thin pen
// strokes are not aliased away.
// Normalisation per output channel j: (v / 255 - mean[j]) / std[j]; when
// `bgr = TRUE` channel 0 is Blue (Paddle convention).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector ocr_warp_crops_cpp(SEXP img,
                                       Rcpp::NumericMatrix quads,
                                       int out_h,
                                       Rcpp::IntegerVector out_w,
                                       int pad_w,
                                       Rcpp::NumericVector mean,
                                       Rcpp::NumericVector std,
                                       bool bgr = true,
                                       double scale = 1.0,
                                       double pad_value = 0.0) {
  ImgView v = make_view(img, scale);
  const int n = quads.nrow();
  if (quads.ncol() != 8) Rcpp::stop("ocr_warp_crops_cpp: `quads` must have 8 columns.");
  if (out_w.size() != n) Rcpp::stop("ocr_warp_crops_cpp: `out_w` must have one entry per quad.");
  if (mean.size() < 3 || std.size() < 3) Rcpp::stop("ocr_warp_crops_cpp: `mean`/`std` need 3 values.");

  const size_t plane = (size_t)out_h * pad_w;
  const size_t per = 3 * plane;
  Rcpp::NumericVector out(per * (size_t)n, pad_value);
  double* po = out.begin();

  std::vector<double> q(quads.begin(), quads.end());  // column-major N x 8
  std::vector<int> ow(out_w.begin(), out_w.end());
  const double m[3] = {mean[0], mean[1], mean[2]};
  const double s[3] = {std[0], std[1], std[2]};

  for (int i = 0; i < n; ++i) {
    const double x1 = q[i + 0 * n], y1 = q[i + 1 * n];
    const double x2 = q[i + 2 * n], y2 = q[i + 3 * n];
    const double x3 = q[i + 4 * n], y3 = q[i + 5 * n];
    const double x4 = q[i + 6 * n], y4 = q[i + 7 * n];
    const int cw = std::max(1, std::min(ow[i], pad_w));

    const double len_w = std::max(std::hypot(x2 - x1, y2 - y1), std::hypot(x3 - x4, y3 - y4));
    const double len_h = std::max(std::hypot(x4 - x1, y4 - y1), std::hypot(x3 - x2, y3 - y2));
    const int ks = std::max(1, std::min(6, (int)std::ceil(len_w / cw - 0.25)));
    const int kt = std::max(1, std::min(6, (int)std::ceil(len_h / out_h - 0.25)));
    const double inv = 1.0 / (ks * kt);
    double* base = po + per * (size_t)i;

    #pragma omp parallel for schedule(static) if((size_t)out_h * cw > 20000)
    for (int yy = 0; yy < out_h; ++yy) {
      for (int xx = 0; xx < cw; ++xx) {
        double acc[3] = {0, 0, 0};
        for (int b = 0; b < kt; ++b) {
          const double t = (yy + (b + 0.5) / kt) / out_h;
          for (int a = 0; a < ks; ++a) {
            const double u = (xx + (a + 0.5) / ks) / cw;
            // Bilinear interpolation of the quad corners (edges at -0.5 offsets
            // are encoded by the caller).
            const double px = (1 - u) * (1 - t) * x1 + u * (1 - t) * x2 + u * t * x3 + (1 - u) * t * x4;
            const double py = (1 - u) * (1 - t) * y1 + u * (1 - t) * y2 + u * t * y3 + (1 - u) * t * y4;
            double rgb[3];
            v.sample(px, py, rgb);
            acc[0] += rgb[0]; acc[1] += rgb[1]; acc[2] += rgb[2];
          }
        }
        const size_t off = (size_t)yy * pad_w + xx;
        for (int j = 0; j < 3; ++j) {
          const int src = bgr ? (2 - j) : j;
          base[(size_t)j * plane + off] = (acc[src] * inv / 255.0 - m[j]) / s[j];
        }
      }
    }
  }
  return out;
}

// -----------------------------------------------------------------------------
// ocr_components_cpp
//
// Single-pass statistics of the labelled connected components of a text
// probability map (replaces an O(pixels x components) R loop). For each label
// it returns the pixel count, mean / max score, axis-aligned bounds and the
// oriented rectangle derived from second-order moments:
//   u = (ux, uy): principal direction (text direction), v = (-uy, ux)
//   [u_min, u_max] x [v_min, v_max]: extents of the pixels along u and v,
//   measured from the centroid (cx, cy).
// Components whose principal axis is closer to vertical are re-oriented
// horizontally unless they are strongly elongated (`vertical_ratio`), so
// isolated tall characters (e.g. "1", "7") are not read sideways.
// Coordinates are 0-based, in the map's pixel grid (x = row, y = column).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::DataFrame ocr_components_cpp(Rcpp::IntegerMatrix labels,
                                   Rcpp::NumericMatrix score,
                                   int n_labels,
                                   double vertical_ratio = 3.0,
                                   double snap_deg = 2.0) {
  const int W = labels.nrow(), H = labels.ncol();
  if (score.nrow() != W || score.ncol() != H) Rcpp::stop("ocr_components_cpp: dimension mismatch.");
  const int K = std::max(0, n_labels);
  std::vector<double> cnt(K + 1, 0), ssum(K + 1, 0), smax(K + 1, -1e300);
  std::vector<double> sx(K + 1, 0), sy(K + 1, 0), sxx(K + 1, 0), syy(K + 1, 0), sxy(K + 1, 0);
  std::vector<int> xmin(K + 1, W), xmax(K + 1, -1), ymin(K + 1, H), ymax(K + 1, -1);
  const int* L = labels.begin();
  const double* S = score.begin();

  for (int y = 0; y < H; ++y) {
    const size_t col = (size_t)y * W;
    for (int x = 0; x < W; ++x) {
      const int id = L[col + x];
      if (id <= 0 || id > K) continue;
      const double sc = S[col + x];
      cnt[id] += 1; ssum[id] += sc; if (sc > smax[id]) smax[id] = sc;
      sx[id] += x; sy[id] += y;
      sxx[id] += (double)x * x; syy[id] += (double)y * y; sxy[id] += (double)x * y;
      if (x < xmin[id]) xmin[id] = x; if (x > xmax[id]) xmax[id] = x;
      if (y < ymin[id]) ymin[id] = y; if (y > ymax[id]) ymax[id] = y;
    }
  }

  const double snap = snap_deg * M_PI / 180.0;
  std::vector<double> cx(K + 1, 0), cy(K + 1, 0), ux(K + 1, 1), uy(K + 1, 0);
  for (int k = 1; k <= K; ++k) {
    if (cnt[k] <= 0) continue;
    cx[k] = sx[k] / cnt[k]; cy[k] = sy[k] / cnt[k];
    const double vxx = sxx[k] / cnt[k] - cx[k] * cx[k];
    const double vyy = syy[k] / cnt[k] - cy[k] * cy[k];
    const double vxy = sxy[k] / cnt[k] - cx[k] * cy[k];
    double theta = 0.5 * std::atan2(2.0 * vxy, vxx - vyy);  // principal axis angle
    const double tr = vxx + vyy;
    const double dt = std::sqrt(std::max(0.0, (vxx - vyy) * (vxx - vyy) / 4.0 + vxy * vxy));
    const double l1 = tr / 2.0 + dt, l2 = std::max(1e-9, tr / 2.0 - dt);
    const double elong = std::sqrt(l1 / l2);
    if (elong < 1.15) theta = 0.0;  // near-isotropic blob: keep horizontal
    if (std::fabs(theta) > M_PI / 4.0 && elong < vertical_ratio) {
      theta += (theta > 0 ? -M_PI / 2.0 : M_PI / 2.0);  // read tall glyphs upright
    }
    if (std::fabs(theta) < snap) theta = 0.0;
    double c = std::cos(theta), s = std::sin(theta);
    if (c < 0 || (std::fabs(c) < 1e-9 && s < 0)) { c = -c; s = -s; }
    ux[k] = c; uy[k] = s;
  }

  // Second pass: extents along the oriented axes.
  std::vector<double> umin(K + 1, 1e300), umax(K + 1, -1e300), vmin(K + 1, 1e300), vmax(K + 1, -1e300);
  for (int y = 0; y < H; ++y) {
    const size_t col = (size_t)y * W;
    for (int x = 0; x < W; ++x) {
      const int id = L[col + x];
      if (id <= 0 || id > K) continue;
      const double dx = x - cx[id], dy = y - cy[id];
      const double pu = dx * ux[id] + dy * uy[id];
      const double pv = -dx * uy[id] + dy * ux[id];
      if (pu < umin[id]) umin[id] = pu; if (pu > umax[id]) umax[id] = pu;
      if (pv < vmin[id]) vmin[id] = pv; if (pv > vmax[id]) vmax[id] = pv;
    }
  }

  std::vector<int> ids;
  for (int k = 1; k <= K; ++k) if (cnt[k] > 0) ids.push_back(k);
  const int m = (int)ids.size();
  Rcpp::IntegerVector o_id(m), o_n(m), o_xmin(m), o_xmax(m), o_ymin(m), o_ymax(m);
  Rcpp::NumericVector o_mean(m), o_max(m), o_cx(m), o_cy(m), o_ux(m), o_uy(m),
                      o_umin(m), o_umax(m), o_vmin(m), o_vmax(m);
  for (int i = 0; i < m; ++i) {
    const int k = ids[i];
    o_id[i] = k; o_n[i] = (int)cnt[k];
    o_mean[i] = ssum[k] / cnt[k]; o_max[i] = smax[k];
    o_xmin[i] = xmin[k]; o_xmax[i] = xmax[k]; o_ymin[i] = ymin[k]; o_ymax[i] = ymax[k];
    o_cx[i] = cx[k]; o_cy[i] = cy[k]; o_ux[i] = ux[k]; o_uy[i] = uy[k];
    o_umin[i] = umin[k]; o_umax[i] = umax[k]; o_vmin[i] = vmin[k]; o_vmax[i] = vmax[k];
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("label") = o_id, Rcpp::Named("n") = o_n,
    Rcpp::Named("mean_score") = o_mean, Rcpp::Named("max_score") = o_max,
    Rcpp::Named("xmin") = o_xmin, Rcpp::Named("xmax") = o_xmax,
    Rcpp::Named("ymin") = o_ymin, Rcpp::Named("ymax") = o_ymax,
    Rcpp::Named("cx") = o_cx, Rcpp::Named("cy") = o_cy,
    Rcpp::Named("ux") = o_ux, Rcpp::Named("uy") = o_uy,
    Rcpp::Named("u_min") = o_umin, Rcpp::Named("u_max") = o_umax,
    Rcpp::Named("v_min") = o_vmin, Rcpp::Named("v_max") = o_vmax
  );
}
