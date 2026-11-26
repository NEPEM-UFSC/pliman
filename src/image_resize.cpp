// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <type_traits>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static const double LANCZOS_A = 3.0;

static inline double sinc(double x) {
  if (std::fabs(x) < 1e-12) return 1.0;
  const double px = M_PI * x;
  return std::sin(px) / px;
}

static inline double lanczos3(double x) {
  if (std::fabs(x) >= LANCZOS_A) return 0.0;
  return sinc(x) * sinc(x / LANCZOS_A);
}

struct FilterEntry {
  int    start;
  int    count;
};

struct FilterTable {
  std::vector<FilterEntry> entries;
  std::vector<double>      weights;
  std::vector<int>         offsets;
};

static FilterTable build_filter_table(int src_len, int dst_len, bool use_lanczos) {
  FilterTable ft;
  ft.entries.resize(dst_len);
  ft.offsets.resize(dst_len);

  const double ratio  = (double)src_len / (double)dst_len;
  const double scale  = std::max(1.0, ratio);
  const double support = use_lanczos ? (LANCZOS_A * scale) : (1.0 * scale);
  const double inv_scale = 1.0 / scale;

  int total_weights = 0;
  for (int j = 0; j < dst_len; ++j) {
    const double centre = ((double)j + 0.5) * ratio - 0.5;
    int start = (int)std::floor(centre - support);
    int end   = (int)std::ceil(centre + support);
    if (start < 0) start = 0;
    if (end >= src_len) end = src_len - 1;
    const int count = end - start + 1;

    ft.entries[j].start = start;
    ft.entries[j].count = count;
    ft.offsets[j] = total_weights;
    total_weights += count;
  }

  ft.weights.resize(total_weights);

  for (int j = 0; j < dst_len; ++j) {
    const double centre = ((double)j + 0.5) * ratio - 0.5;
    const int start = ft.entries[j].start;
    const int count = ft.entries[j].count;
    double* w = ft.weights.data() + ft.offsets[j];

    double sum = 0.0;
    for (int k = 0; k < count; ++k) {
      const double x = ((double)(start + k) - centre) * inv_scale;
      double val = use_lanczos ? lanczos3(x) : std::max(0.0, 1.0 - std::fabs(x));
      w[k] = val;
      sum += val;
    }

    if (sum > 0.0) {
      const double inv_sum = 1.0 / sum;
      for (int k = 0; k < count; ++k) w[k] *= inv_sum;
    }
  }

  return ft;
}

template <typename T>
static void resample_horizontal(const T* __restrict__ src,
                                double*  __restrict__ dst,
                                int W, int H, int nch, int oW,
                                const FilterTable& ft) {
  const std::ptrdiff_t src_plane = (std::ptrdiff_t)W  * H;
  const std::ptrdiff_t dst_plane = (std::ptrdiff_t)oW * H;
  const int total_rows = nch * H;

#ifdef _OPENMP
  #pragma omp parallel for if((long long)total_rows * oW > 100000) schedule(static)
#endif
  for (int abs_r = 0; abs_r < total_rows; ++abs_r) {
    const int ch = abs_r / H;
    const int y  = abs_r % H;

    const T*      __restrict__ src_row = src + (std::ptrdiff_t)ch * src_plane + (std::ptrdiff_t)y * W;
    double*       __restrict__ dst_row = dst + (std::ptrdiff_t)ch * dst_plane + (std::ptrdiff_t)y * oW;

    for (int ox = 0; ox < oW; ++ox) {
      const FilterEntry& e = ft.entries[ox];
      const double* w = ft.weights.data() + ft.offsets[ox];
      const T* sp = src_row + e.start;

      double acc = 0.0;
      for (int k = 0; k < e.count; ++k) {
        acc += w[k] * (double)sp[k];
      }
      dst_row[ox] = acc;
    }
  }
}

template <typename T>
static void resample_vertical(const double* __restrict__ src,
                              T*            __restrict__ dst,
                              int oW, int H, int nch, int oH,
                              const FilterTable& ft) {
  const std::ptrdiff_t src_plane = (std::ptrdiff_t)oW * H;
  const std::ptrdiff_t dst_plane = (std::ptrdiff_t)oW * oH;
  const int total_cols = nch * oW;

#ifdef _OPENMP
  #pragma omp parallel for if((long long)total_cols * oH > 100000) schedule(static)
#endif
  for (int abs_c = 0; abs_c < total_cols; ++abs_c) {
    const int ch = abs_c / oW;
    const int x  = abs_c % oW;

    const double* __restrict__ src_base = src + (std::ptrdiff_t)ch * src_plane + x;
    T*            __restrict__ dst_base = dst + (std::ptrdiff_t)ch * dst_plane + x;

    for (int oy = 0; oy < oH; ++oy) {
      const FilterEntry& e = ft.entries[oy];
      const double* w = ft.weights.data() + ft.offsets[oy];

      double acc = 0.0;
      for (int k = 0; k < e.count; ++k) {
        acc += w[k] * src_base[(std::ptrdiff_t)(e.start + k) * oW];
      }
      
      if constexpr (std::is_same_v<T, Rbyte>) {
        if (acc < 0.0) acc = 0.0;
        else if (acc > 255.0) acc = 255.0;
        dst_base[(std::ptrdiff_t)oy * oW] = (T)std::round(acc);
      } else {
        if (acc < 0.0) acc = 0.0;
        else if (acc > 1.0) acc = 1.0;
        dst_base[(std::ptrdiff_t)oy * oW] = (T)acc;
      }
    }
  }
}

template <typename T>
static void resample_nearest(const T* __restrict__ src,
                             T*       __restrict__ dst,
                             int W, int H, int nch, int oW, int oH) {
  const std::ptrdiff_t src_plane = (std::ptrdiff_t)W  * H;
  const std::ptrdiff_t dst_plane = (std::ptrdiff_t)oW * oH;

  std::vector<int> x_map(oW);
  for (int ox = 0; ox < oW; ++ox) {
    x_map[ox] = std::min((int)((ox + 0.5) * W / oW), W - 1);
  }

  const int total_work = nch * oH;

#ifdef _OPENMP
  #pragma omp parallel for if((long long)total_work * oW > 100000) schedule(static)
#endif
  for (int abs_r = 0; abs_r < total_work; ++abs_r) {
    const int ch = abs_r / oH;
    const int oy = abs_r % oH;

    const int sy = std::min((int)((oy + 0.5) * H / oH), H - 1);
    const T* __restrict__ src_row = src + (std::ptrdiff_t)ch * src_plane + (std::ptrdiff_t)sy * W;
    T*       __restrict__ dst_row = dst + (std::ptrdiff_t)ch * dst_plane + (std::ptrdiff_t)oy * oW;

    for (int ox = 0; ox < oW; ++ox) {
      dst_row[ox] = src_row[x_map[ox]];
    }
  }
}

static inline bool is_identity(int W, int H, int oW, int oH) {
  return (W == oW && H == oH);
}

template <int RTYPE>
SEXP do_resize(SEXP img_sexp, int W, int H, int nch, int out_w, int out_h, int filter) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  
  VecType img(img_sexp);
  const T* src = img.begin();
  
  if (is_identity(W, H, out_w, out_h)) {
    return Rcpp::clone(img_sexp);
  }
  
  VecType out = Rcpp::no_init((std::size_t)out_w * out_h * nch);
  T* dst = out.begin();
  
  if (filter == 2) {
    resample_nearest<T>(src, dst, W, H, nch, out_w, out_h);
  } else {
    bool use_lanczos = (filter == 0);
    FilterTable ft_h = build_filter_table(W, out_w, use_lanczos);
    FilterTable ft_v = build_filter_table(H, out_h, use_lanczos);
    
    std::vector<double> tmp((std::size_t)out_w * H * nch);
    resample_horizontal<T>(src, tmp.data(), W, H, nch, out_w, ft_h);
    resample_vertical<T>(tmp.data(), dst, out_w, H, nch, out_h, ft_v);
  }
  
  if (nch == 1) {
    out.attr("dim") = IntegerVector::create(out_w, out_h);
  } else {
    out.attr("dim") = IntegerVector::create(out_w, out_h, nch);
  }
  return out;
}

// [[Rcpp::export]]
SEXP image_resize_cpp(SEXP img_sexp, int out_w, int out_h, int filter = 0) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);
  
  int W, H, nch;
  if (dims.size() == 2) {
    W = dims[0]; H = dims[1]; nch = 1;
  } else if (dims.size() == 3) {
    W = dims[0]; H = dims[1]; nch = dims[2];
  } else {
    Rcpp::stop("image_resize_cpp: input must be a 2-D matrix or 3-D array.");
    return R_NilValue;
  }

  if (out_w <= 0 || out_h <= 0) {
    Rcpp::stop("image_resize_cpp: out_w and out_h must be positive.");
    return R_NilValue;
  }

  if (TYPEOF(img_sexp) == RAWSXP) {
    return do_resize<RAWSXP>(img_sexp, W, H, nch, out_w, out_h, filter);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    return do_resize<REALSXP>(img_sexp, W, H, nch, out_w, out_h, filter);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    return do_resize<INTSXP>(img_sexp, W, H, nch, out_w, out_h, filter);
  } else {
    Rcpp::stop("image_resize_cpp: unsupported input type, must be raw, integer, or numeric.");
    return R_NilValue;
  }
}

// Computes only the window [x0, x0 + rw) x [y0, y0 + rh) (0-based) of
// image_resize_cpp(mat, out_w, out_h, filter = 1) for a single-channel double matrix.
// Uses the same filter tables and accumulation order, so the values are identical to
// cropping the full-resolution resize, at a cost proportional to the window size.
// Used to post-process SAM instance masks per box instead of per full image.
// [[Rcpp::export]]
Rcpp::NumericMatrix resize_roi_bilinear_cpp(Rcpp::NumericMatrix mat, int out_w, int out_h,
                                            int x0, int y0, int rw, int rh) {
  const int W = mat.nrow();
  const int H = mat.ncol();
  if (out_w <= 0 || out_h <= 0) Rcpp::stop("resize_roi_bilinear_cpp: invalid output size.");
  x0 = std::max(0, std::min(out_w - 1, x0));
  y0 = std::max(0, std::min(out_h - 1, y0));
  rw = std::max(1, std::min(rw, out_w - x0));
  rh = std::max(1, std::min(rh, out_h - y0));
  Rcpp::NumericMatrix out(rw, rh);
  const double* src = mat.begin();
  double* dst = out.begin();
  if (W == out_w && H == out_h) {
    for (int y = 0; y < rh; ++y)
      for (int x = 0; x < rw; ++x)
        dst[(size_t)y * rw + x] = src[(size_t)(y0 + y) * W + (x0 + x)];
    return out;
  }
  FilterTable ft_h = build_filter_table(W, out_w, false);
  FilterTable ft_v = build_filter_table(H, out_h, false);
  // Source rows needed by the vertical filter for the requested output rows.
  int r_min = H, r_max = -1;
  for (int oy = y0; oy < y0 + rh; ++oy) {
    const FilterEntry& e = ft_v.entries[oy];
    r_min = std::min(r_min, e.start);
    r_max = std::max(r_max, e.start + e.count - 1);
  }
  if (r_max < r_min) return out;
  const int nr = r_max - r_min + 1;
  std::vector<double> tmp((size_t)rw * nr);
  for (int r = 0; r < nr; ++r) {
    const double* src_row = src + (size_t)(r_min + r) * W;
    for (int x = 0; x < rw; ++x) {
      const FilterEntry& e = ft_h.entries[x0 + x];
      const double* w = ft_h.weights.data() + ft_h.offsets[x0 + x];
      const double* sp = src_row + e.start;
      double acc = 0.0;
      for (int k = 0; k < e.count; ++k) acc += w[k] * sp[k];
      tmp[(size_t)r * rw + x] = acc;
    }
  }
  for (int y = 0; y < rh; ++y) {
    const FilterEntry& e = ft_v.entries[y0 + y];
    const double* w = ft_v.weights.data() + ft_v.offsets[y0 + y];
    for (int x = 0; x < rw; ++x) {
      double acc = 0.0;
      for (int k = 0; k < e.count; ++k) acc += w[k] * tmp[(size_t)(e.start + k - r_min) * rw + x];
      // same [0, 1] clamp as resample_vertical() applies to double output
      if (acc < 0.0) acc = 0.0;
      else if (acc > 1.0) acc = 1.0;
      dst[(size_t)y * rw + x] = acc;
    }
  }
  return out;
}

// C++ version of the R helper `.preprocess_nchw()`: image array [W, H, C] (raw or
// double) -> normalised NCHW float tensor (returned as double) of size 3 x S x S.
// Matches the R implementation step by step (value scale detection, bilinear resize via
// the same filter tables with the same [0, 1] clamp, optional top-left letterbox padded
// with 0 before normalisation, (x - mean) / std) while avoiding the ~10 full-resolution
// double copies the R version allocates.
// [[Rcpp::export]]
Rcpp::NumericVector preprocess_nchw_cpp(SEXP mat, int target_size,
                                        Rcpp::NumericVector mean, Rcpp::NumericVector std_,
                                        bool letterbox = false) {
  SEXP dim_attr = Rf_getAttrib(mat, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 2) Rcpp::stop("preprocess_nchw_cpp: mat must be a matrix or array.");
  const int* d = INTEGER(dim_attr);
  const int W = d[0], H = d[1];
  const int nch = (Rf_length(dim_attr) >= 3) ? d[2] : 1;
  const int S = target_size;
  if (S <= 0) Rcpp::stop("preprocess_nchw_cpp: invalid target_size.");
  if (mean.size() < 3 || std_.size() < 3) Rcpp::stop("preprocess_nchw_cpp: mean and std must have length 3.");
  const bool is_raw = (TYPEOF(mat) == RAWSXP);
  if (!is_raw && TYPEOF(mat) != REALSXP && TYPEOF(mat) != INTSXP) Rcpp::stop("preprocess_nchw_cpp: unsupported storage.");

  const size_t plane = (size_t)W * H;
  double val_scale = 1.0;
  if (is_raw) {
    val_scale = 1.0 / 255.0;
  } else {
    size_t nchk = std::min((size_t)1000, (size_t)Rf_xlength(mat));
    double mx = -INFINITY;
    for (size_t i = 0; i < nchk; ++i) {
      double v = (TYPEOF(mat) == REALSXP) ? REAL(mat)[i] : (double)INTEGER(mat)[i];
      if (!ISNAN(v) && v > mx) mx = v;
    }
    val_scale = (mx > 1.5) ? (1.0 / 255.0) : 1.0;
  }

  int ow = S, oh = S;
  if (letterbox) {
    double sc = (double)S / (double)std::max(W, H);
    ow = std::max(1, (int)std::nearbyint(W * sc));
    oh = std::max(1, (int)std::nearbyint(H * sc));
  }
  const bool identity = (ow == W && oh == H);

  Rcpp::NumericVector out((R_xlen_t)3 * S * S);
  double* out_ptr = out.begin();
  std::vector<double> src_plane(plane);
  std::vector<double> tmp;
  std::vector<double> res;
  FilterTable ft_h, ft_v;
  if (!identity) {
    ft_h = build_filter_table(W, ow, false);
    ft_v = build_filter_table(H, oh, false);
    tmp.resize((size_t)ow * H);
    res.resize((size_t)ow * oh);
  }

  for (int c = 0; c < 3; ++c) {
    const int sc_idx = (nch >= 3) ? c : 0;
    const size_t off = (size_t)sc_idx * plane;
    if (is_raw) {
      const Rbyte* p = RAW(mat) + off;
      for (size_t i = 0; i < plane; ++i) src_plane[i] = (double)p[i] * val_scale;
    } else if (TYPEOF(mat) == REALSXP) {
      const double* p = REAL(mat) + off;
      for (size_t i = 0; i < plane; ++i) src_plane[i] = p[i] * val_scale;
    } else {
      const int* p = INTEGER(mat) + off;
      for (size_t i = 0; i < plane; ++i) src_plane[i] = (double)p[i] * val_scale;
    }
    const double* resized = src_plane.data();
    if (!identity) {
      resample_horizontal<double>(src_plane.data(), tmp.data(), W, H, 1, ow, ft_h);
      resample_vertical<double>(tmp.data(), res.data(), ow, H, 1, oh, ft_v);
      resized = res.data();
    }
    const double mu = mean[c], sd = std_[c];
    double* dst = out_ptr + (size_t)c * S * S;
    if (letterbox) {
      const double pad_val = (0.0 - mu) / sd;
      for (size_t i = 0; i < (size_t)S * S; ++i) dst[i] = pad_val;
      for (int y = 0; y < oh; ++y)
        for (int x = 0; x < ow; ++x)
          dst[(size_t)y * S + x] = (resized[(size_t)y * ow + x] - mu) / sd;
    } else {
      for (size_t i = 0; i < (size_t)S * S; ++i) dst[i] = (resized[i] - mu) / sd;
    }
  }
  out.attr("dims") = Rcpp::IntegerVector::create(1, 3, S, S);
  if (letterbox) out.attr("letterbox_dims") = Rcpp::NumericVector::create(ow, oh);
  return out;
}
