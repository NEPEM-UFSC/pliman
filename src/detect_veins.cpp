// detect_veins.cpp — Advanced resolution-independent leaf vein detection
//
// Improvements:
//   1. Dynamic Resolution-Independent Boundary Erosion (BFS Distance Transform):
//      - The outer leaf boundary causes massive DoG gradient responses.
//      - To prevent detecting the leaf edge as a vein, we dynamically compute
//        an erosion radius r = round(min(W, H) * rel_erode) (default 1.0% of min(W,H))
//        or a specified erode_size in pixels.
//      - Multi-source BFS distance transform trims the outer r pixels of each object mask,
//        completely eliminating boundary edge artifacts without hardcoding pixel values.
//
//   2. Guo-Hall Skeleton Thinning Support (thinning = TRUE):
//      - Option to apply Guo-Hall thinning on detected vein pixels.
//      - Thinning converts veins to 1-pixel-wide skeletons.
//      - Ideal for measuring total vein length density (skeletal length per unit leaf area)
//        independent of vein thickness.
//      - Without thinning (thinning = FALSE), measures vein area fraction (% of leaf area).
//
//   3. DoG + Per-Object Adaptive Otsu Thresholding:
//      - Difference of Gaussians (DoG) bandpass filter responds selectively to vein scales.
//      - Otsu thresholding performed per object on inner leaf pixels.

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <queue>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// ─────────────────────────────────────────────────────────────────────────────
// 1D Gaussian kernel
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double> make_gaussian_kernel(double sigma) {
  int r = std::max(1, (int)std::ceil(3.0 * sigma));
  int ksize = 2 * r + 1;
  std::vector<double> k(ksize);
  double s = 0.0;
  for (int i = 0; i < ksize; ++i) {
    double x = i - r;
    k[i] = std::exp(-x * x / (2.0 * sigma * sigma));
    s += k[i];
  }
  for (int i = 0; i < ksize; ++i) k[i] /= s;
  return k;
}

// ─────────────────────────────────────────────────────────────────────────────
// Separable 2D Gaussian blur (column-major, reflect-101 padding)
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double>
gaussian_blur_cm(const double* src, int W, int H, double sigma) {
  const std::vector<double> kernel = make_gaussian_kernel(sigma);
  const int r = (int)(kernel.size() / 2);

  std::vector<double> tmp(W * H), out(W * H);

  for (int j = 0; j < H; ++j) {
    const double* row = src + (std::ptrdiff_t)j * W;
    double* dst_row = tmp.data() + (std::ptrdiff_t)j * W;
    for (int i = 0; i < W; ++i) {
      double acc = 0.0;
      for (int k = -r; k <= r; ++k) {
        int ii = i + k;
        if (ii < 0) ii = -ii;
        if (ii >= W) ii = 2 * W - 2 - ii;
        acc += kernel[k + r] * row[ii];
      }
      dst_row[i] = acc;
    }
  }

  for (int i = 0; i < W; ++i) {
    for (int j = 0; j < H; ++j) {
      double acc = 0.0;
      for (int k = -r; k <= r; ++k) {
        int jj = j + k;
        if (jj < 0) jj = -jj;
        if (jj >= H) jj = 2 * H - 2 - jj;
        acc += kernel[k + r] * tmp[i + (std::ptrdiff_t)jj * W];
      }
      out[i + (std::ptrdiff_t)j * W] = acc;
    }
  }

  return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-object Otsu thresholding
// ─────────────────────────────────────────────────────────────────────────────
static double otsu_threshold(const double* vals, int n, int n_bins = 256) {
  if (n == 0) return 0.0;
  double vmin = vals[0], vmax = vals[0];
  for (int i = 1; i < n; ++i) {
    if (vals[i] < vmin) vmin = vals[i];
    if (vals[i] > vmax) vmax = vals[i];
  }
  if (vmax == vmin) return vmin;

  double range = vmax - vmin;
  double scale = (n_bins - 1) / range;

  std::vector<double> hist(n_bins, 0.0);
  for (int i = 0; i < n; ++i)
    hist[(int)((vals[i] - vmin) * scale)]++;

  double total   = n;
  double sum_all = 0.0;
  for (int i = 0; i < n_bins; ++i) sum_all += i * hist[i];

  double sum_bg = 0.0, w_bg = 0.0, best_var = -1.0;
  int best_t = 0;
  for (int t = 0; t < n_bins; ++t) {
    w_bg += hist[t];
    if (w_bg == 0.0) continue;
    double w_fg = total - w_bg;
    if (w_fg == 0.0) break;
    sum_bg += t * hist[t];
    double mb = sum_bg / w_bg;
    double mf = (sum_all - sum_bg) / w_fg;
    double var = w_bg * w_fg * (mb - mf) * (mb - mf);
    if (var > best_var) { best_var = var; best_t = t; }
  }

  return vmin + best_t / scale;
}

// ─────────────────────────────────────────────────────────────────────────────
// Fast Guo-Hall thinning algorithm on binary column-major map
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<int> guo_hall_thinning(const std::vector<int>& binary_map, int W, int H) {
  std::vector<int> grid = binary_map;

  bool changed = true;
  while (changed) {
    changed = false;
    for (int sub = 0; sub < 2; ++sub) {
      bool even = (sub == 1);
      std::vector<int> to_remove;

      for (int j = 1; j < H - 1; ++j) {
        for (int i = 1; i < W - 1; ++i) {
          int k = i + j * W;
          if (grid[k] == 0) continue;

          int p2 = grid[(i - 1) + j * W];
          int p3 = grid[(i - 1) + (j + 1) * W];
          int p4 = grid[i + (j + 1) * W];
          int p5 = grid[(i + 1) + (j + 1) * W];
          int p6 = grid[(i + 1) + j * W];
          int p7 = grid[(i + 1) + (j - 1) * W];
          int p8 = grid[i + (j - 1) * W];
          int p9 = grid[(i - 1) + (j - 1) * W];

          int C = ((!p2) & (p3 | p4)) + ((!p4) & (p5 | p6)) + ((!p6) & (p7 | p8)) + ((!p8) & (p9 | p2));
          if (C != 1) continue;

          int N1 = (p9 | p2) + (p3 | p4) + (p5 | p6) + (p7 | p8);
          int N2 = (p2 | p3) + (p4 | p5) + (p6 | p7) + (p8 | p9);
          int N = (N1 < N2) ? N1 : N2;
          if (N < 2 || N > 3) continue;

          int m = even ? ((p6 | p7 | (!p9)) & p8) : ((p2 | p3 | (!p5)) & p4);
          if (m == 0) {
            to_remove.push_back(k);
          }
        }
      }

      if (!to_remove.empty()) {
        changed = true;
        for (int idx : to_remove) grid[idx] = 0;
      }
    }
  }

  return grid;
}

// ─────────────────────────────────────────────────────────────────────────────
// detect_veins_cpp — main exported function
//
// Parameters:
//   R, G, B      : NumericMatrix (W × H), pixel intensities [0, 1]
//   labels_sexp  : W × H integer label matrix (0 = background)
//   sigma1       : fine Gaussian scale (default 0.75)
//   sigma2       : coarse Gaussian scale (default 3.75)
//   threshold    : fixed vein threshold; if < 0 → per-object Otsu adaptive
//   channel      : 0 = green (default), 1 = grayscale luminance
//   erode_size   : fixed boundary erosion radius in pixels; if < 0 uses rel_erode
//   rel_erode    : relative boundary erosion fraction of min(W,H) (default 0.01 = 1%)
//   thinning     : if TRUE, applies Guo-Hall thinning to measure skeleton length density
//   return_map   : if TRUE, returns W×H binary vein map
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List detect_veins_cpp(SEXP R_sexp,
                      SEXP G_sexp,
                      SEXP B_sexp,
                      SEXP   labels_sexp,
                      double sigma1      = 0.75,
                      double sigma2      = 3.75,
                      double threshold   = -1.0,
                      int    channel     = 0,
                      int    erode_size  = -1,
                      double rel_erode   = 0.01,
                      bool   thinning    = false,
                      bool   return_map  = false) {

  // ── Coerce labels ──────────────────────────────────────────────────────────
  IntegerMatrix labels;
  switch (TYPEOF(labels_sexp)) {
    case INTSXP:  labels = as<IntegerMatrix>(labels_sexp); break;
    case REALSXP: labels = as<IntegerMatrix>(NumericMatrix(labels_sexp)); break;
    default: Rcpp::stop("detect_veins_cpp: 'labels' must be integer or numeric.");
  }

  IntegerVector dims = Rf_getAttrib(R_sexp, R_DimSymbol);
  const int W = dims[0], H = dims[1];
  const int N = W * H;
  
  IntegerVector dims_g = Rf_getAttrib(G_sexp, R_DimSymbol);
  IntegerVector dims_b = Rf_getAttrib(B_sexp, R_DimSymbol);

  if (dims_g[0] != W || dims_g[1] != H || dims_b[0] != W || dims_b[1] != H ||
      labels.nrow() != W || labels.ncol() != H) {
    Rcpp::stop("detect_veins_cpp: all matrices must have the same dimensions.");
  }
  if (sigma1 >= sigma2)
    Rcpp::stop("detect_veins_cpp: sigma1 must be less than sigma2.");

  const int* lp = labels.begin();

  // ── Step 1: Dynamic Resolution-Independent Boundary Erosion ────────────────
  // Calculate erosion radius r
  int r = erode_size;
  if (r < 0) {
    r = (int)std::round(std::min(W, H) * rel_erode);
  }
  if (r < 1) r = 1;

  // Multi-source BFS to compute distance from object boundaries
  std::vector<int> dist(N, 999999);
  std::queue<int> q;

  for (int j = 0; j < H; ++j) {
    for (int i = 0; i < W; ++i) {
      int k = i + j * W;
      int lab = lp[k];
      if (lab == 0) continue;

      bool is_edge = (i == 0 || i == W - 1 || j == 0 || j == H - 1);
      if (!is_edge) {
        if (lp[(i - 1) + j * W] != lab || lp[(i + 1) + j * W] != lab ||
            lp[i + (j - 1) * W] != lab || lp[i + (j + 1) * W] != lab) {
          is_edge = true;
        }
      }
      if (is_edge) {
        dist[k] = 1;
        q.push(k);
      }
    }
  }

  const int dx[4] = {-1, 1, 0, 0};
  const int dy[4] = {0, 0, -1, 1};
  while (!q.empty()) {
    int curr = q.front();
    q.pop();
    int d = dist[curr];
    if (d >= r) continue;

    int cx = curr % W;
    int cy = curr / W;
    for (int dir = 0; dir < 4; ++dir) {
      int nx = cx + dx[dir];
      int ny = cy + dy[dir];
      if (nx >= 0 && nx < W && ny >= 0 && ny < H) {
        int nk = nx + ny * W;
        if (lp[nk] > 0 && dist[nk] > d + 1) {
          dist[nk] = d + 1;
          q.push(nk);
        }
      }
    }
  }

  // Valid inner pixels: lp[k] > 0 AND dist[k] > r (excluding border artifacts)

  // ── Step 2: Build vein channel ─────────────────────────────────────────────
  std::vector<double> ch_vec(N);
  
  if (TYPEOF(R_sexp) == RAWSXP) {
    Rbyte* rp = RAW(R_sexp);
    Rbyte* gp = RAW(G_sexp);
    Rbyte* bp = RAW(B_sexp);
    if (channel == 0) {
      for (int k = 0; k < N; ++k) ch_vec[k] = gp[k] / 255.0;
    } else {
      for (int k = 0; k < N; ++k)
        ch_vec[k] = (0.299 * rp[k] + 0.587 * gp[k] + 0.114 * bp[k]) / 255.0;
    }
  } else if (TYPEOF(R_sexp) == REALSXP) {
    double* rp = REAL(R_sexp);
    double* gp = REAL(G_sexp);
    double* bp = REAL(B_sexp);
    if (channel == 0) {
      for (int k = 0; k < N; ++k) ch_vec[k] = gp[k];
    } else {
      for (int k = 0; k < N; ++k)
        ch_vec[k] = 0.299 * rp[k] + 0.587 * gp[k] + 0.114 * bp[k];
    }
  } else if (TYPEOF(R_sexp) == INTSXP) {
    int* rp = INTEGER(R_sexp);
    int* gp = INTEGER(G_sexp);
    int* bp = INTEGER(B_sexp);
    if (channel == 0) {
      for (int k = 0; k < N; ++k) ch_vec[k] = gp[k] / 255.0;
    } else {
      for (int k = 0; k < N; ++k)
        ch_vec[k] = (0.299 * rp[k] + 0.587 * gp[k] + 0.114 * bp[k]) / 255.0;
    }
  } else {
    Rcpp::stop("detect_veins_cpp: unsupported image type.");
  }


  // ── Step 3: DoG Bandpass Filter ───────────────────────────────────────────
  std::vector<double> blur1 = gaussian_blur_cm(ch_vec.data(), W, H, sigma1);
  std::vector<double> blur2 = gaussian_blur_cm(ch_vec.data(), W, H, sigma2);

  std::vector<double> dog(N);
  for (int k = 0; k < N; ++k)
    dog[k] = std::fabs(blur1[k] - blur2[k]);

  // ── Step 4: Collect per-object inner pixel DoG values ─────────────────────
  int max_label = 0;
  for (int k = 0; k < N; ++k)
    if (lp[k] > max_label) max_label = lp[k];

  const int nlab = max_label + 1;

  std::vector<int> inner_cnt(nlab, 0);
  for (int k = 0; k < N; ++k) {
    if (lp[k] > 0 && dist[k] > r) {
      ++inner_cnt[lp[k]];
    }
  }

  std::vector<std::vector<double>> obj_dog(nlab);
  for (int lab = 1; lab < nlab; ++lab)
    obj_dog[lab].reserve(inner_cnt[lab]);

  for (int j = 0; j < H; ++j) {
    const int base = j * W;
    for (int i = 0; i < W; ++i) {
      const int k   = base + i;
      const int lab = lp[k];
      if (lab > 0 && dist[k] > r) {
        obj_dog[lab].push_back(dog[k]);
      }
    }
  }

  // ── Step 5: Per-object Otsu / Thresholding & Binary Vein Map ───────────────
  std::vector<double> thr_vec(nlab, 0.0);
  for (int lab = 1; lab < nlab; ++lab) {
    const auto& dv = obj_dog[lab];
    if (dv.empty()) continue;
    if (threshold < 0.0) {
      thr_vec[lab] = otsu_threshold(dv.data(), (int)dv.size());
    } else {
      thr_vec[lab] = threshold;
    }
  }

  std::vector<int> binary_vein_map(N, 0);
  for (int k = 0; k < N; ++k) {
    const int lab = lp[k];
    if (lab > 0 && dist[k] > r) {
      if (dog[k] > thr_vec[lab]) {
        binary_vein_map[k] = 1;
      }
    }
  }

  // ── Step 6: Optional Thinning (Guo-Hall) ──────────────────────────────────
  std::vector<int> final_vein_map;
  if (thinning) {
    final_vein_map = guo_hall_thinning(binary_vein_map, W, H);
  } else {
    final_vein_map = binary_vein_map;
  }

  // ── Step 7: Calculate per-object proportions ──────────────────────────────
  NumericVector proportion(nlab, NA_REAL);
  for (int lab = 1; lab < nlab; ++lab) {
    if (inner_cnt[lab] == 0) {
      proportion[lab] = 0.0;
      continue;
    }
    int n_vein = 0;
    // Count vein / thinned pixels belonging to object lab
    for (int j = 0; j < H; ++j) {
      const int base = j * W;
      for (int i = 0; i < W; ++i) {
        const int k = base + i;
        if (lp[k] == lab && dist[k] > r) {
          if (final_vein_map[k] != 0) ++n_vein;
        }
      }
    }
    proportion[lab] = (double)n_vein / (double)inner_cnt[lab];
  }

  // ── Return ─────────────────────────────────────────────────────────────────
  if (return_map) {
    IntegerMatrix vmap(W, H);
    std::copy(final_vein_map.begin(), final_vein_map.end(), vmap.begin());
    return List::create(Named("proportion") = proportion,
                         Named("vein_map")   = vmap);
  }

  return List::create(Named("proportion") = proportion,
                       Named("vein_map")   = R_NilValue);
}
