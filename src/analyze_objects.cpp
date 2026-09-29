#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Forward declarations for functions defined in other C++ translation units
SEXP compute_single_index_cpp(SEXP img_sexp, std::string ind, int r = 1, int g = 2, int b = 3, int re = 4, int nir = 5, int swir = 6, std::string storage = "auto");
double help_otsu(SEXP img_sexp);
SEXP cpp_binary_threshold(SEXP img_sexp, double threshold, int op = 1, bool return_raw = true);
Rcpp::NumericMatrix threshold_adaptive(Rcpp::NumericMatrix mat, double k, int windowsize, double maxsd);
LogicalMatrix help_binary_filters_cpp(SEXP img, int erode = 0, int dilate = 0, int opening = 0, int closing = 0, int filter = 0, bool fill_hull = false, CharacterVector filter_order = CharacterVector(), double max_size = -1.0, double min_neck_dist = 0.0);
LogicalMatrix erode_cpp(LogicalMatrix img, int raio = 3, std::string forma = "disc");
IntegerMatrix bwlabel_cpp(SEXP img_sexp);
IntegerMatrix watershed_cpp(SEXP img_r, double tolerance = 1.0, int ext = 1);
List extract_contours_cpp(IntegerMatrix labels);
List smoothContours(List contours, int window_size = 3);
DataFrame poly_measures_cpp(List contours, bool calc_pcv = false);
List compute_chulls_cpp(List contours);
NumericMatrix haralick_features_cpp(IntegerMatrix labels, SEXP ref_sexp, int nc = 32);

// [[Rcpp::export]]
void filter_labels_cpp(IntegerMatrix labels, IntegerVector keep_ids) {
  int N = labels.nrow() * labels.ncol();
  int* p = INTEGER(labels);
  int max_id = 0;
  #pragma omp parallel for reduction(max:max_id) schedule(static)
  for (int i = 0; i < N; i++) {
    if (p[i] > max_id) max_id = p[i];
  }
  if (max_id <= 0) return;

  std::vector<uint8_t> keep(max_id + 1, 0);
  for (int i = 0; i < keep_ids.size(); i++) {
    int id = keep_ids[i];
    if (id >= 0 && id <= max_id) keep[id] = 1;
  }

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    int id = p[i];
    if (id > 0 && (id > max_id || !keep[id])) {
      p[i] = 0;
    }
  }
}

// Otsu threshold computation on a 256-bin histogram
static inline int otsu_from_histogram(const int* hist, int totalPixels) {
  double sum = 0;
  for (int i = 0; i < 256; i++) sum += (double)i * hist[i];
  double sumBackground = 0;
  int backgroundPixels = 0;
  double maxVariance = 0;
  int threshold = 0;

  for (int i = 0; i < 256; i++) {
    backgroundPixels += hist[i];
    if (backgroundPixels == 0) continue;
    int foregroundPixels = totalPixels - backgroundPixels;
    if (foregroundPixels == 0) break;

    sumBackground += (double)i * hist[i];
    double meanBackground = sumBackground / backgroundPixels;
    double meanForeground = (sum - sumBackground) / foregroundPixels;
    double variance = (double)backgroundPixels * (double)foregroundPixels * (meanBackground - meanForeground) * (meanBackground - meanForeground);

    if (variance > maxVariance) {
      maxVariance = variance;
      threshold = i;
    }
  }
  return threshold;
}

// [[Rcpp::export]]
List analyze_objects_cpp(SEXP img_sexp,
                         SEXP bin_sexp = R_NilValue,
                         std::string index_str = "NB",
                         int r = 1, int g = 2, int b = 3,
                         int re = 4, int nir = 5, int swir = 6,
                         std::string threshold_method = "Otsu",
                         double threshold_val = 0.5,
                         double k_adj = 0.1,
                         int windowsize = 15,
                         bool invert = false,
                         int erode_sz = 0,
                         int dilate_sz = 0,
                         int opening_sz = 0,
                         int closing_sz = 0,
                         int filter_sz = 0,
                         bool fill_hull = false,
                         CharacterVector filter_order = CharacterVector(),
                         bool return_exact = false,
                         bool watershed = true,
                         double tolerance = 1.0,
                         int ext = 1,
                         bool haralick = false,
                         int har_nbins = 32,
                         int har_band = 1,
                         int smooth = 0) {

  // Step 1: Obtain Binary Mask
  LogicalMatrix bin_mat;
  int nrow = 0, ncol = 0, N = 0;

  if (!Rf_isNull(bin_sexp)) {
    bin_mat = as<LogicalMatrix>(bin_sexp);
    nrow = bin_mat.nrow();
    ncol = bin_mat.ncol();
    N = nrow * ncol;
  } else {
    SEXP dims = Rf_getAttrib(img_sexp, R_DimSymbol);
    nrow = INTEGER(dims)[0];
    ncol = INTEGER(dims)[1];
    N = nrow * ncol;

    // Ultra-fast fused path for Otsu on raw uint8 RGB
    if (TYPEOF(img_sexp) == RAWSXP && threshold_method == "Otsu" && (index_str == "NB" || index_str == "NR" || index_str == "NG" || index_str == "GRAY" || index_str == "R" || index_str == "G" || index_str == "B")) {
      bin_mat = LogicalMatrix(nrow, ncol);
      int* p_out = LOGICAL(bin_mat);
      const uint8_t* ptr = RAW(img_sexp);
      const uint8_t* pR = ptr + (r - 1) * N;
      const uint8_t* pG = ptr + (g - 1) * N;
      const uint8_t* pB = ptr + (b - 1) * N;

      int hist[256] = {0};

      if (index_str == "NB") {
        #pragma omp parallel
        {
          int lhist[256] = {0};
          #pragma omp for schedule(static)
          for (int i = 0; i < N; i++) {
            int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
            int bin = (s > 0) ? (((int)pB[i] * 255) / s) : 0;
            lhist[bin]++;
          }
          #pragma omp critical
          {
            for (int k = 0; k < 256; k++) hist[k] += lhist[k];
          }
        }
        int t = otsu_from_histogram(hist, N);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
          int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
          int p255 = (int)pB[i] * 255;
          int ts = t * s;
          p_out[i] = invert ? (s > 0 && p255 > ts) : (s > 0 && p255 < ts);
        }
      } else if (index_str == "NR") {
        #pragma omp parallel
        {
          int lhist[256] = {0};
          #pragma omp for schedule(static)
          for (int i = 0; i < N; i++) {
            int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
            int bin = (s > 0) ? (((int)pR[i] * 255) / s) : 0;
            lhist[bin]++;
          }
          #pragma omp critical
          {
            for (int k = 0; k < 256; k++) hist[k] += lhist[k];
          }
        }
        int t = otsu_from_histogram(hist, N);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
          int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
          int p255 = (int)pR[i] * 255;
          int ts = t * s;
          p_out[i] = invert ? (s > 0 && p255 > ts) : (s > 0 && p255 < ts);
        }
      } else if (index_str == "NG") {
        #pragma omp parallel
        {
          int lhist[256] = {0};
          #pragma omp for schedule(static)
          for (int i = 0; i < N; i++) {
            int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
            int bin = (s > 0) ? (((int)pG[i] * 255) / s) : 0;
            lhist[bin]++;
          }
          #pragma omp critical
          {
            for (int k = 0; k < 256; k++) hist[k] += lhist[k];
          }
        }
        int t = otsu_from_histogram(hist, N);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
          int s = (int)pR[i] + (int)pG[i] + (int)pB[i];
          int p255 = (int)pG[i] * 255;
          int ts = t * s;
          p_out[i] = invert ? (s > 0 && p255 > ts) : (s > 0 && p255 < ts);
        }
      } else if (index_str == "GRAY") {
        #pragma omp parallel
        {
          int lhist[256] = {0};
          #pragma omp for schedule(static)
          for (int i = 0; i < N; i++) {
            int v = (299 * (int)pR[i] + 587 * (int)pG[i] + 114 * (int)pB[i]) / 1000;
            if (v < 0) v = 0; else if (v > 255) v = 255;
            lhist[v]++;
          }
          #pragma omp critical
          {
            for (int k = 0; k < 256; k++) hist[k] += lhist[k];
          }
        }
        int t = otsu_from_histogram(hist, N);
        int t1000 = t * 1000;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
          int v = 299 * (int)pR[i] + 587 * (int)pG[i] + 114 * (int)pB[i];
          p_out[i] = invert ? (v > t1000) : (v < t1000);
        }
      } else {
        const uint8_t* pChan = (index_str == "R") ? pR : ((index_str == "G") ? pG : pB);
        #pragma omp parallel
        {
          int lhist[256] = {0};
          #pragma omp for schedule(static)
          for (int i = 0; i < N; i++) lhist[pChan[i]]++;
          #pragma omp critical
          {
            for (int k = 0; k < 256; k++) hist[k] += lhist[k];
          }
        }
        int t = otsu_from_histogram(hist, N);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) p_out[i] = invert ? (pChan[i] > t) : (pChan[i] < t);
      }
    } else {
      SEXP idx_sexp = PROTECT(compute_single_index_cpp(img_sexp, index_str, r, g, b, re, nir, swir, "auto"));
      if (threshold_method == "Otsu") {
        double otsu_t = help_otsu(idx_sexp);
        SEXP bin_sexp_tmp = PROTECT(cpp_binary_threshold(idx_sexp, otsu_t, invert ? 3 : 1, false));
        bin_mat = as<LogicalMatrix>(bin_sexp_tmp);
        UNPROTECT(2);
      } else if (threshold_method == "adaptive" || threshold_method == "Adaptive") {
        NumericMatrix idx_num = as<NumericMatrix>(idx_sexp);
        int wsize = windowsize;
        if (wsize <= 2) {
          wsize = std::min(idx_num.nrow(), idx_num.ncol()) / 3;
          if (wsize % 2 == 0) wsize++;
        }
        if (wsize < 3) wsize = 3;
        NumericMatrix ad_bin = threshold_adaptive(idx_num, k_adj, wsize, 1.0);
        bin_mat = LogicalMatrix(ad_bin.nrow(), ad_bin.ncol());
        const double* p_ad = REAL(ad_bin);
        int* p_bin = LOGICAL(bin_mat);
        for (int i = 0; i < N; i++) {
          p_bin[i] = invert ? (p_ad[i] == 0.0) : (p_ad[i] != 0.0);
        }
        UNPROTECT(1);
      } else {
        SEXP bin_sexp_tmp = PROTECT(cpp_binary_threshold(idx_sexp, threshold_val, invert ? 3 : 1, false));
        bin_mat = as<LogicalMatrix>(bin_sexp_tmp);
        UNPROTECT(2);
      }
    }
  }

  // Step 2: Apply Morphological Filters
  if (erode_sz > 0 || dilate_sz > 0 || opening_sz > 0 || closing_sz > 0 || filter_sz > 0 || fill_hull) {
    if (filter_order.size() == 0) {
      filter_order = CharacterVector::create("erode", "dilate", "opening", "closing", "filter", "fill_hull");
    }

    bool do_exact_opening = return_exact && (opening_sz > 0);

    bin_mat = help_binary_filters_cpp(
      bin_mat,
      erode_sz,
      dilate_sz,
      do_exact_opening ? 0 : opening_sz,
      closing_sz,
      filter_sz,
      fill_hull,
      filter_order
    );

    if (do_exact_opening) {
      LogicalMatrix eroded_seed = erode_cpp(bin_mat, opening_sz, "disc");
      IntegerMatrix temp_labels = bwlabel_cpp(bin_mat);

      if (N > 0 && temp_labels.size() == N && eroded_seed.size() == N) {
        const int* p_lab = INTEGER(temp_labels);
        const int* p_seed = LOGICAL(eroded_seed);

        int max_lab = 0;
        for (int i = 0; i < N; i++) {
          if (p_lab[i] > max_lab) max_lab = p_lab[i];
        }

        if (max_lab > 0) {
          std::vector<uint8_t> surviving(max_lab + 1, 0);
          for (int i = 0; i < N; i++) {
            if (p_seed[i] != 0 && p_seed[i] != NA_LOGICAL && p_lab[i] > 0 && p_lab[i] <= max_lab) {
              surviving[p_lab[i]] = 1;
            }
          }

          int* p_bin = LOGICAL(bin_mat);
          for (int i = 0; i < N; i++) {
            if (p_lab[i] > 0 && p_lab[i] <= max_lab && !surviving[p_lab[i]]) {
              p_bin[i] = 0;
            }
          }
        }
      }
    }
  }

  // Step 3: Labeling (Watershed vs BWLabel)
  IntegerMatrix labels;
  if (watershed) {
    labels = watershed_cpp(bin_mat, tolerance, ext);
  } else {
    labels = bwlabel_cpp(bin_mat);
  }

  int* p_labels = INTEGER(labels);

  // Step 4: Extract Contours & Compute Shape Measures & Convex Hulls
  int max_label = 0;
  #pragma omp parallel for reduction(max:max_label) schedule(static)
  for (int i = 0; i < N; i++) {
    if (p_labels[i] > max_label) max_label = p_labels[i];
  }

  if (max_label == 0) {
    return List::create(
      _["labels"] = labels,
      _["shape"] = DataFrame(),
      _["contours"] = List(),
      _["chull"] = List(),
      _["haralick"] = R_NilValue
    );
  }

  List ocont;
  DataFrame shape;
  List ch_list;

  if (smooth > 0) {
    ocont = extract_contours_cpp(labels);
    ocont = smoothContours(ocont, smooth);
    shape = poly_measures_cpp(ocont, false);
    ch_list = compute_chulls_cpp(ocont);
  } else {
    std::vector<int> start_r(max_label + 1, -1);
    std::vector<int> start_c(max_label + 1, -1);
    std::vector<int> start_idx(max_label + 1, -1);

    int idx = 0;
    for (int c = 0; c < ncol; c++) {
      for (int r = 0; r < nrow; r++, idx++) {
        int id = p_labels[idx];
        if (id > 0 && start_idx[id] == -1) {
          start_r[id] = r;
          start_c[id] = c;
          start_idx[id] = idx;
        }
      }
    }

    const int dr[] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dc[] = { 0,  1,  1,  1, 0,-1,-1, -1};
    const int doff[] = {-1, nrow - 1, nrow, nrow + 1, 1, -nrow + 1, -nrow, -nrow - 1};

    std::vector<std::vector<int>> all_bx(max_label + 1);
    std::vector<std::vector<int>> all_by(max_label + 1);

    for (int id = 1; id <= max_label; id++) {
      if (start_idx[id] == -1) continue;

      int sr = start_r[id], sc = start_c[id], sidx = start_idx[id];
      auto& b_x = all_bx[id];
      auto& b_y = all_by[id];
      b_x.reserve(256);
      b_y.reserve(256);

      b_x.push_back(sr);
      b_y.push_back(sc);

      int curr_r = sr, curr_c = sc, curr_idx = sidx;
      int backtrack = 6;
      int next_r = -1, next_c = -1, next_idx = -1, next_backtrack = -1;
      bool found = false;

      for (int i = 1; i <= 8; i++) {
        int dir = backtrack + i;
        if (dir >= 8) dir -= 8;
        int nr = curr_r + dr[dir], nc = curr_c + dc[dir];
        if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
            static_cast<unsigned>(nc) < static_cast<unsigned>(ncol)) {
          int nidx = curr_idx + doff[dir];
          if (p_labels[nidx] == id) {
            next_r = nr; next_c = nc; next_idx = nidx;
            next_backtrack = dir + 4;
            if (next_backtrack >= 8) next_backtrack -= 8;
            found = true;
            break;
          }
        }
      }

      if (found && (next_idx != sidx)) {
        int second_idx = next_idx;
        curr_r = next_r; curr_c = next_c; curr_idx = next_idx; backtrack = next_backtrack;
        b_x.push_back(curr_r); b_y.push_back(curr_c);
        int max_iter = N;
        int iter = 0;

        while (iter++ < max_iter) {
          bool step_found = false;
          for (int i = 1; i <= 8; i++) {
            int dir = backtrack + i;
            if (dir >= 8) dir -= 8;
            int nr = curr_r + dr[dir], nc = curr_c + dc[dir];
            if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
                static_cast<unsigned>(nc) < static_cast<unsigned>(ncol)) {
              int nidx = curr_idx + doff[dir];
              if (p_labels[nidx] == id) {
                next_r = nr; next_c = nc; next_idx = nidx;
                next_backtrack = dir + 4;
                if (next_backtrack >= 8) next_backtrack -= 8;
                step_found = true;
                break;
              }
            }
          }
          if (!step_found) break;
          if (curr_idx == sidx && next_idx == second_idx) break;
          curr_r = next_r; curr_c = next_c; curr_idx = next_idx; backtrack = next_backtrack;
          b_x.push_back(curr_r); b_y.push_back(curr_c);
        }
      }
    }

    struct ObjMetric {
      double mass_x = NA_REAL, mass_y = NA_REAL;
      double area = NA_REAL, area_ch = NA_REAL;
      double perimeter = NA_REAL;
      double radius_mean = NA_REAL, radius_min = NA_REAL, radius_max = NA_REAL, radius_sd = NA_REAL;
      double radius_ratio = NA_REAL;
      double diam_mean = NA_REAL, diam_min = NA_REAL, diam_max = NA_REAL;
      double caliper = NA_REAL, length_m = NA_REAL, width_m = NA_REAL;
      double solidity = NA_REAL, convexity = NA_REAL, elongation = NA_REAL;
      double circularity = NA_REAL, circularity_haralick = NA_REAL, circularity_norm = NA_REAL;
      double eccentricity = NA_REAL, maj_axis = NA_REAL, min_axis = NA_REAL, theta = NA_REAL;
      double coverage = NA_REAL, form_factor = NA_REAL, narrow_factor = NA_REAL;
      double asp_ratio = NA_REAL, rectangularity = NA_REAL, pd_ratio = NA_REAL, plw_ratio = NA_REAL;
      std::vector<std::pair<double, double>> hull;
      bool valid = false;
    };

    std::vector<ObjMetric> metrics(max_label);
    double total_pixels = static_cast<double>(N);

    // Pure C++ thread-safe parallel computations (No R API / SEXP allocations)
    #pragma omp parallel for schedule(dynamic) if(max_label > 10)
    for (int id = 1; id <= max_label; id++) {
      const auto& px_i = all_bx[id];
      const auto& py_i = all_by[id];
      int n_pts = px_i.size();
      if (n_pts < 3) continue;

      int i = id - 1;
      auto& m = metrics[i];
      m.valid = true;

      std::vector<double> px(n_pts), py(n_pts);
      for (int j = 0; j < n_pts; j++) {
        px[j] = static_cast<double>(px_i[j] + 2);
        py[j] = static_cast<double>(py_i[j] + 2);
      }

      // Shoelace
      double a = 0.0, cm_x = 0.0, cm_y = 0.0;
      for (int j = 0; j < n_pts; ++j) {
        int next_j = (j + 1) % n_pts;
        double x1 = px[j], y1 = py[j];
        double x2 = px[next_j], y2 = py[next_j];
        double cross = x1 * y2 - x2 * y1;
        a += cross;
        cm_x += (x1 + x2) * cross;
        cm_y += (y1 + y2) * cross;
      }
      double abs_area = std::abs(a / 2.0);
      m.area = abs_area;
      if (abs_area > 0.0 && a != 0.0) {
        m.mass_x = cm_x / (6.0 * (a / 2.0));
        m.mass_y = cm_y / (6.0 * (a / 2.0));
      } else {
        m.mass_x = 0.0; m.mass_y = 0.0;
      }

      // Convex Hull
      std::vector<std::pair<double, double>> pts(n_pts);
      for (int j = 0; j < n_pts; j++) pts[j] = {px[j], py[j]};
      std::sort(pts.begin(), pts.end());

      std::vector<std::pair<double, double>> hull(2 * n_pts);
      int k = 0;
      for (int j = 0; j < n_pts; ++j) {
        while (k >= 2) {
          double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                         (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
          if (cross <= 0) k--; else break;
        }
        hull[k++] = pts[j];
      }
      for (int j = n_pts - 2, t = k + 1; j >= 0; j--) {
        while (k >= t) {
          double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                         (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
          if (cross <= 0) k--; else break;
        }
        hull[k++] = pts[j];
      }

      int n_ch = k - 1;
      double a_ch = 0.0, p_ch = 0.0, max_d_sq = 0.0;

      m.hull.resize(n_ch + 1);
      for (int j = 0; j < n_ch; ++j) {
        m.hull[j] = hull[j];
        int next_j = (j + 1) % n_ch;
        a_ch += hull[j].first * hull[next_j].second - hull[next_j].first * hull[j].second;
        double dx = hull[next_j].first - hull[j].first;
        double dy = hull[next_j].second - hull[j].second;
        p_ch += std::sqrt(dx * dx + dy * dy);

        for (int q = j + 1; q < n_ch; q++) {
          double dxx = hull[j].first - hull[q].first;
          double dyy = hull[j].second - hull[q].second;
          double d_sq = dxx * dxx + dyy * dyy;
          if (d_sq > max_d_sq) max_d_sq = d_sq;
        }
      }
      m.hull[n_ch] = hull[0];

      m.area_ch = std::abs(a_ch / 2.0);
      double cal_val = std::sqrt(max_d_sq);
      m.caliper = cal_val;

      // Perimeter
      double p = 0.0;
      for (int j = 0; j < n_pts - 1; j++) {
        double dx = px[j+1] - px[j], dy = py[j+1] - py[j];
        p += std::sqrt(dx * dx + dy * dy);
      }
      m.perimeter = p;

      // Centroid distance
      double cent_x = 0.0, cent_y = 0.0;
      for (int j = 0; j < n_pts; j++) { cent_x += px[j]; cent_y += py[j]; }
      cent_x /= n_pts; cent_y /= n_pts;

      double c_mean = 0.0, c_min = 1e15, c_max = -1e15;
      std::vector<double> cdists(n_pts);
      for (int j = 0; j < n_pts; j++) {
        double dx = px[j] - cent_x, dy = py[j] - cent_y;
        double d = std::sqrt(dx * dx + dy * dy);
        cdists[j] = d;
        c_mean += d;
        if (d < c_min) c_min = d;
        if (d > c_max) c_max = d;
      }
      c_mean /= n_pts;

      double c_var = 0.0;
      for (int j = 0; j < n_pts; j++) c_var += (cdists[j] - c_mean) * (cdists[j] - c_mean);
      double c_sd = std::sqrt(c_var / (n_pts > 1 ? (n_pts - 1) : 1));

      m.radius_mean = c_mean; m.radius_min = c_min; m.radius_max = c_max; m.radius_sd = c_sd;
      m.radius_ratio = (c_min > 0) ? (c_max / c_min) : 0.0;
      m.diam_mean = c_mean * 2.0; m.diam_min = c_min * 2.0; m.diam_max = c_max * 2.0;

      // Moments & Axes
      double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
      for (int j = 0; j < n_pts; j++) {
        double x = px[j], y = py[j];
        sum_x += x; sum_y += y;
        sum_x2 += x * x; sum_y2 += y * y;
        sum_xy += x * y;
      }
      double cov_xy = sum_xy / n_pts - sum_x * sum_y / n_pts / n_pts;
      double var_x = sum_x2 / n_pts - sum_x * sum_x / n_pts / n_pts;
      double var_y = sum_y2 / n_pts - sum_y * sum_y / n_pts / n_pts;
      double t = 0.5 * std::atan2(2.0 * cov_xy, var_x - var_y);

      double cos_t = std::cos(t), sin_t = std::sin(t);
      double min_u = 1e15, max_u = -1e15, min_v = 1e15, max_v = -1e15;
      for (int j = 0; j < n_pts; j++) {
        double u = px[j] * cos_t + py[j] * sin_t;
        double v = -px[j] * sin_t + py[j] * cos_t;
        if (u < min_u) min_u = u;
        if (u > max_u) max_u = u;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
      }
      double dim1 = max_u - min_u, dim2 = max_v - min_v;
      double l = std::max(dim1, dim2), w = std::min(dim1, dim2);
      m.length_m = l; m.width_m = w;

      double a_maj = std::sqrt(std::max(0.0, 0.5 * (var_x + var_y + std::sqrt(std::pow(var_x - var_y, 2) + 4.0 * std::pow(cov_xy, 2)))));
      double b_min = std::sqrt(std::max(0.0, 0.5 * (var_x + var_y - std::sqrt(std::pow(var_x - var_y, 2) + 4.0 * std::pow(cov_xy, 2)))));
      m.maj_axis = std::fmax(a_maj, b_min); m.min_axis = std::fmin(a_maj, b_min);
      m.eccentricity = (m.maj_axis > 0) ? std::sqrt(std::max(0.0, 1.0 - std::pow(m.min_axis / m.maj_axis, 2))) : 0.0;
      m.theta = t;

      m.solidity = (m.area_ch > 0) ? (abs_area / m.area_ch) : 0.0;
      m.convexity = (p > 0) ? (p_ch / p) : 0.0;
      m.elongation = (l > 0) ? (1.0 - (w / l)) : 0.0;
      m.circularity = (abs_area > 0) ? ((p * p) / abs_area) : 0.0;
      m.circularity_haralick = (c_sd > 0) ? (c_mean / c_sd) : 0.0;
      m.circularity_norm = (p > 0) ? ((abs_area * 4.0 * M_PI) / (p * p)) : 0.0;

      m.coverage = abs_area / total_pixels;
      m.form_factor = (p > 0) ? (4.0 * M_PI * abs_area / (p * p)) : 0.0;
      m.narrow_factor = (l > 0) ? (cal_val / l) : 0.0;
      m.asp_ratio = (w > 0) ? (l / w) : 0.0;
      m.rectangularity = (abs_area > 0) ? (l * w / abs_area) : 0.0;
      m.pd_ratio = (cal_val > 0) ? (p / cal_val) : 0.0;
      m.plw_ratio = ((l + w) > 0) ? (p / (l + w)) : 0.0;
    }

    // Serial Rcpp Object Construction (Thread-safe)
    ocont = List(max_label);
    ch_list = List(max_label);
    CharacterVector col_names = CharacterVector::create("x", "y");

    std::vector<double> mass_x(max_label), mass_y(max_label);
    std::vector<double> area(max_label), area_ch(max_label);
    std::vector<double> perimeter(max_label);
    std::vector<double> radius_mean(max_label), radius_min(max_label), radius_max(max_label), radius_sd(max_label);
    std::vector<double> radius_ratio(max_label);
    std::vector<double> diam_mean(max_label), diam_min(max_label), diam_max(max_label);
    std::vector<double> caliper(max_label), length_m(max_label), width_m(max_label);
    std::vector<double> solidity(max_label), convexity(max_label), elongation(max_label);
    std::vector<double> circularity(max_label), circularity_haralick(max_label), circularity_norm(max_label);
    std::vector<double> eccentricity(max_label), maj_axis(max_label), min_axis(max_label), theta(max_label);
    std::vector<double> coverage(max_label), form_factor(max_label), narrow_factor(max_label);
    std::vector<double> asp_ratio(max_label), rectangularity(max_label), pd_ratio(max_label), plw_ratio(max_label);

    for (int id = 1; id <= max_label; id++) {
      int i = id - 1;
      const auto& m = metrics[i];
      if (!m.valid) {
        mass_x[i] = NA_REAL; mass_y[i] = NA_REAL;
        area[i] = NA_REAL; area_ch[i] = NA_REAL; perimeter[i] = NA_REAL;
        radius_mean[i] = NA_REAL; radius_min[i] = NA_REAL; radius_max[i] = NA_REAL; radius_sd[i] = NA_REAL;
        radius_ratio[i] = NA_REAL; diam_mean[i] = NA_REAL; diam_min[i] = NA_REAL; diam_max[i] = NA_REAL;
        caliper[i] = NA_REAL; length_m[i] = NA_REAL; width_m[i] = NA_REAL;
        solidity[i] = NA_REAL; convexity[i] = NA_REAL; elongation[i] = NA_REAL;
        circularity[i] = NA_REAL; circularity_haralick[i] = NA_REAL; circularity_norm[i] = NA_REAL;
        eccentricity[i] = NA_REAL; maj_axis[i] = NA_REAL; min_axis[i] = NA_REAL; theta[i] = NA_REAL;
        coverage[i] = NA_REAL; form_factor[i] = NA_REAL; narrow_factor[i] = NA_REAL;
        asp_ratio[i] = NA_REAL; rectangularity[i] = NA_REAL; pd_ratio[i] = NA_REAL; plw_ratio[i] = NA_REAL;
        continue;
      }

      const auto& px_i = all_bx[id];
      const auto& py_i = all_by[id];
      int n_pts = px_i.size();

      IntegerMatrix cmat(n_pts, 2);
      int* cmat_ptr = INTEGER(cmat);
      for (int j = 0; j < n_pts; j++) {
        cmat_ptr[j] = px_i[j] + 2;
        cmat_ptr[j + n_pts] = py_i[j] + 2;
      }
      cmat.attr("dimnames") = List::create(R_NilValue, col_names);
      ocont[i] = cmat;

      int n_ch = m.hull.size();
      if (n_ch > 0) {
        NumericMatrix ch_mat(n_ch, 2);
        double* ch_ptr = REAL(ch_mat);
        for (int j = 0; j < n_ch; j++) {
          ch_ptr[j] = m.hull[j].first;
          ch_ptr[j + n_ch] = m.hull[j].second;
        }
        ch_mat.attr("dimnames") = List::create(R_NilValue, col_names);
        ch_list[i] = ch_mat;
      }

      mass_x[i] = m.mass_x; mass_y[i] = m.mass_y;
      area[i] = m.area; area_ch[i] = m.area_ch; perimeter[i] = m.perimeter;
      radius_mean[i] = m.radius_mean; radius_min[i] = m.radius_min; radius_max[i] = m.radius_max; radius_sd[i] = m.radius_sd;
      radius_ratio[i] = m.radius_ratio; diam_mean[i] = m.diam_mean; diam_min[i] = m.diam_min; diam_max[i] = m.diam_max;
      caliper[i] = m.caliper; length_m[i] = m.length_m; width_m[i] = m.width_m;
      solidity[i] = m.solidity; convexity[i] = m.convexity; elongation[i] = m.elongation;
      circularity[i] = m.circularity; circularity_haralick[i] = m.circularity_haralick; circularity_norm[i] = m.circularity_norm;
      eccentricity[i] = m.eccentricity; maj_axis[i] = m.maj_axis; min_axis[i] = m.min_axis; theta[i] = m.theta;
      coverage[i] = m.coverage; form_factor[i] = m.form_factor; narrow_factor[i] = m.narrow_factor;
      asp_ratio[i] = m.asp_ratio; rectangularity[i] = m.rectangularity; pd_ratio[i] = m.pd_ratio; plw_ratio[i] = m.plw_ratio;
    }

    std::vector<int> id_seq(max_label);
    for (int i = 0; i < max_label; i++) id_seq[i] = i + 1;

    shape = DataFrame::create(
      Named("id") = wrap(id_seq),
      Named("x") = wrap(mass_x),
      Named("y") = wrap(mass_y),
      Named("area") = wrap(area),
      Named("area_ch") = wrap(area_ch),
      Named("perimeter") = wrap(perimeter),
      Named("radius_mean") = wrap(radius_mean),
      Named("radius_min") = wrap(radius_min),
      Named("radius_max") = wrap(radius_max),
      Named("radius_sd") = wrap(radius_sd),
      Named("diam_mean") = wrap(diam_mean),
      Named("diam_min") = wrap(diam_min),
      Named("diam_max") = wrap(diam_max),
      Named("major_axis") = wrap(maj_axis),
      Named("minor_axis") = wrap(min_axis),
      Named("caliper") = wrap(caliper),
      Named("length") = wrap(length_m),
      Named("width") = wrap(width_m),
      Named("radius_ratio") = wrap(radius_ratio),
      Named("theta") = wrap(theta),
      Named("eccentricity") = wrap(eccentricity),
      Named("form_factor") = wrap(form_factor),
      Named("narrow_factor") = wrap(narrow_factor),
      Named("asp_ratio") = wrap(asp_ratio),
      Named("rectangularity") = wrap(rectangularity),
      Named("pd_ratio") = wrap(pd_ratio),
      Named("plw_ratio") = wrap(plw_ratio),
      Named("solidity") = wrap(solidity),
      Named("convexity") = wrap(convexity),
      Named("elongation") = wrap(elongation),
      Named("circularity") = wrap(circularity),
      Named("circularity_haralick") = wrap(circularity_haralick),
      Named("circularity_norm") = wrap(circularity_norm),
      Named("coverage") = wrap(coverage)
    );
  }

  // Step 5: Haralick Features (Optional)
  NumericMatrix haralick_mat;
  if (haralick && ocont.size() > 0 && !Rf_isNull(img_sexp)) {
    SEXP ref_chan = R_NilValue;
    SEXP dims = Rf_getAttrib(img_sexp, R_DimSymbol);
    if (!Rf_isNull(dims) && Rf_length(dims) >= 3) {
      int w = INTEGER(dims)[0];
      int h = INTEGER(dims)[1];
      int nch = INTEGER(dims)[2];
      int npix = w * h;
      int target_ch = (har_band >= 1 && har_band <= nch) ? (har_band - 1) : 0;

      if (TYPEOF(img_sexp) == RAWSXP) {
        SEXP chan = PROTECT(Rf_allocMatrix(RAWSXP, w, h));
        std::memcpy(RAW(chan), RAW(img_sexp) + target_ch * npix, npix);
        ref_chan = chan;
        UNPROTECT(1);
      } else if (TYPEOF(img_sexp) == REALSXP) {
        SEXP chan = PROTECT(Rf_allocMatrix(REALSXP, w, h));
        std::memcpy(REAL(chan), REAL(img_sexp) + target_ch * npix, npix * sizeof(double));
        ref_chan = chan;
        UNPROTECT(1);
      }
    } else {
      ref_chan = img_sexp;
    }
    if (!Rf_isNull(ref_chan)) {
      haralick_mat = haralick_features_cpp(labels, ref_chan, har_nbins);
    }
  }

  return List::create(
    _["labels"] = labels,
    _["shape"] = shape,
    _["contours"] = ocont,
    _["chull"] = ch_list,
    _["haralick"] = haralick_mat
  );
}

