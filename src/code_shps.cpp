#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// =============================================================================
// NEW OPTIMIZED ENGINE (OpenMP Multi-threading + Binary Search + Low-level R C API)
// =============================================================================

inline double dist_sq(double x1, double y1, double x2, double y2) {
  double dx = x2 - x1;
  double dy = y2 - y1;
  return dx * dx + dy * dy;
}

inline double dist_eucl(double x1, double y1, double x2, double y2) {
  return std::sqrt(dist_sq(x1, y1, x2, y2));
}

inline void get_point_at_dist_fast(const double* coords, const double* cum_dist, int n, double target_dist, double* out_xy) {
  if (target_dist <= 0.0) {
    out_xy[0] = coords[0];
    out_xy[1] = coords[n];
    return;
  }
  if (target_dist >= cum_dist[n - 1]) {
    out_xy[0] = coords[n - 1];
    out_xy[1] = coords[n - 1 + n];
    return;
  }
  const double* it = std::upper_bound(cum_dist, cum_dist + n, target_dist);
  int i = std::distance(cum_dist, it) - 1;
  if (i < 0) i = 0;
  if (i >= n - 1) i = n - 2;

  double seg_len = cum_dist[i + 1] - cum_dist[i];
  if (seg_len <= 0.0) {
    out_xy[0] = coords[i];
    out_xy[1] = coords[i + n];
    return;
  }
  double t = (target_dist - cum_dist[i]) / seg_len;
  out_xy[0] = coords[i] * (1.0 - t) + coords[i + 1] * t;
  out_xy[1] = coords[i + n] * (1.0 - t) + coords[i + 1 + n] * t;
}

// [[Rcpp::export]]
NumericVector get_point_at_dist(NumericMatrix coords, NumericVector cum_dist, double target_dist) {
  NumericVector out(2);
  int n = coords.nrow();
  const double* c_ptr = REAL(coords);
  const double* d_ptr = REAL(cum_dist);
  get_point_at_dist_fast(c_ptr, d_ptr, n, target_dist, REAL(out));
  return out;
}

// [[Rcpp::export]]
int get_closest_idx_forward(NumericMatrix rail, double tx, double ty, int start_idx) {
  int n = rail.nrow();
  const double* r_ptr = REAL(rail);
  int best_idx = start_idx;
  double min_d2 = dist_sq(r_ptr[start_idx], r_ptr[start_idx + n], tx, ty);
  int search_limit = std::min(n, start_idx + 1000);

  for (int i = start_idx + 1; i < search_limit; i++) {
    double d2 = dist_sq(r_ptr[i], r_ptr[i + n], tx, ty);
    if (d2 < min_d2) {
      min_d2 = d2;
      best_idx = i;
    }
  }
  return best_idx;
}

// [[Rcpp::export]]
SEXP make_grid_structure(NumericMatrix rail1,
                         NumericMatrix rail2,
                         int nrow,
                         int ncol,
                         double buffer_col,
                         double buffer_row,
                         Nullable<double> plot_width_opt,
                         Nullable<double> plot_height_opt) {
  int n1 = rail1.nrow();
  const double* r1_ptr = REAL(rail1);
  std::vector<double> cum_dist1(n1, 0.0);
  for (int i = 1; i < n1; i++) {
    cum_dist1[i] = cum_dist1[i - 1] + dist_eucl(r1_ptr[i - 1], r1_ptr[i - 1 + n1], r1_ptr[i], r1_ptr[i + n1]);
  }
  double total_len1 = cum_dist1[n1 - 1];

  int n2 = rail2.nrow();
  const double* r2_ptr = REAL(rail2);
  std::vector<double> cum_dist2(n2, 0.0);
  for (int i = 1; i < n2; i++) {
    cum_dist2[i] = cum_dist2[i - 1] + dist_eucl(r2_ptr[i - 1], r2_ptr[i - 1 + n2], r2_ptr[i], r2_ptr[i + n2]);
  }
  double total_len2 = cum_dist2[n2 - 1];

  double cell_along_avg = (total_len1 / ncol + total_len2 / ncol) / 2.0;
  double margin_along = 0.0;
  if (plot_width_opt.isNotNull()) {
    double pw = as<double>(plot_width_opt);
    margin_along = std::max(0.0, (cell_along_avg - pw) / 2.0);
  } else if (buffer_col > 0) {
    margin_along = buffer_col / 2.0;
  }

  int total_plots = nrow * ncol;
  std::vector<std::array<double, 10>> ring_data(total_plots);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ncol; i++) {
    double t_base_start = (double)i / ncol;
    double t_base_end   = (double)(i + 1) / ncol;
    double dist_start_l1 = t_base_start * total_len1 + margin_along;
    double dist_end_l1   = t_base_end * total_len1   - margin_along;
    double dist_start_l2 = t_base_start * total_len2 + margin_along;
    double dist_end_l2   = t_base_end * total_len2   - margin_along;

    double p_r1_start[2], p_r1_end[2], p_r2_start[2], p_r2_end[2];
    get_point_at_dist_fast(r1_ptr, cum_dist1.data(), n1, dist_start_l1, p_r1_start);
    get_point_at_dist_fast(r1_ptr, cum_dist1.data(), n1, dist_end_l1,   p_r1_end);
    get_point_at_dist_fast(r2_ptr, cum_dist2.data(), n2, dist_start_l2, p_r2_start);
    get_point_at_dist_fast(r2_ptr, cum_dist2.data(), n2, dist_end_l2,   p_r2_end);

    double w_top = dist_eucl(p_r1_start[0], p_r1_start[1], p_r2_start[0], p_r2_start[1]);
    double w_bot = dist_eucl(p_r1_end[0],   p_r1_end[1],   p_r2_end[0],   p_r2_end[1]);
    double cell_cross_avg = (w_top + w_bot) / 2.0 / nrow;

    double margin_cross = 0.0;
    if (plot_height_opt.isNotNull()) {
      double ph = as<double>(plot_height_opt);
      margin_cross = std::max(0.0, (cell_cross_avg - ph) / 2.0);
    } else if (buffer_row > 0) {
      margin_cross = buffer_row / 2.0;
    }

    int col_offset = i * nrow;

    for (int j = 0; j < nrow; j++) {
      double u_base_start = (double)j / nrow;
      double u_base_end   = (double)(j + 1) / nrow;
      double u_margin_top = margin_cross / w_top;
      double u_margin_bot = margin_cross / w_bot;
      double u_start_top = u_base_start + u_margin_top;
      double u_end_top   = u_base_end   - u_margin_top;
      double u_start_bot = u_base_start + u_margin_bot;
      double u_end_bot   = u_base_end   - u_margin_bot;

      double x1 = p_r1_start[0] * (1 - u_start_top) + p_r2_start[0] * u_start_top;
      double y1 = p_r1_start[1] * (1 - u_start_top) + p_r2_start[1] * u_start_top;
      double x2 = p_r1_start[0] * (1 - u_end_top)   + p_r2_start[0] * u_end_top;
      double y2 = p_r1_start[1] * (1 - u_end_top)   + p_r2_start[1] * u_end_top;
      double x3 = p_r1_end[0] * (1 - u_end_bot)   + p_r2_end[0] * u_end_bot;
      double y3 = p_r1_end[1] * (1 - u_end_bot)   + p_r2_end[1] * u_end_bot;
      double x4 = p_r1_end[0] * (1 - u_start_bot) + p_r2_end[0] * u_start_bot;
      double y4 = p_r1_end[1] * (1 - u_start_bot) + p_r2_end[1] * u_start_bot;

      auto& r = ring_data[col_offset + j];
      r[0] = x1; r[5] = y1;
      r[1] = x2; r[6] = y2;
      r[2] = x3; r[7] = y3;
      r[3] = x4; r[8] = y4;
      r[4] = x1; r[9] = y1;
    }
  }

  // Fast R C API allocation
  SEXP sfg_class_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(sfg_class_sexp, 0, Rf_mkChar("XY"));
  SET_STRING_ELT(sfg_class_sexp, 1, Rf_mkChar("POLYGON"));
  SET_STRING_ELT(sfg_class_sexp, 2, Rf_mkChar("sfg"));

  SEXP out_list_sexp = PROTECT(Rf_allocVector(VECSXP, total_plots));

  for (int p = 0; p < total_plots; p++) {
    const auto& r = ring_data[p];
    SEXP ring = PROTECT(Rf_allocMatrix(REALSXP, 5, 2));
    double* r_ptr = REAL(ring);
    r_ptr[0] = r[0]; r_ptr[5] = r[5];
    r_ptr[1] = r[1]; r_ptr[6] = r[6];
    r_ptr[2] = r[2]; r_ptr[7] = r[7];
    r_ptr[3] = r[3]; r_ptr[8] = r[8];
    r_ptr[4] = r[4]; r_ptr[9] = r[9];

    SEXP sfg = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(sfg, 0, ring);
    Rf_setAttrib(sfg, R_ClassSymbol, sfg_class_sexp);

    SET_VECTOR_ELT(out_list_sexp, p, sfg);
    UNPROTECT(2);
  }

  UNPROTECT(2);
  return out_list_sexp;
}

// [[Rcpp::export]]
SEXP make_grid_curved(NumericMatrix rail1,
                      NumericMatrix rail2,
                      int nrow,
                      int ncol,
                      bool curved = true,
                      int density = 20) {
  int n_dense = rail1.nrow();
  const double* r1_ptr = REAL(rail1);
  const double* r2_ptr = REAL(rail2);

  std::vector<double> cl_x(n_dense);
  std::vector<double> cl_y(n_dense);
  std::vector<double> cl_cum_dist(n_dense, 0.0);

  cl_x[0] = (r1_ptr[0] + r2_ptr[0]) / 2.0;
  cl_y[0] = (r1_ptr[n_dense] + r2_ptr[n_dense]) / 2.0;

  for (int i = 1; i < n_dense; i++) {
    cl_x[i] = (r1_ptr[i] + r2_ptr[i]) / 2.0;
    cl_y[i] = (r1_ptr[i + n_dense] + r2_ptr[i + n_dense]) / 2.0;
    cl_cum_dist[i] = cl_cum_dist[i - 1] + dist_eucl(cl_x[i - 1], cl_y[i - 1], cl_x[i], cl_y[i]);
  }
  double total_cl_len = cl_cum_dist[n_dense - 1];

  int total_plots = nrow * ncol;
  int n_steps = curved ? density : 1;
  int ring_rows = 2 * (n_steps + 1) + 1;

  std::vector<std::vector<double>> ring_data(total_plots, std::vector<double>(ring_rows * 2));
  std::vector<int> r1_s(ncol), r2_s(ncol), r1_e(ncol), r2_e(ncol);
  int last_idx_r1 = 0, last_idx_r2 = 0;

  for (int i = 0; i < ncol; i++) {
    if (i == 0) {
      r1_s[i] = 0; r2_s[i] = 0;
    } else {
      double dist_s = ((double)i / ncol) * total_cl_len;
      const double* c_ptr = cl_cum_dist.data();
      const double* it = std::upper_bound(c_ptr, c_ptr + n_dense, dist_s);

      int idx_cl = std::distance(c_ptr, it) - 1;
      if (idx_cl < 0) idx_cl = 0;

      r1_s[i] = get_closest_idx_forward(rail1, cl_x[idx_cl], cl_y[idx_cl], last_idx_r1);
      r2_s[i] = get_closest_idx_forward(rail2, cl_x[idx_cl], cl_y[idx_cl], last_idx_r2);
    }
    if (i == ncol - 1) {
      r1_e[i] = n_dense - 1; r2_e[i] = n_dense - 1;
    } else {
      double dist_e = ((double)(i + 1) / ncol) * total_cl_len;
      const double* c_ptr = cl_cum_dist.data();
      const double* it = std::upper_bound(c_ptr, c_ptr + n_dense, dist_e);
      int idx_cl = std::distance(c_ptr, it) - 1;

      if (idx_cl < 0) idx_cl = 0;

      r1_e[i] = get_closest_idx_forward(rail1, cl_x[idx_cl], cl_y[idx_cl], r1_s[i]);
      r2_e[i] = get_closest_idx_forward(rail2, cl_x[idx_cl], cl_y[idx_cl], r2_s[i]);
    }
    last_idx_r1 = r1_s[i];
    last_idx_r2 = r2_s[i];
  }

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ncol; i++) {
    int idx_r1_s = r1_s[i], idx_r2_s = r2_s[i];
    int idx_r1_e = r1_e[i], idx_r2_e = r2_e[i];
    int col_offset = i * nrow;

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      auto& r_vec = ring_data[col_offset + j];
      int pt_idx = 0;

      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_s + (int)(t * (idx_r1_e - idx_r1_s));
        int i2 = idx_r2_s + (int)(t * (idx_r2_e - idx_r2_s));

        double x_r1 = r1_ptr[i1], y_r1 = r1_ptr[i1 + n_dense];
        double x_r2 = r2_ptr[i2], y_r2 = r2_ptr[i2 + n_dense];

        r_vec[pt_idx]             = x_r1 * (1.0 - u_top) + x_r2 * u_top;
        r_vec[pt_idx + ring_rows] = y_r1 * (1.0 - u_top) + y_r2 * u_top;
        pt_idx++;
      }
      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_s + (int)(t * (idx_r1_e - idx_r1_s));
        int i2 = idx_r2_s + (int)(t * (idx_r2_e - idx_r2_s));

        double x_r1 = r1_ptr[i1], y_r1 = r1_ptr[i1 + n_dense];
        double x_r2 = r2_ptr[i2], y_r2 = r2_ptr[i2 + n_dense];

        r_vec[pt_idx]             = x_r1 * (1.0 - u_bot) + x_r2 * u_bot;
        r_vec[pt_idx + ring_rows] = y_r1 * (1.0 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }
      r_vec[pt_idx]             = r_vec[0];
      r_vec[pt_idx + ring_rows] = r_vec[ring_rows];
    }
  }

  SEXP sfg_class_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(sfg_class_sexp, 0, Rf_mkChar("XY"));
  SET_STRING_ELT(sfg_class_sexp, 1, Rf_mkChar("POLYGON"));
  SET_STRING_ELT(sfg_class_sexp, 2, Rf_mkChar("sfg"));

  SEXP out_list_sexp = PROTECT(Rf_allocVector(VECSXP, total_plots));

  for (int p = 0; p < total_plots; p++) {
    const auto& r_vec = ring_data[p];
    SEXP ring = PROTECT(Rf_allocMatrix(REALSXP, ring_rows, 2));
    double* r_ptr = REAL(ring);
    std::copy(r_vec.begin(), r_vec.end(), r_ptr);

    SEXP sfg = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(sfg, 0, ring);
    Rf_setAttrib(sfg, R_ClassSymbol, sfg_class_sexp);

    SET_VECTOR_ELT(out_list_sexp, p, sfg);
    UNPROTECT(2);
  }

  UNPROTECT(2);
  return out_list_sexp;
}

// [[Rcpp::export]]
SEXP make_grid_landmarks(NumericMatrix rail1,
                         NumericMatrix rail2,
                         IntegerVector anchors1,
                         IntegerVector anchors2,
                         int nrow,
                         bool curved = true,
                         int density = 30) {
  int n_cols = anchors1.size() - 1;
  if (anchors2.size() != anchors1.size()) {
    stop("Rail 1 and Rail 2 must have the same number of control points for manual mode.");
  }

  int n_dense = rail1.nrow();
  const double* r1_ptr = REAL(rail1);
  const double* r2_ptr = REAL(rail2);

  int total_plots = nrow * n_cols;
  int n_steps = curved ? density : 1;
  int ring_rows = 2 * (n_steps + 1) + 1;

  std::vector<std::vector<double>> ring_data(total_plots, std::vector<double>(ring_rows * 2));

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_cols; i++) {
    int idx_r1_start = anchors1[i];
    int idx_r1_end   = anchors1[i + 1];
    int idx_r2_start = anchors2[i];
    int idx_r2_end   = anchors2[i + 1];

    double diff_r1 = (double)(idx_r1_end - idx_r1_start);
    double diff_r2 = (double)(idx_r2_end - idx_r2_start);
    int col_offset = i * nrow;

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      auto& r_vec = ring_data[col_offset + j];
      int pt_idx = 0;

      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);
        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }

        double x_r1 = r1_ptr[i1], y_r1 = r1_ptr[i1 + n_dense];
        double x_r2 = r2_ptr[i2], y_r2 = r2_ptr[i2 + n_dense];

        r_vec[pt_idx]             = x_r1 * (1.0 - u_top) + x_r2 * u_top;
        r_vec[pt_idx + ring_rows] = y_r1 * (1.0 - u_top) + y_r2 * u_top;
        pt_idx++;
      }

      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);
        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }

        double x_r1 = r1_ptr[i1], y_r1 = r1_ptr[i1 + n_dense];
        double x_r2 = r2_ptr[i2], y_r2 = r2_ptr[i2 + n_dense];

        r_vec[pt_idx]             = x_r1 * (1.0 - u_bot) + x_r2 * u_bot;
        r_vec[pt_idx + ring_rows] = y_r1 * (1.0 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }

      r_vec[pt_idx]             = r_vec[0];
      r_vec[pt_idx + ring_rows] = r_vec[ring_rows];
    }
  }

  SEXP sfg_class_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(sfg_class_sexp, 0, Rf_mkChar("XY"));
  SET_STRING_ELT(sfg_class_sexp, 1, Rf_mkChar("POLYGON"));
  SET_STRING_ELT(sfg_class_sexp, 2, Rf_mkChar("sfg"));

  SEXP out_list_sexp = PROTECT(Rf_allocVector(VECSXP, total_plots));

  for (int p = 0; p < total_plots; p++) {
    const auto& r_vec = ring_data[p];
    SEXP ring = PROTECT(Rf_allocMatrix(REALSXP, ring_rows, 2));
    double* r_ptr = REAL(ring);
    std::copy(r_vec.begin(), r_vec.end(), r_ptr);

    SEXP sfg = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(sfg, 0, ring);
    Rf_setAttrib(sfg, R_ClassSymbol, sfg_class_sexp);

    SET_VECTOR_ELT(out_list_sexp, p, sfg);
    UNPROTECT(2);
  }

  UNPROTECT(2);
  return out_list_sexp;
}

// =============================================================================
// OLD ENGINE (Single-threaded baseline for comparison)
// =============================================================================

double dist_eucl_old(double x1, double y1, double x2, double y2) {
  return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

double dist_sq_old(double x1, double y1, double x2, double y2) {
  return std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2);
}

NumericVector get_point_at_dist_old(NumericMatrix coords, NumericVector cum_dist, double target_dist) {
  int n = coords.nrow();
  if (target_dist <= 0) return coords(0, _);
  if (target_dist >= cum_dist[n-1]) return coords(n-1, _);
  int i = 0;
  while(i < n - 1 && cum_dist[i+1] < target_dist) i++;
  double segment_len = cum_dist[i+1] - cum_dist[i];
  if (segment_len == 0) return coords(i, _);
  double t = (target_dist - cum_dist[i]) / segment_len;
  NumericVector out(2);
  out[0] = coords(i, 0) * (1 - t) + coords(i+1, 0) * t;
  out[1] = coords(i, 1) * (1 - t) + coords(i+1, 1) * t;
  return out;
}

int get_closest_idx_forward_old(NumericMatrix rail, double tx, double ty, int start_idx) {
  int n = rail.nrow();
  int best_idx = start_idx;
  double min_d2 = dist_sq_old(rail(start_idx, 0), rail(start_idx, 1), tx, ty);
  int search_limit = std::min(n, start_idx + 1000);

  for(int i = start_idx + 1; i < search_limit; i++) {
    double d2 = dist_sq_old(rail(i, 0), rail(i, 1), tx, ty);
    if (d2 < min_d2) {
      min_d2 = d2;
      best_idx = i;
    }
  }
  return best_idx;
}

// [[Rcpp::export]]
List make_grid_structure_old(NumericMatrix rail1,
                             NumericMatrix rail2,
                             int nrow,
                             int ncol,
                             double buffer_col,
                             double buffer_row,
                             Nullable<double> plot_width_opt,
                             Nullable<double> plot_height_opt) {
  int n1 = rail1.nrow();
  NumericVector cum_dist1(n1); cum_dist1[0] = 0;
  for(int i=1; i<n1; i++) cum_dist1[i] = cum_dist1[i-1] + dist_eucl_old(rail1(i-1,0), rail1(i-1,1), rail1(i,0), rail1(i,1));
  double total_len1 = cum_dist1[n1-1];

  int n2 = rail2.nrow();
  NumericVector cum_dist2(n2); cum_dist2[0] = 0;
  for(int i=1; i<n2; i++) cum_dist2[i] = cum_dist2[i-1] + dist_eucl_old(rail2(i-1,0), rail2(i-1,1), rail2(i,0), rail2(i,1));
  double total_len2 = cum_dist2[n2-1];

  double cell_along_avg = (total_len1 / ncol + total_len2 / ncol) / 2.0;
  double margin_along = 0.0;
  if (plot_width_opt.isNotNull()) {
    double pw = as<double>(plot_width_opt);
    margin_along = std::max(0.0, (cell_along_avg - pw) / 2.0);
  } else if (buffer_col > 0) {
    margin_along = buffer_col / 2.0;
  }

  List out_list(nrow * ncol);
  int idx = 0;

  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");

  for (int i = 0; i < ncol; i++) {
    double t_base_start = (double)i / ncol;
    double t_base_end   = (double)(i + 1) / ncol;
    double dist_start_l1 = t_base_start * total_len1 + margin_along;
    double dist_end_l1   = t_base_end * total_len1   - margin_along;
    double dist_start_l2 = t_base_start * total_len2 + margin_along;
    double dist_end_l2   = t_base_end * total_len2   - margin_along;

    NumericVector p_r1_start = get_point_at_dist_old(rail1, cum_dist1, dist_start_l1);
    NumericVector p_r1_end   = get_point_at_dist_old(rail1, cum_dist1, dist_end_l1);
    NumericVector p_r2_start = get_point_at_dist_old(rail2, cum_dist2, dist_start_l2);
    NumericVector p_r2_end   = get_point_at_dist_old(rail2, cum_dist2, dist_end_l2);

    double w_top = dist_eucl_old(p_r1_start[0], p_r1_start[1], p_r2_start[0], p_r2_start[1]);
    double w_bot = dist_eucl_old(p_r1_end[0],   p_r1_end[1],   p_r2_end[0],   p_r2_end[1]);
    double cell_cross_avg = (w_top + w_bot) / 2.0 / nrow;

    double margin_cross = 0.0;
    if (plot_height_opt.isNotNull()) {
      double ph = as<double>(plot_height_opt);
      margin_cross = std::max(0.0, (cell_cross_avg - ph) / 2.0);
    } else if (buffer_row > 0) {
      margin_cross = buffer_row / 2.0;
    }
    for (int j = 0; j < nrow; j++) {
      double u_base_start = (double)j / nrow;
      double u_base_end   = (double)(j + 1) / nrow;
      double u_margin_top = margin_cross / w_top;
      double u_margin_bot = margin_cross / w_bot;
      double u_start_top = u_base_start + u_margin_top;
      double u_end_top   = u_base_end   - u_margin_top;
      double u_start_bot = u_base_start + u_margin_bot;
      double u_end_bot   = u_base_end   - u_margin_bot;
      double x1 = p_r1_start[0] * (1-u_start_top) + p_r2_start[0] * u_start_top;
      double y1 = p_r1_start[1] * (1-u_start_top) + p_r2_start[1] * u_start_top;
      double x2 = p_r1_start[0] * (1-u_end_top)   + p_r2_start[0] * u_end_top;
      double y2 = p_r1_start[1] * (1-u_end_top)   + p_r2_start[1] * u_end_top;
      double x3 = p_r1_end[0] * (1-u_end_bot)   + p_r2_end[0] * u_end_bot;
      double y3 = p_r1_end[1] * (1-u_end_bot)   + p_r2_end[1] * u_end_bot;
      double x4 = p_r1_end[0] * (1-u_start_bot) + p_r2_end[0] * u_start_bot;
      double y4 = p_r1_end[1] * (1-u_start_bot) + p_r2_end[1] * u_start_bot;
      NumericMatrix ring(5, 2);
      ring(0,0) = x1; ring(0,1) = y1;
      ring(1,0) = x2; ring(1,1) = y2;
      ring(2,0) = x3; ring(2,1) = y3;
      ring(3,0) = x4; ring(3,1) = y4;
      ring(4,0) = x1; ring(4,1) = y1;
      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;

      out_list[idx] = polygon_sfg;
      idx++;
    }
  }
  return out_list;
}

// [[Rcpp::export]]
List make_grid_curved_old(NumericMatrix rail1,
                          NumericMatrix rail2,
                          int nrow,
                          int ncol,
                          bool curved = true,
                          int density = 20) {

  int n_dense = rail1.nrow();
  NumericMatrix centerline(n_dense, 2);
  NumericVector cl_cum_dist(n_dense);
  cl_cum_dist[0] = 0;

  centerline(0, _) = (rail1(0, _) + rail2(0, _)) / 2.0;

  for(int i=1; i<n_dense; i++) {
    centerline(i, 0) = (rail1(i, 0) + rail2(i, 0)) / 2.0;
    centerline(i, 1) = (rail1(i, 1) + rail2(i, 1)) / 2.0;
    cl_cum_dist[i] = cl_cum_dist[i-1] + dist_eucl_old(centerline(i-1,0), centerline(i-1,1), centerline(i,0), centerline(i,1));
  }
  double total_cl_len = cl_cum_dist[n_dense-1];

  List out_list(nrow * ncol);
  int idx = 0;
  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");

  int last_idx_r1 = 0;
  int last_idx_r2 = 0;
  int n_steps = curved ? density : 1;

  for (int i = 0; i < ncol; i++) {
    int idx_r1_s, idx_r2_s, idx_r1_e, idx_r2_e;
    if (i == 0) {
      idx_r1_s = 0; idx_r2_s = 0;
    } else {
      double dist_s = ((double)i / ncol) * total_cl_len;
      NumericVector p_cl_s = get_point_at_dist_old(centerline, cl_cum_dist, dist_s);
      idx_r1_s = get_closest_idx_forward_old(rail1, p_cl_s[0], p_cl_s[1], last_idx_r1);
      idx_r2_s = get_closest_idx_forward_old(rail2, p_cl_s[0], p_cl_s[1], last_idx_r2);
    }
    if (i == ncol - 1) {
      idx_r1_e = n_dense - 1; idx_r2_e = n_dense - 1;
    } else {
      double dist_e = ((double)(i + 1) / ncol) * total_cl_len;
      NumericVector p_cl_e = get_point_at_dist_old(centerline, cl_cum_dist, dist_e);
      idx_r1_e = get_closest_idx_forward_old(rail1, p_cl_e[0], p_cl_e[1], idx_r1_s);
      idx_r2_e = get_closest_idx_forward_old(rail2, p_cl_e[0], p_cl_e[1], idx_r2_s);
    }
    last_idx_r1 = idx_r1_s;
    last_idx_r2 = idx_r2_s;

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      NumericMatrix ring(2 * (n_steps + 1) + 1, 2);
      int pt_idx = 0;

      // Top Edge
      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;
        double curr_idx_r1 = idx_r1_s + t * (idx_r1_e - idx_r1_s);
        double curr_idx_r2 = idx_r2_s + t * (idx_r2_e - idx_r2_s);

        int i1 = (int)curr_idx_r1;
        int i2 = (int)curr_idx_r2;

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_top) + x_r2 * u_top;
        ring(pt_idx, 1) = y_r1 * (1 - u_top) + y_r2 * u_top;
        pt_idx++;
      }
      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        double curr_idx_r1 = idx_r1_s + t * (idx_r1_e - idx_r1_s);
        double curr_idx_r2 = idx_r2_s + t * (idx_r2_e - idx_r2_s);

        int i1 = (int)curr_idx_r1;
        int i2 = (int)curr_idx_r2;

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_bot) + x_r2 * u_bot;
        ring(pt_idx, 1) = y_r1 * (1 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }
      ring(pt_idx, 0) = ring(0, 0);
      ring(pt_idx, 1) = ring(0, 1);

      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;
      out_list[idx++] = polygon_sfg;
    }
  }
  return out_list;
}

// [[Rcpp::export]]
List make_grid_landmarks_old(NumericMatrix rail1,
                             NumericMatrix rail2,
                             IntegerVector anchors1,
                             IntegerVector anchors2,
                             int nrow,
                             bool curved = true,
                             int density = 30) {

  int n_cols = anchors1.size() - 1;
  if (anchors2.size() != anchors1.size()) {
    stop("Rail 1 and Rail 2 must have the same number of control points for manual mode.");
  }

  List out_list(nrow * n_cols);
  int idx = 0;
  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");
  int n_steps = curved ? density : 1;

  for (int i = 0; i < n_cols; i++) {
    int idx_r1_start = anchors1[i];
    int idx_r1_end   = anchors1[i+1];

    int idx_r2_start = anchors2[i];
    int idx_r2_end   = anchors2[i+1];

    double diff_r1 = (double)(idx_r1_end - idx_r1_start);
    double diff_r2 = (double)(idx_r2_end - idx_r2_start);

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      NumericMatrix ring(2 * (n_steps + 1) + 1, 2);
      int pt_idx = 0;
      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;

        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);

        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_top) + x_r2 * u_top;
        ring(pt_idx, 1) = y_r1 * (1 - u_top) + y_r2 * u_top;
        pt_idx++;
      }

      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);
        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }
        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);
        ring(pt_idx, 0) = x_r1 * (1 - u_bot) + x_r2 * u_bot;
        ring(pt_idx, 1) = y_r1 * (1 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }
      ring(pt_idx, 0) = ring(0, 0);
      ring(pt_idx, 1) = ring(0, 1);

      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;
      out_list[idx++] = polygon_sfg;
    }
  }
  return out_list;
}

// [[Rcpp::export]]
NumericMatrix cpp_shapefile_measures(List geoms) {
  int n_geoms = geoms.size();
  NumericMatrix out(n_geoms, 6);
  double* out_ptr = REAL(out);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_geoms; i++) {
    SEXP g_sexp = geoms[i];
    if (TYPEOF(g_sexp) != VECSXP || LENGTH(g_sexp) == 0) {
      continue;
    }
    SEXP ring_sexp = VECTOR_ELT(g_sexp, 0);
    if (TYPEOF(ring_sexp) != REALSXP) {
      continue;
    }
    int n = Rf_nrows(ring_sexp);
    int cols = Rf_ncols(ring_sexp);
    if (n < 3 || cols < 2) {
      continue;
    }

    const double* r_ptr = REAL(ring_sexp);

    double area2 = 0.0;
    double cx_sum = 0.0;
    double cy_sum = 0.0;
    double perim = 0.0;

    double min_x = r_ptr[0], max_x = r_ptr[0];
    double min_y = r_ptr[n], max_y = r_ptr[n];

    for (int k = 0; k < n - 1; k++) {
      double x0 = r_ptr[k];
      double y0 = r_ptr[k + n];
      double x1 = r_ptr[k + 1];
      double y1 = r_ptr[k + 1 + n];

      if (x0 < min_x) min_x = x0;
      if (x0 > max_x) max_x = x0;
      if (y0 < min_y) min_y = y0;
      if (y0 > max_y) max_y = y0;

      double cross = (x0 * y1 - x1 * y0);
      area2 += cross;
      cx_sum += (x0 + x1) * cross;
      cy_sum += (y0 + y1) * cross;

      double dx = x1 - x0;
      double dy = y1 - y0;
      perim += std::sqrt(dx * dx + dy * dy);
    }

    double area = std::abs(area2) / 2.0;
    double cx = 0.0;
    double cy = 0.0;

    if (std::abs(area2) > 1e-12) {
      cx = cx_sum / (3.0 * area2);
      cy = cy_sum / (3.0 * area2);
    } else {
      double sum_x = 0.0, sum_y = 0.0;
      for (int k = 0; k < n - 1; k++) {
        sum_x += r_ptr[k];
        sum_y += r_ptr[k + n];
      }
      cx = sum_x / (n - 1);
      cy = sum_y / (n - 1);
    }

    double width_val = 0.0;
    double height_val = 0.0;

    if (n == 5) {
      double x1 = r_ptr[0], y1 = r_ptr[n];
      double x2 = r_ptr[1], y2 = r_ptr[1 + n];
      double x4 = r_ptr[3], y4 = r_ptr[3 + n];

      double side_A = dist_eucl(x1, y1, x2, y2);
      double side_B = dist_eucl(x1, y1, x4, y4);

      double angle_B_rad = std::atan2(y4 - y1, x4 - x1);
      double cos_angle_B = std::abs(std::cos(angle_B_rad));

      if (cos_angle_B > 0.707) {
        width_val  = std::round(side_B * 1000.0) / 1000.0;
        height_val = std::round(side_A * 1000.0) / 1000.0;
      } else {
        width_val  = std::round(side_A * 1000.0) / 1000.0;
        height_val = std::round(side_B * 1000.0) / 1000.0;
      }
    } else {
      width_val  = std::round((max_x - min_x) * 1000.0) / 1000.0;
      height_val = std::round((max_y - min_y) * 1000.0) / 1000.0;
    }

    out_ptr[i]               = cx;
    out_ptr[i + n_geoms]     = cy;
    out_ptr[i + 2 * n_geoms] = area;
    out_ptr[i + 3 * n_geoms] = perim;
    out_ptr[i + 4 * n_geoms] = width_val;
    out_ptr[i + 5 * n_geoms] = height_val;
  }

  return out;
}
