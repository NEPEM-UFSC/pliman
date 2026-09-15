// index_aggregator.cpp — Fast per-object index mean computation
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// ─────────────────────────────────────────────────────────────────────────────
// compute_index_means_cpp
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix compute_index_means_cpp(const NumericMatrix& idx_mat,
                                      SEXP                 labels_sexp,
                                      const IntegerVector& valid_ids) {

  IntegerVector labels_flat;
  switch (TYPEOF(labels_sexp)) {
    case INTSXP:  labels_flat = as<IntegerVector>(labels_sexp); break;
    case REALSXP: labels_flat = as<IntegerVector>(NumericVector(labels_sexp)); break;
    default:
      Rcpp::stop("compute_index_means_cpp: 'labels' must be integer or numeric.");
  }

  const int N     = labels_flat.size();
  const int n_idx = idx_mat.ncol();
  const int n_obj = valid_ids.size();

  if (idx_mat.nrow() != N)
    Rcpp::stop("compute_index_means_cpp: nrow(idx_mat) must equal length(labels_flat).");

  int max_id = 0;
  for (int i = 0; i < n_obj; ++i)
    if (valid_ids[i] > max_id) max_id = valid_ids[i];

  std::vector<int> id_map(max_id + 1, -1);
  for (int i = 0; i < n_obj; ++i)
    id_map[valid_ids[i]] = i;

  std::vector<double> sums((std::size_t)n_obj * n_idx, 0.0);
  std::vector<int>    counts(n_obj, 0);

  const int* lp = labels_flat.begin();

  for (int k = 0; k < N; ++k) {
    const int lab = lp[k];
    if (lab > 0 && lab <= max_id) {
      const int pos = id_map[lab];
      if (pos >= 0) ++counts[pos];
    }
  }

  for (int j = 0; j < n_idx; ++j) {
    const double* col = idx_mat.begin() + (std::ptrdiff_t)j * N;
    double*       acc = sums.data()     + (std::ptrdiff_t)j * n_obj;

    for (int k = 0; k < N; ++k) {
      const int lab = lp[k];
      if (lab > 0 && lab <= max_id) {
        const int pos = id_map[lab];
        if (pos >= 0) acc[pos] += col[k];
      }
    }
  }

  NumericMatrix result(n_obj, n_idx);
  for (int j = 0; j < n_idx; ++j) {
    const double* acc = sums.data() + (std::ptrdiff_t)j * n_obj;
    for (int i = 0; i < n_obj; ++i) {
      result(i, j) = counts[i] > 0 ? acc[i] / counts[i] : NA_REAL;
    }
  }

  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// rgb_to_hsb_cpp 
// ─────────────────────────────────────────────────────────────────────────────

template <int RTYPE>
List do_rgb_to_hsb(SEXP R_sexp, SEXP G_sexp, SEXP B_sexp, const IntegerVector& dims, int N) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  VecType r_vec(R_sexp), g_vec(G_sexp), b_vec(B_sexp);
  const T* rp = r_vec.begin();
  const T* gp = g_vec.begin();
  const T* bp = b_vec.begin();

  NumericVector H(N), S(N), Br(N);

  for (int k = 0; k < N; ++k) {
    double r, g, b;
    if constexpr (std::is_same_v<T, Rbyte> || std::is_same_v<T, int>) {
      r = (double)rp[k] / 255.0;
      g = (double)gp[k] / 255.0;
      b = (double)bp[k] / 255.0;
    } else {
      r = (double)rp[k];
      g = (double)gp[k];
      b = (double)bp[k];
    }
    const double mx = r > g ? (r > b ? r : b) : (g > b ? g : b);
    const double mn = r < g ? (r < b ? r : b) : (g < b ? g : b);
    const double diff = mx - mn;

    Br[k] = mx * 100.0;
    S[k]  = mx > 0.0 ? diff / mx * 100.0 : 0.0;

    if (diff == 0.0) {
      H[k] = 0.0;
    } else if (mx == r) {
      H[k] = 60.0 * ((g - b) / diff);
    } else if (mx == g) {
      H[k] = 60.0 * (2.0 + (b - r) / diff);
    } else {
      H[k] = 60.0 * (4.0 + (r - g) / diff);
    }
    if (H[k] < 0.0) H[k] += 360.0;
  }

  H.attr("dim")  = Rcpp::clone(dims);
  S.attr("dim")  = Rcpp::clone(dims);
  Br.attr("dim") = Rcpp::clone(dims);

  return List::create(_["H"] = H, _["S"] = S, _["B"] = Br);
}

// [[Rcpp::export]]
List rgb_to_hsb_cpp(SEXP R_sexp, SEXP G_sexp, SEXP B_sexp) {
  IntegerVector dims = Rf_getAttrib(R_sexp, R_DimSymbol);
  int N = dims.size() > 1 ? dims[0] * dims[1] : Rf_length(R_sexp);

  if (TYPEOF(R_sexp) == RAWSXP) {
    return do_rgb_to_hsb<RAWSXP>(R_sexp, G_sexp, B_sexp, dims, N);
  } else if (TYPEOF(R_sexp) == REALSXP) {
    return do_rgb_to_hsb<REALSXP>(R_sexp, G_sexp, B_sexp, dims, N);
  } else if (TYPEOF(R_sexp) == INTSXP) {
    return do_rgb_to_hsb<INTSXP>(R_sexp, G_sexp, B_sexp, dims, N);
  } else {
    Rcpp::stop("rgb_to_hsb_cpp: unsupported input type.");
    return R_NilValue;
  }
}
