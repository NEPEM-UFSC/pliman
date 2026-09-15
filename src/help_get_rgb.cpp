// help_get_rgb.cpp — Optimized pixel extraction by label
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <vector>
#include <type_traits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static inline IntegerMatrix coerce_labels(SEXP labels_sexp, const char* caller) {
  switch (TYPEOF(labels_sexp)) {
    case INTSXP:  return as<IntegerMatrix>(labels_sexp);
    case REALSXP: return as<IntegerMatrix>(NumericMatrix(labels_sexp));
    default:
      Rcpp::stop("%s: 'labels' must be an integer or numeric matrix (got %s).",
                 caller, Rf_type2char(TYPEOF(labels_sexp)));
  }
}

template <int RTYPE>
std::vector<std::vector<double>> do_help_get_rgb(SEXP R_sexp, SEXP G_sexp, SEXP B_sexp,
                                                 const IntegerMatrix& labels) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  VecType r_vec(R_sexp), g_vec(G_sexp), b_vec(B_sexp);
  const T* rp = r_vec.begin();
  const T* gp = g_vec.begin();
  const T* bp = b_vec.begin();

  const int nrow = labels.nrow();
  const int ncol = labels.ncol();
  const int N    = nrow * ncol;
  const int* lp  = labels.begin();

  int max_label = 0;
  for (int i = 0; i < N; ++i)
    if (lp[i] > max_label) max_label = lp[i];

  const int nlab = max_label + 1;
  std::vector<int> cnt(nlab, 0);
  for (int k = 0; k < N; ++k)
    if (lp[k] > 0) ++cnt[lp[k]];

  std::vector<std::vector<double>> result(nlab);
  for (int lab = 1; lab < nlab; ++lab)
    result[lab].reserve((std::size_t)cnt[lab] * 4);

  for (int j = 0; j < ncol; ++j) {
    const int base = j * nrow;
    for (int i = 0; i < nrow; ++i) {
      const int k   = base + i;
      const int lab = lp[k];
      if (lab > 0) {
        auto& v = result[lab];
        v.push_back(static_cast<double>(lab));
        v.push_back(static_cast<double>(rp[k]));
        v.push_back(static_cast<double>(gp[k]));
        v.push_back(static_cast<double>(bp[k]));
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
help_get_rgb(SEXP R_sexp, SEXP G_sexp, SEXP B_sexp, SEXP labels_sexp) {
  const IntegerMatrix labels = coerce_labels(labels_sexp, "help_get_rgb");
  if (TYPEOF(R_sexp) == RAWSXP) {
    return do_help_get_rgb<RAWSXP>(R_sexp, G_sexp, B_sexp, labels);
  } else if (TYPEOF(R_sexp) == REALSXP) {
    return do_help_get_rgb<REALSXP>(R_sexp, G_sexp, B_sexp, labels);
  } else if (TYPEOF(R_sexp) == INTSXP) {
    return do_help_get_rgb<INTSXP>(R_sexp, G_sexp, B_sexp, labels);
  } else {
    Rcpp::stop("help_get_rgb: unsupported image type.");
    return std::vector<std::vector<double>>();
  }
}

template <int RTYPE>
std::vector<std::vector<double>> do_help_get_renir(SEXP RE_sexp, SEXP NIR_sexp,
                                                   const IntegerMatrix& labels) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  VecType re_vec(RE_sexp), nir_vec(NIR_sexp);
  const T* rep = re_vec.begin();
  const T* nip = nir_vec.begin();

  const int nrow = labels.nrow();
  const int ncol = labels.ncol();
  const int N    = nrow * ncol;
  const int* lp  = labels.begin();

  int max_label = 0;
  for (int i = 0; i < N; ++i)
    if (lp[i] > max_label) max_label = lp[i];

  const int nlab = max_label + 1;
  std::vector<int> cnt(nlab, 0);
  for (int k = 0; k < N; ++k)
    if (lp[k] > 0) ++cnt[lp[k]];

  std::vector<std::vector<double>> result(nlab);
  for (int lab = 1; lab < nlab; ++lab)
    result[lab].reserve((std::size_t)cnt[lab] * 3);

  for (int j = 0; j < ncol; ++j) {
    const int base = j * nrow;
    for (int i = 0; i < nrow; ++i) {
      const int k   = base + i;
      const int lab = lp[k];
      if (lab > 0) {
        auto& v = result[lab];
        v.push_back(static_cast<double>(lab));
        v.push_back(static_cast<double>(rep[k]));
        v.push_back(static_cast<double>(nip[k]));
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
help_get_renir(SEXP RE_sexp, SEXP NIR_sexp, SEXP labels_sexp) {
  const IntegerMatrix labels = coerce_labels(labels_sexp, "help_get_renir");
  if (TYPEOF(RE_sexp) == RAWSXP) {
    return do_help_get_renir<RAWSXP>(RE_sexp, NIR_sexp, labels);
  } else if (TYPEOF(RE_sexp) == REALSXP) {
    return do_help_get_renir<REALSXP>(RE_sexp, NIR_sexp, labels);
  } else if (TYPEOF(RE_sexp) == INTSXP) {
    return do_help_get_renir<INTSXP>(RE_sexp, NIR_sexp, labels);
  } else {
    Rcpp::stop("help_get_renir: unsupported image type.");
    return std::vector<std::vector<double>>();
  }
}
