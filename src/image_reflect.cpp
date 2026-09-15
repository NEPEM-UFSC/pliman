// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cstring>
#include <type_traits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static inline void parse_dims(const IntegerVector& dims,
                               int& nrow, int& ncol, int& nch,
                               const char* fname) {
  if (dims.size() == 2) {
    nrow = dims[0]; ncol = dims[1]; nch = 1;
  } else if (dims.size() == 3) {
    nrow = dims[0]; ncol = dims[1]; nch = dims[2];
  } else {
    Rcpp::stop("%s: input must be a 2-D matrix or 3-D array.", fname);
  }
}

template <int RTYPE>
SEXP do_hreflect(SEXP img_sexp, const IntegerVector& dims, int nrow, int ncol, int nch) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  
  VecType img(img_sexp);
  VecType out = Rcpp::no_init(img.size());

  const T* __restrict__ src = img.begin();
  T*       __restrict__ dst = out.begin();

  const std::ptrdiff_t plane_size = (std::ptrdiff_t)nrow * ncol;
  const int total_cols            = nch * ncol;
  const int nrow_m1               = nrow - 1;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int abs_c = 0; abs_c < total_cols; ++abs_c) {
    const int ch = abs_c / ncol;
    const int c  = abs_c % ncol;

    const T* __restrict__ col_src =
        src + (std::ptrdiff_t)ch * plane_size + (std::ptrdiff_t)c * nrow;
    T*       __restrict__ col_dst =
        dst + (std::ptrdiff_t)ch * plane_size + (std::ptrdiff_t)c * nrow;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int r = 0; r < nrow; ++r) {
      col_dst[r] = col_src[nrow_m1 - r];
    }
  }

  out.attr("dim") = Rcpp::clone(dims);
  return out;
}

// [[Rcpp::export]]
SEXP image_hreflect_cpp(SEXP img_sexp) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);

  int nrow, ncol, nch;
  parse_dims(dims, nrow, ncol, nch, "image_hreflect_cpp");

  if (TYPEOF(img_sexp) == RAWSXP) {
    return do_hreflect<RAWSXP>(img_sexp, dims, nrow, ncol, nch);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    return do_hreflect<REALSXP>(img_sexp, dims, nrow, ncol, nch);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    return do_hreflect<INTSXP>(img_sexp, dims, nrow, ncol, nch);
  } else {
    Rcpp::stop("image_hreflect_cpp: unsupported input type.");
    return R_NilValue;
  }
}

template <int RTYPE>
SEXP do_vreflect(SEXP img_sexp, const IntegerVector& dims, int nrow, int ncol, int nch) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;
  
  VecType img(img_sexp);
  VecType out = Rcpp::no_init(img.size());

  const T* __restrict__ src = img.begin();
  T*       __restrict__ dst = out.begin();

  const std::size_t  col_bytes  = (std::size_t)nrow * sizeof(T);
  const std::ptrdiff_t plane_size = (std::ptrdiff_t)nrow * ncol;
  const int total_cols            = nch * ncol;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int abs_c = 0; abs_c < total_cols; ++abs_c) {
    const int ch = abs_c / ncol;
    const int c  = abs_c % ncol;

    const int dst_c = (ncol - 1) - c;

    const T* __restrict__ col_src =
        src + (std::ptrdiff_t)ch * plane_size + (std::ptrdiff_t)c * nrow;
    T*       __restrict__ col_dst =
        dst + (std::ptrdiff_t)ch * plane_size + (std::ptrdiff_t)dst_c * nrow;

    std::memcpy(col_dst, col_src, col_bytes);
  }

  out.attr("dim") = Rcpp::clone(dims);
  return out;
}

// [[Rcpp::export]]
SEXP image_vreflect_cpp(SEXP img_sexp) {
  IntegerVector dims = Rf_getAttrib(img_sexp, R_DimSymbol);

  int nrow, ncol, nch;
  parse_dims(dims, nrow, ncol, nch, "image_vreflect_cpp");

  if (TYPEOF(img_sexp) == RAWSXP) {
    return do_vreflect<RAWSXP>(img_sexp, dims, nrow, ncol, nch);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    return do_vreflect<REALSXP>(img_sexp, dims, nrow, ncol, nch);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    return do_vreflect<INTSXP>(img_sexp, dims, nrow, ncol, nch);
  } else {
    Rcpp::stop("image_vreflect_cpp: unsupported input type.");
    return R_NilValue;
  }
}
