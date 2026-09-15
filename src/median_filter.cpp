// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// =============================================================================
// Filtro de Mediana Hiper-Otimizado (Huang 2-Level Coarse-Fine Histogram)
// =============================================================================

template <int RTYPE>
SEXP do_median_filter(SEXP img_sexp, int nrow, int ncol, int nch, int radius) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  int N = nrow * ncol;
  
  std::vector<uint8_t> buffer;
  const uint8_t* in_ptr;
  if constexpr (std::is_same_v<T, Rbyte>) {
    in_ptr = RAW(img_sexp);
  } else if constexpr (std::is_same_v<T, double>) {
    buffer.resize((size_t)N * nch);
    double* dptr = REAL(img_sexp);
    for (int i = 0; i < N * nch; i++) buffer[i] = (uint8_t)(dptr[i] * 255.0 + 0.5);
    in_ptr = buffer.data();
  } else {
    buffer.resize((size_t)N * nch);
    int* iptr = INTEGER(img_sexp);
    for (int i = 0; i < N * nch; i++) {
      buffer[i] = (uint8_t)(std::min(255, std::max(0, iptr[i])));
    }
    in_ptr = buffer.data();
  }

  VecType out = Rcpp::no_init(N * nch);
  T* ou_d = out.begin();

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
  #endif
  for (int ch = 0; ch < nch; ch++) {
    const uint8_t* img_ch = in_ptr + (size_t)ch * N;
    T* out_ch = ou_d + (size_t)ch * N;

    std::vector<int> wfine(256, 0);
    std::vector<int> wcoarse(16, 0);

    for (int r = 0; r < nrow; r++) {
      std::fill(wfine.begin(), wfine.end(), 0);
      std::fill(wcoarse.begin(), wcoarse.end(), 0);
      int win_count = 0;

      int r_start = std::max(0, r - radius);
      int r_end = std::min(nrow - 1, r + radius);

      auto get_median_raw = [&]() -> int {
        int target = (win_count + 1) >> 1;
        int cum = 0, k = 0;
        while (k < 15 && cum + wcoarse[k] < target) {
          cum += wcoarse[k++];
        }
        int base = k << 4;
        for (int b = base; b < base + 16; b++) {
          cum += wfine[b];
          if (cum >= target) return b;
        }
        return 255;
      };

      auto set_out = [&](int r_idx, int c_idx) {
        int m = get_median_raw();
        if constexpr (std::is_same_v<T, Rbyte>) {
          out_ch[(size_t)c_idx * nrow + r_idx] = (T)m;
        } else if constexpr (std::is_same_v<T, double>) {
          out_ch[(size_t)c_idx * nrow + r_idx] = (T)(m / 255.0);
        } else {
          out_ch[(size_t)c_idx * nrow + r_idx] = (T)m;
        }
      };

      if (ncol <= 2 * radius) {
        for (int c = -radius; c <= radius; c++) {
          int clamped_c = std::max(0, std::min(ncol - 1, c));
          const uint8_t* col_ptr = img_ch + (size_t)clamped_c * nrow;
          for (int k = r_start; k <= r_end; k++) {
            uint8_t v = col_ptr[k];
            wfine[v]++; wcoarse[v >> 4]++; win_count++;
          }
        }
        for (int c = 0; c < ncol; c++) {
          if (c > 0) {
            int col_rem = std::max(0, std::min(ncol - 1, c - radius - 1));
            const uint8_t* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
            for (int k = r_start; k <= r_end; k++) {
              uint8_t v = col_rem_ptr[k];
              wfine[v]--; wcoarse[v >> 4]--; win_count--;
            }
            int col_add = std::max(0, std::min(ncol - 1, c + radius));
            const uint8_t* col_add_ptr = img_ch + (size_t)col_add * nrow;
            for (int k = r_start; k <= r_end; k++) {
              uint8_t v = col_add_ptr[k];
              wfine[v]++; wcoarse[v >> 4]++; win_count++;
            }
          }
          set_out(r, c);
        }
        continue;
      }
      
      for (int c = -radius; c <= radius; c++) {
        int clamped_c = std::max(0, c);
        const uint8_t* col_ptr = img_ch + (size_t)clamped_c * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_ptr[k];
          wfine[v]++; wcoarse[v >> 4]++; win_count++;
        }
      }

      for (int c = 0; c <= radius; c++) {
        if (c > 0) {
          const uint8_t* col_rem_ptr = img_ch;
          for (int k = r_start; k <= r_end; k++) {
            uint8_t v = col_rem_ptr[k];
            wfine[v]--; wcoarse[v >> 4]--; win_count--;
          }
          int col_add = c + radius;
          const uint8_t* col_add_ptr = img_ch + (size_t)col_add * nrow;
          for (int k = r_start; k <= r_end; k++) {
            uint8_t v = col_add_ptr[k];
            wfine[v]++; wcoarse[v >> 4]++; win_count++;
          }
        }
        set_out(r, c);
      }

      for (int c = radius + 1; c < ncol - radius; c++) {
        int col_rem = c - radius - 1;
        const uint8_t* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_rem_ptr[k];
          wfine[v]--; wcoarse[v >> 4]--; win_count--;
        }
        int col_add = c + radius;
        const uint8_t* col_add_ptr = img_ch + (size_t)col_add * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_add_ptr[k];
          wfine[v]++; wcoarse[v >> 4]++; win_count++;
        }
        set_out(r, c);
      }

      for (int c = ncol - radius; c < ncol; c++) {
        int col_rem = c - radius - 1;
        const uint8_t* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_rem_ptr[k];
          wfine[v]--; wcoarse[v >> 4]--; win_count--;
        }
        const uint8_t* col_add_ptr = img_ch + (size_t)(ncol - 1) * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_add_ptr[k];
          wfine[v]++; wcoarse[v >> 4]++; win_count++;
        }
        set_out(r, c);
      }
    }
  }

  out.attr("dim") = Rf_getAttrib(img_sexp, R_DimSymbol);
  return out;
}

// [[Rcpp::export]]
SEXP median_filter_cpp(SEXP img_sexp, int nrow, int ncol,
                       int nch, int radius) {
  if (TYPEOF(img_sexp) == RAWSXP) {
    return do_median_filter<RAWSXP>(img_sexp, nrow, ncol, nch, radius);
  } else if (TYPEOF(img_sexp) == REALSXP) {
    return do_median_filter<REALSXP>(img_sexp, nrow, ncol, nch, radius);
  } else if (TYPEOF(img_sexp) == INTSXP) {
    return do_median_filter<INTSXP>(img_sexp, nrow, ncol, nch, radius);
  } else {
    Rcpp::stop("median_filter_cpp: unsupported input type.");
    return R_NilValue;
  }
}

// =============================================================================
// Filtro de Mediana Binária (Majority Filter) - Hiper-Otimizado para Logical
// =============================================================================
// [[Rcpp::export]]
LogicalVector median_filter_binary_cpp(SEXP img_sexp, int nrow, int ncol,
                                       int nch, int radius) {
  int* img = LOGICAL(img_sexp); 
  int N = nrow * ncol;
  LogicalVector out(N * nch);
  int* out_ptr = LOGICAL(out);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
  #endif
  for (int ch = 0; ch < nch; ch++) {
    const int* img_ch = img + (size_t)ch * N;
    int* out_ch = out_ptr + (size_t)ch * N;

    for (int r = 0; r < nrow; r++) {
      int win_ones = 0;
      int win_count = 0;

      int r_start = std::max(0, r - radius);
      int r_end = std::min(nrow - 1, r + radius);

      if (ncol <= 2 * radius) {
        for (int c = -radius; c <= radius; c++) {
          int clamped_c = std::max(0, std::min(ncol - 1, c));
          const int* col_ptr = img_ch + (size_t)clamped_c * nrow;
          for (int k = r_start; k <= r_end; k++) {
            if (col_ptr[k]) win_ones++;
            win_count++;
          }
        }
        for (int c = 0; c < ncol; c++) {
          if (c > 0) {
            int col_rem = std::max(0, std::min(ncol - 1, c - radius - 1));
            const int* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
            for (int k = r_start; k <= r_end; k++) {
              if (col_rem_ptr[k]) win_ones--;
              win_count--;
            }
            int col_add = std::max(0, std::min(ncol - 1, c + radius));
            const int* col_add_ptr = img_ch + (size_t)col_add * nrow;
            for (int k = r_start; k <= r_end; k++) {
              if (col_add_ptr[k]) win_ones++;
              win_count++;
            }
          }
          int target = (win_count + 1) >> 1;
          out_ch[(size_t)c * nrow + r] = (win_ones >= target) ? 1 : 0;
        }
        continue;
      }

      for (int c = -radius; c <= radius; c++) {
        int clamped_c = std::max(0, c);
        const int* col_ptr = img_ch + (size_t)clamped_c * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_ptr[k]) win_ones++;
          win_count++;
        }
      }

      for (int c = 0; c <= radius; c++) {
        if (c > 0) {
          const int* col_rem_ptr = img_ch;
          for (int k = r_start; k <= r_end; k++) {
            if (col_rem_ptr[k]) win_ones--;
            win_count--;
          }
          int col_add = c + radius;
          const int* col_add_ptr = img_ch + (size_t)col_add * nrow;
          for (int k = r_start; k <= r_end; k++) {
            if (col_add_ptr[k]) win_ones++;
            win_count++;
          }
        }
        int target = (win_count + 1) >> 1;
        out_ch[(size_t)c * nrow + r] = (win_ones >= target) ? 1 : 0;
      }

      for (int c = radius + 1; c < ncol - radius; c++) {
        int col_rem = c - radius - 1;
        const int* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_rem_ptr[k]) win_ones--;
          win_count--;
        }
        int col_add = c + radius;
        const int* col_add_ptr = img_ch + (size_t)col_add * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_add_ptr[k]) win_ones++;
          win_count++;
        }
        int target = (win_count + 1) >> 1;
        out_ch[(size_t)c * nrow + r] = (win_ones >= target) ? 1 : 0;
      }

      for (int c = ncol - radius; c < ncol; c++) {
        int col_rem = c - radius - 1;
        const int* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_rem_ptr[k]) win_ones--;
          win_count--;
        }
        const int* col_add_ptr = img_ch + (size_t)(ncol - 1) * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_add_ptr[k]) win_ones++;
          win_count++;
        }
        int target = (win_count + 1) >> 1;
        out_ch[(size_t)c * nrow + r] = (win_ones >= target) ? 1 : 0;
      }
    }
  }

  out.attr("dim") = RObject(img_sexp).attr("dim");
  return out;
}
