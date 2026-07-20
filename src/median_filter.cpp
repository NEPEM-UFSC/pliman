// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// =============================================================================
// Filtro de Mediana Hiper-Otimizado (Huang 2-Level Coarse-Fine Histogram)
//
// Otimizações:
//   1. Algoritmo de Huang (O(N*radius)) em vez de Perreault-Hébert (O(N*256)).
//      Para raios típicos (r=3, 5, 9), Huang faz até 20x MENOS atualizações no hist!
//   2. Loop Unswitching completo: removemos todos os 'clamping' (min/max) 
//      do miolo da imagem.
//   3. Histograma Coarse-Fine (O(32) busca de mediana).
//   4. OpenMP multi-threading por canal de cor.
// =============================================================================

// [[Rcpp::export]]
NumericVector median_filter_cpp(NumericVector img, int nrow, int ncol,
                                int nch, int radius) {
  int N = nrow * ncol;
  NumericVector out(img.size());

  const double* in_d = img.begin();
  double*       ou_d = out.begin();

  // 1. Pré-conversão vetorizada e contígua para uint8 (reduz largura de banda)
  std::vector<uint8_t> img8((size_t)N * nch);
  for (int i = 0; i < N * nch; i++) {
    img8[i] = (uint8_t)(in_d[i] * 255.0 + 0.5);
  }

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
  #endif
  for (int ch = 0; ch < nch; ch++) {
    const uint8_t* img_ch = img8.data() + (size_t)ch * N;
    double* out_ch = ou_d + (size_t)ch * N;

    // Histograma de janela thread-local (Fine: 256 bins, Coarse: 16 bins)
    std::vector<int> wfine(256, 0);
    std::vector<int> wcoarse(16, 0);

    for (int r = 0; r < nrow; r++) {
      std::fill(wfine.begin(), wfine.end(), 0);
      std::fill(wcoarse.begin(), wcoarse.end(), 0);
      int win_count = 0;

      int r_start = std::max(0, r - radius);
      int r_end = std::min(nrow - 1, r + radius);

      // Lambda para busca rápida da mediana (O(32) operações)
      auto get_median = [&]() -> double {
        int target = (win_count + 1) >> 1;
        int cum = 0, k = 0;
        while (k < 15 && cum + wcoarse[k] < target) {
          cum += wcoarse[k++];
        }
        int base = k << 4;
        for (int b = base; b < base + 16; b++) {
          cum += wfine[b];
          if (cum >= target) return b / 255.0;
        }
        return 1.0;
      };

      // Se a imagem for extremamente estreita, usa o fallback com clamping
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
          out_ch[(size_t)c * nrow + r] = get_median();
        }
        continue;
      }

      // =======================================================================
      // ALGORITMO OTIMIZADO: Sem Clamping no Corpo Principal
      // =======================================================================
      
      // Inicializa janela para c = 0 (colunas -radius até radius)
      // Colunas negativas são replicadas da coluna 0
      for (int c = -radius; c <= radius; c++) {
        int clamped_c = std::max(0, c);
        const uint8_t* col_ptr = img_ch + (size_t)clamped_c * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_ptr[k];
          wfine[v]++; wcoarse[v >> 4]++; win_count++;
        }
      }

      // --- 1. Fronteira Esquerda (c = 0 até radius) ---
      for (int c = 0; c <= radius; c++) {
        if (c > 0) {
          // col_rem < 0 -> sempre mapeia para coluna 0
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
        out_ch[(size_t)c * nrow + r] = get_median();
      }

      // --- 2. Corpo Principal (c = radius + 1 até ncol - radius - 1) ---
      // NENHUM IF, NENHUM MIN/MAX AQUI! Velocidade máxima de cache.
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
        out_ch[(size_t)c * nrow + r] = get_median();
      }

      // --- 3. Fronteira Direita (c = ncol - radius até ncol - 1) ---
      for (int c = ncol - radius; c < ncol; c++) {
        int col_rem = c - radius - 1;
        const uint8_t* col_rem_ptr = img_ch + (size_t)col_rem * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_rem_ptr[k];
          wfine[v]--; wcoarse[v >> 4]--; win_count--;
        }
        // col_add >= ncol -> sempre mapeia para a última coluna (ncol - 1)
        const uint8_t* col_add_ptr = img_ch + (size_t)(ncol - 1) * nrow;
        for (int k = r_start; k <= r_end; k++) {
          uint8_t v = col_add_ptr[k];
          wfine[v]++; wcoarse[v >> 4]++; win_count++;
        }
        out_ch[(size_t)c * nrow + r] = get_median();
      }
    }
  }

  out.attr("dim") = img.attr("dim");
  return out;
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

      // Fallback para imagens extremamente estreitas
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

      // Inicializa janela para c = 0 (colunas -radius até radius)
      for (int c = -radius; c <= radius; c++) {
        int clamped_c = std::max(0, c);
        const int* col_ptr = img_ch + (size_t)clamped_c * nrow;
        for (int k = r_start; k <= r_end; k++) {
          if (col_ptr[k]) win_ones++;
          win_count++;
        }
      }

      // --- 1. Fronteira Esquerda (c = 0 até radius) ---
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

      // --- 2. Corpo Principal (c = radius + 1 até ncol - radius - 1) ---
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

      // --- 3. Fronteira Direita (c = ncol - radius até ncol - 1) ---
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

