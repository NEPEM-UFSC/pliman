#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// Extração de Contorno Altamente Otimizada
// =============================================================================
// Melhorias:
// 1. Acesso à matriz cache-friendly (col-major scan em vez de row-major).
// 2. Uso de ponteiro 1D direto (`INTEGER(labels)`) evitando overhead de Rcpp::Matrix.
// 3. Verificação de limites unificada (unsigned cast) e remoção de módulo (`% 8`).
// 4. Pré-reserva de memória (`reserve`) para evitar realocações durante o traçado.
//
// [[Rcpp::export]]
List extract_contours_cpp(IntegerMatrix labels) {
  int nrow = labels.nrow(), ncol = labels.ncol();
  int N = nrow * ncol;
  const int* ptr = INTEGER(labels);

  int max_label = 0;
  for (int i = 0; i < N; i++) {
    if (ptr[i] > max_label) max_label = ptr[i];
  }

  std::vector<int> start_r(max_label + 1, -1);
  std::vector<int> start_c(max_label + 1, -1);
  std::vector<int> start_idx(max_label + 1, -1);

  // Scan O(N) Cache-Friendly (Column-Major)
  int idx = 0;
  for (int c = 0; c < ncol; c++) {
    for (int r = 0; r < nrow; r++, idx++) {
      int id = ptr[idx];
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

  List out_list(max_label);
  CharacterVector col_names = CharacterVector::create("x", "y");

  for (int id = 1; id <= max_label; id++) {
    if (start_idx[id] == -1) continue;

    int sr = start_r[id];
    int sc = start_c[id];
    int sidx = start_idx[id];

    std::vector<int> b_x;
    std::vector<int> b_y;
    b_x.reserve(512); // Previne múltiplas realocações
    b_y.reserve(512);

    b_x.push_back(sr);
    b_y.push_back(sc);

    int curr_r = sr;
    int curr_c = sc;
    int curr_idx = sidx;
    int backtrack = 6; // Background garantido a oeste devido ao col-major scan

    int next_r = -1, next_c = -1, next_idx = -1, next_backtrack = -1;
    bool found = false;

    for (int i = 1; i <= 8; i++) {
      int dir = backtrack + i;
      if (dir >= 8) dir -= 8;

      int nr = curr_r + dr[dir];
      int nc = curr_c + dc[dir];

      if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
          static_cast<unsigned>(nc) < static_cast<unsigned>(ncol)) {
        int nidx = curr_idx + doff[dir];
        if (ptr[nidx] == id) {
          next_r = nr;
          next_c = nc;
          next_idx = nidx;
          next_backtrack = dir + 4;
          if (next_backtrack >= 8) next_backtrack -= 8;
          found = true;
          break;
        }
      }
    }

    if (found && (next_idx != sidx)) {
      int second_idx = next_idx;

      curr_r = next_r;
      curr_c = next_c;
      curr_idx = next_idx;
      backtrack = next_backtrack;

      b_x.push_back(curr_r);
      b_y.push_back(curr_c);

      int max_iter = nrow * ncol;
      int iter = 0;

      while (iter < max_iter) {
        iter++;
        bool step_found = false;

        for (int i = 1; i <= 8; i++) {
          int dir = backtrack + i;
          if (dir >= 8) dir -= 8;

          int nr = curr_r + dr[dir];
          int nc = curr_c + dc[dir];

          if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
              static_cast<unsigned>(nc) < static_cast<unsigned>(ncol)) {
            int nidx = curr_idx + doff[dir];
            if (ptr[nidx] == id) {
              next_r = nr;
              next_c = nc;
              next_idx = nidx;
              next_backtrack = dir + 4;
              if (next_backtrack >= 8) next_backtrack -= 8;
              step_found = true;
              break;
            }
          }
        }

        if (!step_found) break;

        if (curr_idx == sidx && next_idx == second_idx) {
          break;
        }

        curr_r = next_r;
        curr_c = next_c;
        curr_idx = next_idx;
        backtrack = next_backtrack;

        b_x.push_back(curr_r);
        b_y.push_back(curr_c);
      }
    }

    int n_points = b_x.size();
    if (n_points > 0) {
      IntegerMatrix mat(n_points, 2);
      int* mat_ptr = INTEGER(mat);
      for (int j = 0; j < n_points; j++) {
        mat_ptr[j] = b_x[j] + 2;            // x = 1st dim + 1
        mat_ptr[j + n_points] = b_y[j] + 2; // y = 2nd dim + 1
      }
      colnames(mat) = col_names;
      out_list[id - 1] = mat;
    }
  }

  return out_list;
}
