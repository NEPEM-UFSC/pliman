// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

static inline int find_root(int x, int* p) {
  while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
  return x;
}

// =============================================================================
// bwlabel_cpp: CCL RLE Column-Major com Union-Find + OpenMP
//
// Estratégia de alocação:
//   Passe 0 (gratuito, O(N)):  Conta os runs exatos por coluna. Zero alocação.
//   Passe 1 (RLE + UF):        Alocação EXATA baseada na contagem real.
//   Passe 2 (OpenMP render):   Preenchimento paralelo por coluna.
// =============================================================================
// [[Rcpp::export]]
IntegerMatrix bwlabel_cpp(SEXP img_sexp) {
  int* img = INTEGER(img_sexp);
  int nrow = Rf_nrows(img_sexp);
  int ncol = Rf_ncols(img_sexp);

  // ---------------------------------------------------------------------------
  // PASSE 0: Contagem exata de runs (zero alocação, apenas leitura)
  // ---------------------------------------------------------------------------
  std::vector<int> col_rs(ncol + 1, 0);
  int total_runs = 0;

  for (int c = 0; c < ncol; c++) {
    const int* col = img + (size_t)c * nrow;
    int cnt = 0;
    int r = 0;
    while (r < nrow) {
      if (!col[r]) { r++; continue; }
      while (r < nrow && col[r]) r++;
      cnt++;
    }
    col_rs[c] = total_runs;
    total_runs += cnt;
  }
  col_rs[ncol] = total_runs;

  // ---------------------------------------------------------------------------
  // Alocação EXATA (sem desperdício)
  // R_alloc = stack-like, liberado pelo GC sem custo de destruição
  // ---------------------------------------------------------------------------
  int* rs     = (int*) R_alloc(total_runs + 1, sizeof(int)); // run start
  int* re_arr = (int*) R_alloc(total_runs + 1, sizeof(int)); // run end
  int* rl     = (int*) R_alloc(total_runs + 1, sizeof(int)); // run label (temp)
  std::memset(rl, 0, (total_runs + 1) * sizeof(int));

  // Parent UF: no máximo total_runs componentes distintos
  int* parent = (int*) R_alloc(total_runs + 2, sizeof(int));
  parent[0] = 0;
  int next_label = 1;

  // ---------------------------------------------------------------------------
  // PASSE 1: RLE scan + Union-Find (usa col_rs já calculado)
  // ---------------------------------------------------------------------------
  for (int c = 0; c < ncol; c++) {
    const int* col = img + (size_t)c * nrow;

    int prev_start = (c > 0) ? col_rs[c - 1] : 0;
    int prev_end   = col_rs[c];
    int prev_j     = prev_start;
    int nr         = prev_end; // índice de escrita para esta coluna

    int r = 0;
    while (r < nrow) {
      if (!col[r]) { r++; continue; }

      int s = r;
      while (r < nrow && col[r]) r++;
      int e = r - 1;

      rs[nr]     = s;
      re_arr[nr] = e;

      // Avança runs da coluna anterior que ficaram acima deste
      while (prev_j < prev_end && re_arr[prev_j] < s - 1) prev_j++;

      // Union-Find com runs adjacentes (8-connectivity)
      int m = 0;
      for (int j = prev_j; j < prev_end && rs[j] <= e + 1; j++) {
        int l = rl[j];
        if (!l) continue;
        int root_l = find_root(l, parent);
        if (!m) {
          m = root_l;
        } else if (root_l != m) {
          if (root_l < m) { parent[m] = root_l; m = root_l; }
          else             { parent[root_l] = m; }
        }
      }

      if (!m) {
        parent[next_label] = next_label;
        m = next_label++;
      }
      rl[nr++] = m;
    }
  }

  // ---------------------------------------------------------------------------
  // PASSE 2: Flatten Tree em O(labels)
  // ---------------------------------------------------------------------------
  // Reutiliza 'rs' como lookup table (não precisa mais dos starts durante render)
  // Na verdade usa um vetor menor:
  std::vector<int> final_label(next_label, 0);
  int current_label = 0;

  for (int i = 1; i < next_label; i++) {
    int root = find_root(i, parent);
    if (!final_label[root]) final_label[root] = ++current_label;
    final_label[i] = final_label[root];
  }

  // Aplica flatten aos rótulos de cada run (O(total_runs))
  for (int i = 0; i < total_runs; i++) {
    rl[i] = final_label[rl[i]];
  }

  // ---------------------------------------------------------------------------
  // PASSE 3: Renderização paralela por coluna (sem race conditions)
  // ---------------------------------------------------------------------------
  IntegerMatrix out(nrow, ncol);
  int* p_out = INTEGER(out);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 64)
  #endif
  for (int c = 0; c < ncol; c++) {
    int* dst    = p_out + (size_t)c * nrow;
    int  r_from = col_rs[c];
    int  r_to   = col_rs[c + 1];

    // Zera a coluna (memset é vetorizado nativamente)
    std::memset(dst, 0, nrow * sizeof(int));

    for (int i = r_from; i < r_to; i++) {
      int l   = rl[i];
      int len = re_arr[i] - rs[i] + 1;
      std::fill(dst + rs[i], dst + rs[i] + len, l);
    }
  }

  return out;
}
