// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// =============================================================================
// color_labels_cpp: Coloriza uma matriz de rótulos inteiros (labels) em RGB
//
// Usa o número áureo para gerar cores infinitamente distintas e vibrantes.
// Background (0) é mapeado para Preto (0, 0, 0).
// =============================================================================
// [[Rcpp::export]]
NumericVector color_labels_cpp(SEXP labels_sexp) {
  int nrow = Rf_nrows(labels_sexp);
  int ncol = Rf_ncols(labels_sexp);
  int N = nrow * ncol;

  NumericVector out(N * 3);
  double* out_ptr = out.begin();

  // 1. Encontra o maior label para criar uma Tabela de Cores (Color Lookup Table)
  int max_l = 0;
  
  if (TYPEOF(labels_sexp) == REALSXP) {
    const double* p_labels = REAL(labels_sexp);
    for (int i = 0; i < N; i++) {
      int val = (int)p_labels[i];
      if (val > max_l) max_l = val;
    }
  } else if (TYPEOF(labels_sexp) == INTSXP || TYPEOF(labels_sexp) == LGLSXP) {
    const int* p_labels = INTEGER(labels_sexp);
    for (int i = 0; i < N; i++) {
      int val = p_labels[i];
      if (val > max_l) max_l = val;
    }
  }

  // Se não houver nenhum objeto, retorna tudo preto
  if (max_l <= 0) {
    out.attr("dim") = IntegerVector::create(nrow, ncol, 3);
    return out;
  }

  // 2. Prepara a tabela de cores RGB
  std::vector<double> r_table(max_l + 1, 0.0);
  std::vector<double> g_table(max_l + 1, 0.0);
  std::vector<double> b_table(max_l + 1, 0.0);

  for (int l = 1; l <= max_l; l++) {
    // Distribuição ultra-vibrante usando a Proporção Áurea
    double h = std::fmod(l * 0.618033988749895, 1.0);
    double s = 0.85; // Saturação ideal
    double v = 0.95; // Brilho ideal

    int i = std::floor(h * 6.0);
    double f = h * 6.0 - i;
    double p = v * (1.0 - s);
    double q = v * (1.0 - f * s);
    double t = v * (1.0 - (1.0 - f) * s);

    double r_val = 0.0, g_val = 0.0, b_val = 0.0;
    switch (i % 6) {
      case 0: r_val = v;     g_val = t;     b_val = p;     break;
      case 1: r_val = q;     g_val = v;     b_val = p;     break;
      case 2: r_val = p;     g_val = v;     b_val = t;     break;
      case 3: r_val = p;     g_val = q;     b_val = v;     break;
      case 4: r_val = t;     g_val = p;     b_val = v;     break;
      case 5: r_val = v;     g_val = p;     b_val = q;     break;
    }
    r_table[l] = r_val;
    g_table[l] = g_val;
    b_table[l] = b_val;
  }

  double* out_r = out_ptr;
  double* out_g = out_ptr + N;
  double* out_b = out_ptr + 2 * (size_t)N;

  // 3. Mapeamento paralelo via OpenMP
  if (TYPEOF(labels_sexp) == REALSXP) {
    const double* p_labels = REAL(labels_sexp);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < N; i++) {
      int l = (int)p_labels[i];
      if (l > 0 && l <= max_l) {
        out_r[i] = r_table[l];
        out_g[i] = g_table[l];
        out_b[i] = b_table[l];
      } else {
        out_r[i] = 0.0;
        out_g[i] = 0.0;
        out_b[i] = 0.0;
      }
    }
  } else {
    const int* p_labels = INTEGER(labels_sexp);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < N; i++) {
      int l = p_labels[i];
      if (l > 0 && l <= max_l) {
        out_r[i] = r_table[l];
        out_g[i] = g_table[l];
        out_b[i] = b_table[l];
      } else {
        out_r[i] = 0.0;
        out_g[i] = 0.0;
        out_b[i] = 0.0;
      }
    }
  }

  out.attr("dim") = IntegerVector::create(nrow, ncol, 3);
  return out;
}
