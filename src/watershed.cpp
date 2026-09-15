// watershed.cpp — Watershed morfológico ultra-otimizado (Hyper Engine)
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Vec8 = std::vector<uint8_t>;
using VecI = std::vector<int>;

static constexpr int WS_INF_INT = 1000000000;
static constexpr float WS_INF_FLOAT = 1e15f;

// Tabela estática de consulta (LUT) para raiz quadrada inteira rápida
static float g_ws_sqrt_lut[262145];
static bool g_ws_lut_init = []() {
  for (int i = 0; i <= 262144; i++) {
    g_ws_sqrt_lut[i] = std::sqrt(static_cast<float>(i));
  }
  return true;
}();

static inline float ws_fast_sqrt(int d2) {
  if (d2 <= 262144) return g_ws_sqrt_lut[d2];
  return std::sqrt(static_cast<float>(d2));
}

// 1-D Felzenszwalb Distance Transform com buffer na pilha (Zero alocação no heap)
static inline void ws_dt1d_fast(const int* f_in, int* f_out_base, int out_stride, int* v, float* z, int n) {
  int k = 0;
  v[0] = 0; z[0] = -WS_INF_FLOAT; z[1] = WS_INF_FLOAT;
  for (int q = 1; q < n; q++) {
    int fq = f_in[q];
    float s;
    do {
      int vk = v[k];
      s = static_cast<float>((fq + q * q) - (f_in[vk] + vk * vk)) / static_cast<float>(2 * (q - vk));
      if (s > z[k]) break;
      --k;
    } while (true);
    ++k; v[k] = q; z[k] = s; z[k + 1] = WS_INF_FLOAT;
  }
  k = 0;
  for (int q = 0; q < n; q++) {
    while (z[k + 1] < static_cast<float>(q)) ++k;
    int vk = v[k];
    f_out_base[q * out_stride] = (q - vk) * (q - vk) + f_in[vk];
  }
}

// 2-D EDT Ultra-Rápida para grids globais
static VecI ws_edt2d_fused(const uint8_t* p_bin, int nrow, int ncol) {
  int N = nrow * ncol;
  VecI D_t(N);

  #pragma omp parallel for schedule(static)
  for (int c = 0; c < ncol; c++) {
    const uint8_t* col_in = p_bin + c * nrow;
    int* row_out = D_t.data() + c;

    int temp_d[4096];
    int* p_td = (nrow <= 4096) ? temp_d : (new int[nrow]);

    int d = 100000;
    for (int r = 0; r < nrow; r++) {
      if (col_in[r] == 0) d = 0;
      else if (d < 100000) d++;
      p_td[r] = d;
    }
    d = 100000;
    for (int r = nrow - 1; r >= 0; r--) {
      if (col_in[r] == 0) d = 0;
      else if (d < 100000) d++;
      if (d < p_td[r]) p_td[r] = d;
      int v = p_td[r];
      row_out[r * ncol] = (v >= 100000) ? WS_INF_INT : v * v;
    }
    if (nrow > 4096) delete[] p_td;
  }

  VecI D(N);
  #pragma omp parallel
  {
    int v_buf[4096];
    float z_buf[4096];
    int* v = (ncol + 2 <= 4096) ? v_buf : (new int[ncol + 2]);
    float* z = (ncol + 2 <= 4096) ? z_buf : (new float[ncol + 2]);

    #pragma omp for schedule(static)
    for (int r = 0; r < nrow; r++) {
      const int* p_r = D_t.data() + r * ncol;
      ws_dt1d_fast(p_r, D.data() + r, nrow, v, z, ncol);
    }

    if (ncol + 2 > 4096) { delete[] v; delete[] z; }
  }

  return D;
}

// Fast inlined Union-Find com path compression
static inline int fast_find_root(int i, int* parent) {
  int root = i;
  while (root != parent[root]) root = parent[root];
  int curr = i;
  while (curr != root) {
    int nxt = parent[curr];
    parent[curr] = root;
    curr = nxt;
  }
  return root;
}

// Pool de memória scratch para evitar malloc/free repetidos em cada thread
struct WSThreadScratch {
  std::vector<int> D_vec;
  std::vector<int> D_t_vec;
  std::vector<int> td_vec;
  std::vector<int> v_vec;
  std::vector<float> z_vec;
  std::vector<int> head;
  std::vector<int> next;
  std::vector<int> label;
  std::vector<int> parent;
  std::vector<float> seed_peak_val;
  std::vector<int> seed_peak_r;
  std::vector<int> seed_peak_c;

  void ensure_capacity(int s_N, int s_nrow, int s_ncol, int max_d2) {
    if ((int)D_vec.size() < s_N) D_vec.resize(s_N * 2);
    if ((int)D_t_vec.size() < s_N) D_t_vec.resize(s_N * 2);
    if ((int)td_vec.size() < s_nrow) td_vec.resize(s_nrow * 2);
    if ((int)v_vec.size() < s_ncol + 2) v_vec.resize((s_ncol + 2) * 2);
    if ((int)z_vec.size() < s_ncol + 2) z_vec.resize((s_ncol + 2) * 2);
    if ((int)head.size() < max_d2 + 1) head.resize((max_d2 + 1) * 2, -1);
    if ((int)next.size() < s_N) next.resize(s_N * 2, -1);
    if ((int)label.size() < s_N) label.resize(s_N * 2, 0);
  }
};

// Watershed Regional dentro de uma sub-caixa recortada
static inline int ws_scratch_sub_watershed(const uint8_t* sub_bin, int* sub_out, int s_nrow, int s_ncol, double tolerance, int ext, WSThreadScratch& sc) {
  int s_N = s_nrow * s_ncol;
  sc.ensure_capacity(s_N, s_nrow, s_ncol, 4096);

  int* D = sc.D_vec.data();
  int* D_t = sc.D_t_vec.data();
  int* p_td = sc.td_vec.data();
  int* v_buf = sc.v_vec.data();
  float* z_buf = sc.z_vec.data();

  // 2D EDT no sub-box
  for (int c = 0; c < s_ncol; c++) {
    const uint8_t* col_in = sub_bin + c * s_nrow;
    int* row_out = D_t + c;
    int d = 100000;
    for (int r = 0; r < s_nrow; r++) {
      if (col_in[r] == 0) d = 0;
      else if (d < 100000) d++;
      p_td[r] = d;
    }
    d = 100000;
    for (int r = s_nrow - 1; r >= 0; r--) {
      if (col_in[r] == 0) d = 0;
      else if (d < 100000) d++;
      if (d < p_td[r]) p_td[r] = d;
      int v = p_td[r];
      row_out[r * s_ncol] = (v >= 100000) ? WS_INF_INT : v * v;
    }
  }
  for (int r = 0; r < s_nrow; r++) {
    const int* p_r = D_t + r * s_ncol;
    ws_dt1d_fast(p_r, D + r, s_nrow, v_buf, z_buf, s_ncol);
  }

  int max_d2 = 0;
  for (int i = 0; i < s_N; i++) {
    if (!sub_bin[i]) D[i] = 0;
    else if (D[i] > max_d2) max_d2 = D[i];
  }

  if (max_d2 == 0) return 0;
  sc.ensure_capacity(s_N, s_nrow, s_ncol, max_d2);

  int* head = sc.head.data();
  int* next = sc.next.data();
  std::fill(head, head + max_d2 + 1, -1);
  std::fill(next, next + s_N, -1);

  for (int i = 0; i < s_N; i++) {
    int d2 = D[i];
    if (d2 > 0) {
      next[i] = head[d2];
      head[d2] = i;
    }
  }

  int* p_label = sc.label.data();
  std::fill(p_label, p_label + s_N, 0);

  auto& parent = sc.parent;
  auto& seed_peak_val = sc.seed_peak_val;
  auto& seed_peak_r = sc.seed_peak_r;
  auto& seed_peak_c = sc.seed_peak_c;

  parent.clear();
  seed_peak_val.clear();
  seed_peak_r.clear();
  seed_peak_c.clear();

  parent.push_back(0);
  seed_peak_val.push_back(0.0f);
  seed_peak_r.push_back(0);
  seed_peak_c.push_back(0);

  int n_seeds = 0;
  bool use_relative = (tolerance < 1.0);
  float tol_f = static_cast<float>(tolerance);
  int* p_parent = parent.data();

  const int doff_8[8] = {-1, 1, -s_nrow, s_nrow, -s_nrow - 1, -s_nrow + 1, s_nrow - 1, s_nrow + 1};
  const int dr_8[8]   = {-1, 1,  0,       0,      -1,          1,          -1,         1};
  const int dc_8[8]   = { 0, 0, -1,       1,      -1,         -1,           1,         1};

  for (int d2 = max_d2; d2 >= 1; d2--) {
    int p_idx = head[d2];
    if (p_idx == -1) continue;
    float v = ws_fast_sqrt(d2);

    while (p_idx != -1) {
      int p_c = p_idx / s_nrow;
      int p_r = p_idx - p_c * s_nrow;
      bool is_interior = (p_r >= 1 && p_r < s_nrow - 1 && p_c >= 1 && p_c < s_ncol - 1);

      int n_roots[8];
      int n_roots_sz = 0;

      if (is_interior) {
        const int* pl = p_label + p_idx;
        int l0 = pl[-1], l1 = pl[1], l2 = pl[-s_nrow], l3 = pl[s_nrow];
        int l4 = pl[-s_nrow - 1], l5 = pl[-s_nrow + 1], l6 = pl[s_nrow - 1], l7 = pl[s_nrow + 1];

        int first_lab = (l0 > 0) ? l0 : ((l1 > 0) ? l1 : ((l2 > 0) ? l2 : ((l3 > 0) ? l3 : ((l4 > 0) ? l4 : ((l5 > 0) ? l5 : ((l6 > 0) ? l6 : l7))))));

        if (first_lab > 0) {
          int r_first = (p_parent[first_lab] == first_lab) ? first_lab : fast_find_root(first_lab, p_parent);
          bool all_same = true;
          #define WSH_CHK(l) if (l > 0 && ((p_parent[l] == l) ? l : fast_find_root(l, p_parent)) != r_first) { all_same = false; }
          WSH_CHK(l0); WSH_CHK(l1); WSH_CHK(l2); WSH_CHK(l3);
          WSH_CHK(l4); WSH_CHK(l5); WSH_CHK(l6); WSH_CHK(l7);
          #undef WSH_CHK

          if (all_same) {
            p_label[p_idx] = r_first;
            p_idx = next[p_idx];
            continue;
          }

          n_roots[0] = r_first;
          n_roots_sz = 1;
          #define WSH_AD(l) if (l > 0) { int r = (p_parent[l] == l) ? l : fast_find_root(l, p_parent); bool f = false; for (int k=0; k<n_roots_sz; k++) if (n_roots[k]==r) { f=true; break; } if (!f) n_roots[n_roots_sz++] = r; }
          WSH_AD(l0); WSH_AD(l1); WSH_AD(l2); WSH_AD(l3);
          WSH_AD(l4); WSH_AD(l5); WSH_AD(l6); WSH_AD(l7);
          #undef WSH_AD
        }
      } else {
        for (int i = 0; i < 8; i++) {
          int nr = p_r + dr_8[i];
          int nc = p_c + dc_8[i];
          if (static_cast<unsigned>(nr) < static_cast<unsigned>(s_nrow) &&
              static_cast<unsigned>(nc) < static_cast<unsigned>(s_ncol)) {
            int l = p_label[p_idx + doff_8[i]];
            if (l > 0) {
              int r = (p_parent[l] == l) ? l : fast_find_root(l, p_parent);
              bool f = false;
              for (int k = 0; k < n_roots_sz; k++) {
                if (n_roots[k] == r) { f = true; break; }
              }
              if (!f) n_roots[n_roots_sz++] = r;
            }
          }
        }
      }

      if (n_roots_sz == 0) {
        n_seeds++;
        p_label[p_idx] = n_seeds;
        parent.push_back(n_seeds);
        seed_peak_val.push_back(v);
        seed_peak_r.push_back(p_r);
        seed_peak_c.push_back(p_c);
        p_parent = parent.data();

      } else if (n_roots_sz == 1) {
        p_label[p_idx] = n_roots[0];

      } else {
        int highest_root = n_roots[0];
        float max_peak = seed_peak_val[highest_root];
        for (int k = 1; k < n_roots_sz; k++) {
          int r = n_roots[k];
          if (seed_peak_val[r] > max_peak) {
            max_peak = seed_peak_val[r];
            highest_root = r;
          }
        }

        for (int k = 0; k < n_roots_sz; k++) {
          int root = n_roots[k];
          if (root == highest_root) continue;
          float diff = seed_peak_val[root] - v;
          float thresh = use_relative ? (tol_f * (seed_peak_val[root] + 1e-6f)) : tol_f;
          if (diff < thresh) {
            p_parent[root] = highest_root;
          }
        }

        int res = (p_parent[highest_root] == highest_root) ? highest_root : fast_find_root(highest_root, p_parent);
        float min_dist_spatial = 1e15f;

        for (int k = 0; k < n_roots_sz; k++) {
          int root = fast_find_root(n_roots[k], p_parent);
          int pr = seed_peak_r[root];
          int pc = seed_peak_c[root];
          float dist_spatial = static_cast<float>((pr - p_r)*(pr - p_r) + (pc - p_c)*(pc - p_c));
          if (dist_spatial < min_dist_spatial) {
            min_dist_spatial = dist_spatial;
            res = root;
          }
        }

        p_label[p_idx] = res;
      }

      p_idx = next[p_idx];
    }
  }

  for (size_t i = 1; i < parent.size(); i++) {
    parent[i] = fast_find_root(static_cast<int>(i), p_parent);
  }

  std::vector<int> final_label(parent.size(), 0);
  int current_label = 0;
  for (size_t i = 1; i < parent.size(); i++) {
    if (parent[i] == static_cast<int>(i)) {
      final_label[i] = ++current_label;
    }
  }
  for (size_t i = 1; i < parent.size(); i++) {
    final_label[i] = final_label[parent[i]];
  }

  for (int i = 0; i < s_N; i++) {
    int l = p_label[i];
    sub_out[i] = (l > 0) ? final_label[l] : 0;
  }

  return current_label;
}

// Watershed de Grid Completo (usado como fallback para componentes gigantes > 40% da imagem)
static IntegerMatrix watershed_full_grid_internal(SEXP img_r, double tolerance, int ext) {
  int nrow = Rf_nrows(img_r);
  int ncol = Rf_ncols(img_r);
  int N = nrow * ncol;

  Vec8 img(N);
  if (TYPEOF(img_r) == RAWSXP) {
    const uint8_t* ptr = RAW(img_r);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) img[i] = ptr[i] ? 1u : 0u;
  } else if (TYPEOF(img_r) == LGLSXP || TYPEOF(img_r) == INTSXP) {
    const int* ptr = INTEGER(img_r);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) img[i] = ptr[i] ? 1u : 0u;
  } else if (TYPEOF(img_r) == REALSXP) {
    const double* ptr = REAL(img_r);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) img[i] = (ptr[i] != 0.0) ? 1u : 0u;
  }

  VecI D = ws_edt2d_fused(img.data(), nrow, ncol);

  int max_d2 = 0;
  #pragma omp parallel for reduction(max:max_d2) schedule(static)
  for (int i = 0; i < N; i++) {
    if (!img[i]) {
      D[i] = 0;
    } else if (D[i] > max_d2) {
      max_d2 = D[i];
    }
  }

  if (max_d2 == 0) return IntegerMatrix(nrow, ncol);

  VecI head(max_d2 + 1, -1);
  VecI next(N, -1);
  for (int i = 0; i < N; i++) {
    int d2 = D[i];
    if (d2 > 0) {
      next[i] = head[d2];
      head[d2] = i;
    }
  }

  std::vector<int> label(N, 0);
  int* p_label = label.data();
  std::vector<int> parent;
  parent.reserve(32768);
  std::vector<float> seed_peak_val;
  std::vector<int> seed_peak_r, seed_peak_c;
  seed_peak_val.reserve(32768);
  seed_peak_r.reserve(32768);
  seed_peak_c.reserve(32768);

  parent.push_back(0);
  seed_peak_val.push_back(0.0f);
  seed_peak_r.push_back(0);
  seed_peak_c.push_back(0);
  int n_seeds = 0;

  const int doff_8[8] = {-1, 1, -nrow, nrow, -nrow - 1, -nrow + 1, nrow - 1, nrow + 1};
  const int dr_8[8]   = {-1, 1,  0,    0,    -1,        1,         -1,       1};
  const int dc_8[8]   = { 0, 0, -1,    1,    -1,       -1,          1,       1};

  bool use_relative = (tolerance < 1.0);
  float tol_f = static_cast<float>(tolerance);
  int* p_parent = parent.data();

  for (int d2 = max_d2; d2 >= 1; d2--) {
    int p_idx = head[d2];
    if (p_idx == -1) continue;
    float v = ws_fast_sqrt(d2);

    while (p_idx != -1) {
      int p_c = p_idx / nrow;
      int p_r = p_idx - p_c * nrow;
      bool is_interior = (p_r >= 1 && p_r < nrow - 1 && p_c >= 1 && p_c < ncol - 1);

      int n_roots[8];
      int n_roots_sz = 0;

      if (is_interior) {
        const int* pl = p_label + p_idx;
        int l0 = pl[-1], l1 = pl[1], l2 = pl[-nrow], l3 = pl[nrow];
        int l4 = pl[-nrow - 1], l5 = pl[-nrow + 1], l6 = pl[nrow - 1], l7 = pl[nrow + 1];

        int first_lab = (l0 > 0) ? l0 : ((l1 > 0) ? l1 : ((l2 > 0) ? l2 : ((l3 > 0) ? l3 : ((l4 > 0) ? l4 : ((l5 > 0) ? l5 : ((l6 > 0) ? l6 : l7))))));

        if (first_lab > 0) {
          int r_first = (p_parent[first_lab] == first_lab) ? first_lab : fast_find_root(first_lab, p_parent);
          bool all_same = true;
          #define WSH_CHK_G(l) if (l > 0 && ((p_parent[l] == l) ? l : fast_find_root(l, p_parent)) != r_first) { all_same = false; }
          WSH_CHK_G(l0); WSH_CHK_G(l1); WSH_CHK_G(l2); WSH_CHK_G(l3);
          WSH_CHK_G(l4); WSH_CHK_G(l5); WSH_CHK_G(l6); WSH_CHK_G(l7);
          #undef WSH_CHK_G

          if (all_same) {
            p_label[p_idx] = r_first;
            p_idx = next[p_idx];
            continue;
          }

          n_roots[0] = r_first;
          n_roots_sz = 1;
          #define WSH_AD_G(l) if (l > 0) { int r = (p_parent[l] == l) ? l : fast_find_root(l, p_parent); bool f = false; for (int k=0; k<n_roots_sz; k++) if (n_roots[k]==r) { f=true; break; } if (!f) n_roots[n_roots_sz++] = r; }
          WSH_AD_G(l0); WSH_AD_G(l1); WSH_AD_G(l2); WSH_AD_G(l3);
          WSH_AD_G(l4); WSH_AD_G(l5); WSH_AD_G(l6); WSH_AD_G(l7);
          #undef WSH_AD_G
        }
      } else {
        for (int i = 0; i < 8; i++) {
          int nr = p_r + dr_8[i];
          int nc = p_c + dc_8[i];
          if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
              static_cast<unsigned>(nc) < static_cast<unsigned>(ncol)) {
            int l = p_label[p_idx + doff_8[i]];
            if (l > 0) {
              int r = (p_parent[l] == l) ? l : fast_find_root(l, p_parent);
              bool f = false;
              for (int k = 0; k < n_roots_sz; k++) {
                if (n_roots[k] == r) { f = true; break; }
              }
              if (!f) n_roots[n_roots_sz++] = r;
            }
          }
        }
      }

      if (n_roots_sz == 0) {
        n_seeds++;
        p_label[p_idx] = n_seeds;
        parent.push_back(n_seeds);
        seed_peak_val.push_back(v);
        seed_peak_r.push_back(p_r);
        seed_peak_c.push_back(p_c);
        p_parent = parent.data();

      } else if (n_roots_sz == 1) {
        p_label[p_idx] = n_roots[0];

      } else {
        int highest_root = n_roots[0];
        float max_peak = seed_peak_val[highest_root];
        for (int k = 1; k < n_roots_sz; k++) {
          int r = n_roots[k];
          if (seed_peak_val[r] > max_peak) {
            max_peak = seed_peak_val[r];
            highest_root = r;
          }
        }

        for (int k = 0; k < n_roots_sz; k++) {
          int root = n_roots[k];
          if (root == highest_root) continue;
          float diff = seed_peak_val[root] - v;
          float thresh = use_relative ? (tol_f * (seed_peak_val[root] + 1e-6f)) : tol_f;
          if (diff < thresh) {
            p_parent[root] = highest_root;
          }
        }

        int res = (p_parent[highest_root] == highest_root) ? highest_root : fast_find_root(highest_root, p_parent);
        float min_dist_spatial = 1e15f;

        for (int k = 0; k < n_roots_sz; k++) {
          int root = fast_find_root(n_roots[k], p_parent);
          int pr = seed_peak_r[root];
          int pc = seed_peak_c[root];
          float dist_spatial = static_cast<float>((pr - p_r)*(pr - p_r) + (pc - p_c)*(pc - p_c));
          if (dist_spatial < min_dist_spatial) {
            min_dist_spatial = dist_spatial;
            res = root;
          }
        }

        p_label[p_idx] = res;
      }

      p_idx = next[p_idx];
    }
  }

  for (size_t i = 1; i < parent.size(); i++) {
    parent[i] = fast_find_root(static_cast<int>(i), p_parent);
  }

  std::vector<int> final_label(parent.size(), 0);
  int current_label = 0;
  for (size_t i = 1; i < parent.size(); i++) {
    if (parent[i] == static_cast<int>(i)) {
      final_label[i] = ++current_label;
    }
  }
  for (size_t i = 1; i < parent.size(); i++) {
    final_label[i] = final_label[parent[i]];
  }

  IntegerMatrix out(nrow, ncol);
  int* p_out = INTEGER(out);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    int l = p_label[i];
    p_out[i] = (l > 0) ? final_label[l] : 0;
  }

  return out;
}

// -----------------------------------------------------------------------------
// watershed_cpp: Motor Watershed Principal e Ultra-Otimizado do pliman
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix watershed_cpp(SEXP img_r, double tolerance = 1.0, int ext = 1) {
  int nrow = Rf_nrows(img_r);
  int ncol = Rf_ncols(img_r);
  int N = nrow * ncol;

  const int* p_bin_int = nullptr;
  const uint8_t* p_bin_raw = nullptr;
  const double* p_bin_real = nullptr;
  if (TYPEOF(img_r) == LGLSXP || TYPEOF(img_r) == INTSXP) {
    p_bin_int = INTEGER(img_r);
  } else if (TYPEOF(img_r) == RAWSXP) {
    p_bin_raw = RAW(img_r);
  } else if (TYPEOF(img_r) == REALSXP) {
    p_bin_real = REAL(img_r);
  }

  auto is_fg = [&](int idx) -> bool {
    if (p_bin_int) return p_bin_int[idx] != 0;
    if (p_bin_raw) return p_bin_raw[idx] != 0;
    if (p_bin_real) return p_bin_real[idx] != 0.0 && !std::isnan(p_bin_real[idx]);
    return false;
  };

  struct Run {
    int r0, r1, c;
    int label;
  };

  std::vector<Run> runs;
  runs.reserve(N / 32);
  std::vector<int> col_run_start(ncol + 1, 0);

  for (int c = 0; c < ncol; c++) {
    col_run_start[c] = (int)runs.size();
    int c_off = c * nrow;
    int r = 0;
    while (r < nrow) {
      while (r < nrow && !is_fg(c_off + r)) r++;
      if (r >= nrow) break;
      int r0 = r;
      while (r < nrow && is_fg(c_off + r)) r++;
      runs.push_back({r0, r - 1, c, 0});
    }
  }
  col_run_start[ncol] = (int)runs.size();

  int num_runs = (int)runs.size();
  if (num_runs == 0) {
    return IntegerMatrix(nrow, ncol);
  }

  std::vector<int> parent(num_runs + 1);
  for (int i = 1; i <= num_runs; i++) parent[i] = i;

  for (int i = 0; i < num_runs; i++) {
    runs[i].label = i + 1;
  }

  for (int c = 1; c < ncol; c++) {
    int start_prev = col_run_start[c - 1];
    int end_prev   = col_run_start[c];
    int start_curr = col_run_start[c];
    int end_curr   = col_run_start[c + 1];

    if (start_prev == end_prev || start_curr == end_curr) continue;

    int p = start_prev;
    for (int q = start_curr; q < end_curr; q++) {
      int qr0 = runs[q].r0 - 1;
      int qr1 = runs[q].r1 + 1;

      while (p < end_prev && runs[p].r1 < qr0) p++;

      int cur_p = p;
      while (cur_p < end_prev && runs[cur_p].r0 <= qr1) {
        int root1 = fast_find_root(runs[q].label, parent.data());
        int root2 = fast_find_root(runs[cur_p].label, parent.data());
        if (root1 != root2) {
          parent[root2] = root1;
        }
        cur_p++;
      }
    }
  }

  std::vector<int> comp_r0(num_runs + 1, nrow), comp_r1(num_runs + 1, -1);
  std::vector<int> comp_c0(num_runs + 1, ncol), comp_c1(num_runs + 1, -1);
  std::vector<int> comp_area(num_runs + 1, 0);

  int* p_parent = parent.data();
  for (int i = 0; i < num_runs; i++) {
    int root = fast_find_root(runs[i].label, p_parent);
    runs[i].label = root;

    int r0 = runs[i].r0, r1 = runs[i].r1, c = runs[i].c;
    if (r0 < comp_r0[root]) comp_r0[root] = r0;
    if (r1 > comp_r1[root]) comp_r1[root] = r1;
    if (c < comp_c0[root])  comp_c0[root] = c;
    if (c > comp_c1[root])  comp_c1[root] = c;
    comp_area[root] += (r1 - r0 + 1);
  }

  std::vector<int> active_comps;
  active_comps.reserve(num_runs);
  for (int i = 1; i <= num_runs; i++) {
    if (comp_area[i] > 0) active_comps.push_back(i);
  }

  int num_comps = (int)active_comps.size();

  long long total_bbox_area = 0;
  int max_single_bbox = 0;
  for (int ci = 0; ci < num_comps; ci++) {
    int id = active_comps[ci];
    int s_area = (comp_r1[id] - comp_r0[id] + 1) * (comp_c1[id] - comp_c0[id] + 1);
    total_bbox_area += s_area;
    if (s_area > max_single_bbox) max_single_bbox = s_area;
  }

  if (max_single_bbox > N * 0.3 || total_bbox_area > N * 0.7 || num_comps > 400 || (N > 2000000 && num_comps > 200)) {
    return watershed_full_grid_internal(img_r, tolerance, ext);
  }

  std::vector<std::vector<int>> comp_run_indices(num_runs + 1);
  for (int i = 0; i < num_runs; i++) {
    comp_run_indices[runs[i].label].push_back(i);
  }

  struct SubRes {
    int r0, c0, s_nrow, s_ncol;
    int num_sub_labels;
    std::vector<int> sub_labels;
  };

  std::vector<SubRes> sub_res(num_comps);

  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_max_threads();
  #endif
  std::vector<WSThreadScratch> thread_scratch(max_threads);

  #pragma omp parallel for schedule(dynamic)
  for (int ci = 0; ci < num_comps; ci++) {
    int id = active_comps[ci];
    int r0 = std::max(0, comp_r0[id] - 1);
    int r1 = std::min(nrow - 1, comp_r1[id] + 1);
    int c0 = std::max(0, comp_c0[id] - 1);
    int c1 = std::min(ncol - 1, comp_c1[id] + 1);

    int s_nrow = r1 - r0 + 1;
    int s_ncol = c1 - c0 + 1;
    int s_N = s_nrow * s_ncol;

    auto& res = sub_res[ci];
    res.r0 = r0;
    res.c0 = c0;
    res.s_nrow = s_nrow;
    res.s_ncol = s_ncol;
    res.sub_labels.resize(s_N, 0);

    const auto& r_indices = comp_run_indices[id];

    if (comp_area[id] < 16) {
      res.num_sub_labels = 1;
      for (int ri : r_indices) {
        int sc = runs[ri].c - c0;
        int sc_off = sc * s_nrow;
        for (int sr = runs[ri].r0 - r0; sr <= runs[ri].r1 - r0; sr++) {
          res.sub_labels[sc_off + sr] = 1;
        }
      }
      continue;
    }

    std::vector<uint8_t> sub_bin(s_N, 0);
    for (int ri : r_indices) {
      int sc = runs[ri].c - c0;
      int sc_off = sc * s_nrow;
      for (int sr = runs[ri].r0 - r0; sr <= runs[ri].r1 - r0; sr++) {
        sub_bin[sc_off + sr] = 1;
      }
    }

    int tid = 0;
    #ifdef _OPENMP
      tid = omp_get_thread_num();
    #endif

    res.num_sub_labels = ws_scratch_sub_watershed(sub_bin.data(), res.sub_labels.data(), s_nrow, s_ncol, tolerance, ext, thread_scratch[tid]);
  }

  IntegerMatrix out(nrow, ncol);
  int* p_out = INTEGER(out);

  int global_label_offset = 0;
  for (int ci = 0; ci < num_comps; ci++) {
    const auto& res = sub_res[ci];
    if (res.num_sub_labels == 0) continue;

    int offset = global_label_offset;
    for (int sc = 0; sc < res.s_ncol; sc++) {
      int gc = res.c0 + sc;
      int gc_off = gc * nrow;
      int sc_off = sc * res.s_nrow;
      for (int sr = 0; sr < res.s_nrow; sr++) {
        int gr = res.r0 + sr;
        int l = res.sub_labels[sc_off + sr];
        if (l > 0) {
          p_out[gc_off + gr] = l + offset;
        }
      }
    }
    global_label_offset += res.num_sub_labels;
  }

  return out;
}



// [[Rcpp::export]]
NumericMatrix help_dist_transform(LogicalMatrix bin) {
  int nrow = bin.nrow();
  int ncol = bin.ncol();
  int N = nrow * ncol;

  Vec8 img(N);
  const int* p_bin = LOGICAL(bin);
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) {
    img[i] = p_bin[i] ? 1u : 0u;
  }

  VecI dist_sq = ws_edt2d_fused(img.data(), nrow, ncol);

  NumericMatrix out(nrow, ncol);
  double* p_out = REAL(out);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    p_out[i] = p_bin[i] ? std::sqrt(static_cast<double>(dist_sq[i])) : 0.0;
  }

  return out;
}
