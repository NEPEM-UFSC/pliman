// watershed.cpp — Watershed morfológico otimizado via EDT e Union-Find
// Algoritmo inspirado no EBImage, mas centenas de vezes mais rápido
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using Vec8 = std::vector<uint8_t>;
using VecF = std::vector<float>;
using VecI = std::vector<int>;

static constexpr float WS_INF = 1e15f;

// =============================================================================
// 1-D Felzenszwalb Distance Transform
// =============================================================================
static void ws_dt1d(float* f, float* tmp, int* v, float* z, int n) {
  int k = 0;
  v[0] = 0; z[0] = -WS_INF; z[1] = WS_INF;
  for (int q = 1; q < n; q++) {
    float fq = f[q], s;
    do {
      int vk = v[k];
      s = ((fq + (float)(q*q)) - (f[vk] + (float)(vk*vk))) / (float)(2*(q-vk));
      if (s > z[k]) break;
      --k;
    } while (true);
    ++k; v[k] = q; z[k] = s; z[k+1] = WS_INF;
  }
  k = 0;
  for (int q = 0; q < n; q++) {
    while (z[k+1] < (float)q) ++k;
    int vk = v[k];
    tmp[q] = (float)(q-vk)*(q-vk) + f[vk];
  }
  for (int i = 0; i < n; i++) f[i] = tmp[i];
}

// =============================================================================
// 2-D EDT (column-major)
// =============================================================================
static VecF ws_edt2d_real(const Vec8& img, int nrow, int ncol) {
  int N = nrow * ncol;
  VecF D(N, WS_INF);
  for (int i = 0; i < N; i++) if (!img[i]) D[i] = 0.0f;

  #pragma omp parallel
  {
    VecF tmp_c(nrow); VecI vc(nrow + 1); VecF zc(nrow + 2);
    #pragma omp for schedule(static)
    for (int c = 0; c < ncol; c++)
      ws_dt1d(D.data() + c * nrow, tmp_c.data(), vc.data(), zc.data(), nrow);

    VecF row_buf(ncol), tmp_r(ncol); VecI vr(ncol + 1); VecF zr(ncol + 2);
    #pragma omp for schedule(static)
    for (int r = 0; r < nrow; r++) {
      for (int c = 0; c < ncol; c++) row_buf[c] = D[c * nrow + r];
      ws_dt1d(row_buf.data(), tmp_r.data(), vr.data(), zr.data(), ncol);
      for (int c = 0; c < ncol; c++) D[c * nrow + r] = row_buf[c];
    }
  }

  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) {
    if (!img[i]) D[i] = 0.0f;
    else D[i] = std::sqrt(D[i]);
  }
  return D;
}

// =============================================================================
// Union-Find (Disjoint Set) com path compression (sem recursão)
// =============================================================================
inline int find_root(int i, std::vector<int>& parent) {
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

// =============================================================================
// Watershed Spanning Forest Altamente Otimizado (watershed_cpp2)
// =============================================================================
// Melhorias:
// 1. Array de structs (Pixel) e sort amigável ao cache.
// 2. Remoção de divisões inteiras (/, %) durante a resolução guardando peak_r/c.
// 3. Pré-cálculo de offsets 1D para os vizinhos e cast unificado (unsigned).
// 4. "Fast Stack Array" para n_roots.

struct Pixel {
  float d;
  int idx;
  int r;
  int c;
};

struct SeedInfo {
  float peak_val;
  int peak_r;
  int peak_c;
};

// [[Rcpp::export]]
IntegerMatrix watershed_cpp(LogicalMatrix img_r, double tolerance = 1.0, int ext = 1) {
  int nrow = img_r.nrow(), ncol = img_r.ncol();
  int N = nrow * ncol;

  Vec8 img(N);
  for (int i = 0; i < N; i++) img[i] = img_r[i] ? 1u : 0u;

  // 1. Distância Euclidiana
  VecF dist = ws_edt2d_real(img, nrow, ncol);

  // 2. Coleta de pixels em Struct cache-friendly
  std::vector<Pixel> px;
  px.reserve(N); // Reserva caso extremo (toda a imagem for foreground)
  for (int c = 0; c < ncol; c++) {
    for (int r = 0; r < nrow; r++) {
      int idx = c * nrow + r;
      if (img[idx]) {
        px.push_back({dist[idx], idx, r, c});
      }
    }
  }

  // 3. Ordenação extremamente rápida (Struct locality)
  std::sort(px.begin(), px.end(), [](const Pixel& a, const Pixel& b) {
    if (a.d != b.d) return a.d > b.d;
    return a.idx > b.idx;
  });

  // Pré-cálculo de deslocamentos
  std::vector<int> off_dr, off_dc, off_1d;
  for (int r = -ext; r <= ext; r++) {
    for (int c = -ext; c <= ext; c++) {
      if (r == 0 && c == 0) continue;
      off_dr.push_back(r);
      off_dc.push_back(c);
      off_1d.push_back(c * nrow + r);
    }
  }
  int n_offsets = off_dr.size();

  // 4. Estruturas do Union-Find e Sementes
  VecI label(N, 0);
  VecI parent;
  parent.reserve(100000);
  std::vector<SeedInfo> seed_info;
  seed_info.reserve(100000);

  parent.push_back(0); // Root 0
  seed_info.push_back({0.0f, 0, 0});
  int n_seeds = 0;

  // Espaço de stack super rápido (suporta ext até 15 sem quebrar o laço)
  const int MAX_ROOTS = 1024;

  // 5. Varredura Otimizada
  for (const Pixel& p : px) {
    float v = p.d;
    int p_r = p.r;
    int p_c = p.c;

    int n_roots[MAX_ROOTS];
    int n_roots_sz = 0;

    for (int i = 0; i < n_offsets; i++) {
      int r = p_r + off_dr[i];
      int c = p_c + off_dc[i];

      // Bounds check ultra-veloz de única instrução por eixo
      if (static_cast<unsigned>(r) < static_cast<unsigned>(nrow) &&
          static_cast<unsigned>(c) < static_cast<unsigned>(ncol)) {

        int n_idx = p.idx + off_1d[i];
        if (label[n_idx] > 0) {
          int root = find_root(label[n_idx], parent);
          bool found = false;
          for (int k = 0; k < n_roots_sz; k++) {
            if (n_roots[k] == root) { found = true; break; }
          }
          if (!found && n_roots_sz < MAX_ROOTS) {
            n_roots[n_roots_sz++] = root;
          }
        }
      }
    }

    if (n_roots_sz == 0) {
      n_seeds++;
      label[p.idx] = n_seeds;
      parent.push_back(n_seeds);
      seed_info.push_back({v, p_r, p_c});

    } else if (n_roots_sz == 1) {
      label[p.idx] = n_roots[0];

    } else {
      float maxdiff = -1.0f;
      int res = n_roots[0];
      float min_dist_spatial = 1e15f;

      for (int k = 0; k < n_roots_sz; k++) {
        int root = n_roots[k];
        float diff = seed_info[root].peak_val - v;

        if (diff > maxdiff) {
          maxdiff = diff;
          if (min_dist_spatial == 1e15f) res = root;
        }

        if (diff >= tolerance) {
          int pr = seed_info[root].peak_r;
          int pc = seed_info[root].peak_c;
          float dist_spatial = (float)((pr - p_r)*(pr - p_r) + (pc - p_c)*(pc - p_c));

          if (dist_spatial < min_dist_spatial) {
            min_dist_spatial = dist_spatial;
            res = root;
          }
        }
      }

      label[p.idx] = res;

      for (int k = 0; k < n_roots_sz; k++) {
        int root = n_roots[k];
        if (root == res) continue;
        float diff = seed_info[root].peak_val - v;
        if (diff < tolerance) {
          parent[root] = res;
        }
      }
    }
  }

  // 6. Renumeração em matriz de saída
  IntegerMatrix out(nrow, ncol);
  VecI final_label(parent.size(), 0);
  int current_label = 0;

  for (int i = 0; i < N; i++) {
    if (label[i] > 0) {
      int r = find_root(label[i], parent);
      if (final_label[r] == 0) {
        current_label++;
        final_label[r] = current_label;
      }
      out[i] = final_label[r];
    }
  }

  return out;
}

// =============================================================================
// Distmap Euclidiano Exato e Ultra-Rápido (O(N))
// =============================================================================
// Calcula a Transformada de Distância Euclidiana Exata baseada em Felzenszwalb
//
// [[Rcpp::export]]
NumericMatrix help_dist_transform(LogicalMatrix bin) {
  int nrow = bin.nrow();
  int ncol = bin.ncol();
  int N = nrow * ncol;

  // Pré-aloca e popula buffer de entrada
  Vec8 img(N);
  for (int i = 0; i < N; i++) {
    img[i] = bin[i] ? 1u : 0u;
  }

  // Utiliza a função OpenMP otimizada já disponível no escopo (ws_edt2d_real)
  // Lembrando que ws_edt2d_real enxerga !img[i] (onde bin == 0) como o fundo (Dist = 0)
  VecF dist = ws_edt2d_real(img, nrow, ncol);

  NumericMatrix out(nrow, ncol);
  for (int i = 0; i < N; i++) {
    out[i] = dist[i];
  }

  return out;
}

