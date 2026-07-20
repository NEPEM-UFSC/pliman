// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// Tipos internos: uint8 (1 byte/pixel)
using Vec8 = std::vector<uint8_t>;
using VecF = std::vector<float>;
using VecI = std::vector<int>;

static constexpr float F_INF = 1e15f;

// Struct para gerenciar e reutilizar os buffers temporários de tamanho N.
// Isso evita múltiplos 'malloc' / 'free' lentos e fragmentação de memória.
struct MemoBuffers {
  Vec8 ext;
  Vec8 temp8;
  Vec8 hdil;
  VecF dist;
  VecI q_buf;

  MemoBuffers(int N) {
    ext.resize(N, 0u);
    temp8.resize(N, 0u);
    hdil.resize(N, 0u);
    dist.resize(N, 0.0f);
    q_buf.reserve(N / 4);
  }
};

// =============================================================================
// Conversões LogicalMatrix <-> Vec8 (flat, column-major para respeitar o R)
// =============================================================================

static inline Vec8 lm_to_vec(const LogicalMatrix& m) {
  int N = m.nrow() * m.ncol();
  Vec8 v(N);
  for (int i = 0; i < N; i++) v[i] = m[i] ? 1u : 0u;
  return v;
}

static inline LogicalMatrix vec_to_lm(const Vec8& v, int nrow, int ncol) {
  LogicalMatrix m(nrow, ncol);
  int N = nrow * ncol;
  for (int i = 0; i < N; i++) m[i] = (v[i] != 0u);
  return m;
}

// =============================================================================
// BFS: encontra fundo externo (pixels FALSE conectados às bordas da imagem)
// Indexação nativa em column-major: idx = c * nrow + r
// =============================================================================

static void ext_bg_bfs(const Vec8& img, int nrow, int ncol, Vec8& ext, VecI& q) {
  int N = nrow * ncol;
  std::fill(ext.begin(), ext.end(), 0u);
  q.clear();

  // Helper para semear na fila
  auto push = [&](int r, int c) {
    int idx = c * nrow + r;
    if (!img[idx] && !ext[idx]) {
      ext[idx] = 1u;
      q.push_back(idx);
    }
  };

  // Semeia todas as bordas externas da imagem
  for (int r = 0; r < nrow; r++) {
    push(r, 0);
    push(r, ncol - 1);
  }
  for (int c = 1; c < ncol - 1; c++) {
    push(0, c);
    push(nrow - 1, c);
  }

  // BFS usando vetor como fila para ótima localidade de cache
  for (int qi = 0; qi < (int)q.size(); qi++) {
    int idx = q[qi];
    int r = idx % nrow;
    int c = idx / nrow;

    auto try_push = [&](int nr, int nc) {
      int n_idx = nc * nrow + nr;
      if (!img[n_idx] && !ext[n_idx]) {
        ext[n_idx] = 1u;
        q.push_back(n_idx);
      }
    };

    if (r > 0)        try_push(r - 1, c);
    if (r < nrow - 1) try_push(r + 1, c);
    if (c > 0)        try_push(r, c - 1);
    if (c < ncol - 1) try_push(r, c + 1);
  }
}

// =============================================================================
// 1D Felzenszwalb Distance Transform (in-place)
// =============================================================================

static void dt1d(float* f, float* tmp, int* v, float* z, int n) {
  int k = 0;
  v[0] = 0; z[0] = -F_INF; z[1] = F_INF;

  for (int q = 1; q < n; q++) {
    float fq = f[q];
    float s;
    do {
      int vk = v[k];
      s = ((fq + (float)(q*q)) - (f[vk] + (float)(vk*vk)))
          / (float)(2 * (q - vk));
      if (s > z[k]) break;
      --k;
    } while (true);
    ++k; v[k] = q; z[k] = s; z[k+1] = F_INF;
  }

  k = 0;
  for (int q = 0; q < n; q++) {
    while (z[k+1] < (float)q) ++k;
    int vk = v[k];
    tmp[q] = (float)(q - vk) * (float)(q - vk) + f[vk];
  }
  for (int q = 0; q < n; q++) f[q] = tmp[q];
}

// =============================================================================
// EDT 2D (distância euclidiana ao quadrado em column-major)
// =============================================================================

static void edt2d(const Vec8& src, int nrow, int ncol, VecF& D) {
  int N = nrow * ncol;
  std::fill(D.begin(), D.end(), F_INF);
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) if (src[i]) D[i] = 0.0f;

  #pragma omp parallel
  {
    // Passo 1: DT 1D ao longo das colunas (blocos contíguos na memória!)
    VecF tmp_c(nrow);
    VecI v_c(nrow + 1);
    VecF z_c(nrow + 2);
    #pragma omp for schedule(static)
    for (int c = 0; c < ncol; c++) {
      float* col_ptr = D.data() + c * nrow;
      dt1d(col_ptr, tmp_c.data(), v_c.data(), z_c.data(), nrow);
    }

    // Passo 2: DT 1D ao longo das linhas (elementos esparsos em column-major)
    VecF row(ncol), tmp_r(ncol);
    VecI v_r(ncol + 1);
    VecF z_r(ncol + 2);
    #pragma omp for schedule(static)
    for (int r = 0; r < nrow; r++) {
      for (int c = 0; c < ncol; c++) row[c] = D[c * nrow + r];
      dt1d(row.data(), tmp_r.data(), v_r.data(), z_r.data(), ncol);
      for (int c = 0; c < ncol; c++) D[c * nrow + r] = row[c];
    }
  }
}

// =============================================================================
// Dilatação por Prefix Sums para elemento estruturante quadrado
// =============================================================================

static void dil_square_pfx(const Vec8& src, int nrow, int ncol, int raio, Vec8& dest, Vec8& hdil) {
  // Passo 1: dilatação vertical (ao longo de cada coluna contígua)
  #pragma omp parallel for schedule(static)
  for (int c = 0; c < ncol; c++) {
    const uint8_t* col_in = src.data() + c * nrow;
    uint8_t* col_vd = hdil.data() + c * nrow;
    VecI pfx(nrow + 1, 0);
    for (int r = 0; r < nrow; r++) pfx[r+1] = pfx[r] + col_in[r];
    for (int r = 0; r < nrow; r++) {
      int lo = std::max(0, r - raio);
      int hi = std::min(nrow, r + raio + 1);
      col_vd[r] = (pfx[hi] - pfx[lo]) > 0 ? 1u : 0u;
    }
  }

  // Passo 2: dilatação horizontal (ao longo de cada linha)
  #pragma omp parallel for schedule(static)
  for (int r = 0; r < nrow; r++) {
    VecI pfx(ncol + 1, 0);
    for (int c = 0; c < ncol; c++) pfx[c+1] = pfx[c] + hdil[c * nrow + r];
    for (int c = 0; c < ncol; c++) {
      int lo = std::max(0, c - raio);
      int hi = std::min(ncol, c + raio + 1);
      dest[c * nrow + r] = (pfx[hi] - pfx[lo]) > 0 ? 1u : 0u;
    }
  }
}

// =============================================================================
// Dilatação principal
// =============================================================================

static void dil_fast(const Vec8& src, int nrow, int ncol, int raio, bool disc, Vec8& dest, VecF& dist_buf, Vec8& hdil_buf) {
  int N = nrow * ncol;
  if (disc) {
    edt2d(src, nrow, ncol, dist_buf);
    float r2 = (float)raio * (float)raio;
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) dest[i] = (dist_buf[i] <= r2) ? 1u : 0u;
  } else {
    dil_square_pfx(src, nrow, ncol, raio, dest, hdil_buf);
  }
}

// =============================================================================
// Erosão PADRÃO: NOT dil(NOT img)
// =============================================================================

static void ero_std_fast(const Vec8& img, int nrow, int ncol, int raio, bool disc, Vec8& dest, MemoBuffers& memo) {
  int N = nrow * ncol;
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) memo.temp8[i] = img[i] ^ 1u;

  dil_fast(memo.temp8, nrow, ncol, raio, disc, memo.ext, memo.dist, memo.hdil);

  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) dest[i] = memo.ext[i] ^ 1u;
}

// =============================================================================
// Erosão EXTERNA seletiva: img AND NOT dil(ext_bg)
// =============================================================================

static void ero_ext_fast(const Vec8& img, int nrow, int ncol, int raio, bool disc, Vec8& dest, MemoBuffers& memo) {
  ext_bg_bfs(img, nrow, ncol, memo.ext, memo.q_buf);

  // Adiciona a borda externa virtualmente
  for (int r = 0; r < nrow; r++) { memo.ext[r] |= 1u; memo.ext[(ncol-1)*nrow + r] |= 1u; }
  for (int c = 0; c < ncol; c++) { memo.ext[c * nrow] |= 1u; memo.ext[c * nrow + nrow - 1] |= 1u; }

  dil_fast(memo.ext, nrow, ncol, raio, disc, memo.temp8, memo.dist, memo.hdil);

  int N = nrow * ncol;
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) {
    dest[i] = img[i] & (~memo.temp8[i] & 1u);
  }
}

// =============================================================================
// Preenchimento de buracos: img OR NOT ext_bg
// =============================================================================

static void fill_holes_fast(const Vec8& img, int nrow, int ncol, Vec8& dest, MemoBuffers& memo) {
  ext_bg_bfs(img, nrow, ncol, memo.ext, memo.q_buf);
  int N = nrow * ncol;
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) dest[i] = img[i] | (memo.ext[i] ^ 1u);
}

// =============================================================================
// Funções exportadas para R
// =============================================================================

//' Erosão externa rápida (preserva borda interna de objetos fechados)
// [[Rcpp::export]]
LogicalMatrix erode_external_cpp(LogicalMatrix img,
                                 int raio = 3,
                                 std::string forma = "disc") {
  if (raio < 1) stop("raio deve ser >= 1");
  if (forma != "disc" && forma != "square") stop("forma deve ser 'disc' ou 'square'");
  bool disc = (forma == "disc");
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  Vec8 dest(N);
  MemoBuffers memo(N);
  ero_ext_fast(v, nrow, ncol, raio, disc, dest, memo);
  return vec_to_lm(dest, nrow, ncol);
}

//' Erosão morfológica padrão (todas as bordas)
// [[Rcpp::export]]
LogicalMatrix erode_cpp(LogicalMatrix img,
                        int raio = 3,
                        std::string forma = "disc") {
  if (raio < 1) stop("raio deve ser >= 1");
  if (forma != "disc" && forma != "square") stop("forma deve ser 'disc' ou 'square'");
  bool disc = (forma == "disc");
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  Vec8 dest(N);
  MemoBuffers memo(N);
  ero_std_fast(v, nrow, ncol, raio, disc, dest, memo);
  return vec_to_lm(dest, nrow, ncol);
}

//' Dilatação morfológica rápida
// [[Rcpp::export]]
LogicalMatrix dilate_cpp(LogicalMatrix img,
                         int raio = 3,
                         std::string forma = "disc") {
  if (raio < 1) stop("raio deve ser >= 1");
  if (forma != "disc" && forma != "square") stop("forma deve ser 'disc' ou 'square'");
  bool disc = (forma == "disc");
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  Vec8 dest(N);
  MemoBuffers memo(N);
  dil_fast(v, nrow, ncol, raio, disc, dest, memo.dist, memo.hdil);
  return vec_to_lm(dest, nrow, ncol);
}

//' Preenchimento de buracos internos rápido
// [[Rcpp::export]]
LogicalMatrix fill_holes_cpp(LogicalMatrix img) {
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  Vec8 dest(N);
  MemoBuffers memo(N);
  fill_holes_fast(v, nrow, ncol, dest, memo);
  return vec_to_lm(dest, nrow, ncol);
}

//' Erosão conservativa rápida (pipeline completo)
// [[Rcpp::export]]
LogicalMatrix erosao_conservativa(LogicalMatrix img,
                                  int raio_erosao    = 3,
                                  int raio_dilatacao = 3,
                                  std::string forma  = "disc") {
  if (raio_erosao < 1)    stop("raio_erosao deve ser >= 1");
  if (raio_dilatacao < 1) stop("raio_dilatacao deve ser >= 1");
  if (forma != "disc" && forma != "square")
    stop("forma deve ser 'disc' ou 'square'");

  bool disc = (forma == "disc");
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  MemoBuffers memo(N);

  Vec8 erodida(N);
  ero_ext_fast(v, nrow, ncol, raio_erosao, disc, erodida, memo);

  Vec8 preenchida(N);
  fill_holes_fast(erodida, nrow, ncol, preenchida, memo);

  Vec8 dilatada(N);
  dil_fast(preenchida, nrow, ncol, raio_dilatacao, disc, dilatada, memo.dist, memo.hdil);

  return vec_to_lm(dilatada, nrow, ncol);
}

