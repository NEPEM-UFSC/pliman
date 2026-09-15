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
// =============================================================================
// =============================================================================
// Preenchimento de buracos Inteligente por Distância / Gargalo (Distance-Guided Specular Hole Filling)
// Preenche APENAS buracos de reflexo internos de grãos individuais.
// Cavidades formadas por 3+ grãos se tocando (gargalos com distância < min_neck_dist) NÃO são preenchidas!
// =============================================================================

static void fill_holes_fast(const Vec8& img, int nrow, int ncol, Vec8& dest, MemoBuffers& memo, double max_size = -1.0, double min_neck_dist = 2.5) {
  ext_bg_bfs(img, nrow, ncol, memo.ext, memo.q_buf);
  int N = nrow * ncol;

  if (max_size <= 0.0 && min_neck_dist <= 0.0) {
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) dest[i] = img[i] | (memo.ext[i] ^ 1u);
    return;
  }

  // 1. Transformada de distância do primeiro plano (foreground) ao fundo externo
  #pragma omp parallel for simd schedule(static)
  for (int i = 0; i < N; i++) memo.temp8[i] = img[i] ^ 1u;
  edt2d(memo.temp8, nrow, ncol, memo.dist);

  std::copy(img.begin(), img.end(), dest.begin());
  std::fill(memo.temp8.begin(), memo.temp8.end(), 0u);

  double min_dist_sq = (min_neck_dist > 0.0) ? (min_neck_dist * min_neck_dist) : 0.0;

  std::vector<int> comp;
  comp.reserve(1024);

  for (int i = 0; i < N; i++) {
    if (!img[i] && !memo.ext[i] && !memo.temp8[i]) {
      comp.clear();
      comp.push_back(i);
      memo.temp8[i] = 1u;

      float min_border_dist_sq = F_INF;

      size_t head = 0;
      while (head < comp.size()) {
        int curr = comp[head++];
        int r = curr % nrow;
        int c = curr / nrow;

        auto try_add = [&](int nr, int nc) {
          int idx = nc * nrow + nr;
          if (!img[idx] && !memo.ext[idx]) {
            if (!memo.temp8[idx]) {
              memo.temp8[idx] = 1u;
              comp.push_back(idx);
            }
          } else if (img[idx]) {
            float d_sq = memo.dist[idx];
            if (d_sq < min_border_dist_sq) {
              min_border_dist_sq = d_sq;
            }
          }
        };

        if (r > 0)        try_add(r - 1, c);
        if (r < nrow - 1) try_add(r + 1, c);
        if (c > 0)        try_add(r, c - 1);
        if (c < ncol - 1) try_add(r, c + 1);
      }

      bool fill_it = true;
      if (min_dist_sq > 0.0 && min_border_dist_sq < min_dist_sq) {
        fill_it = false; // Cavidade entre grãos que se tocam! NÃO preenche!
      }
      if (max_size > 0.0 && static_cast<double>(comp.size()) > max_size) {
        fill_it = false;
      }

      if (fill_it) {
        for (int p : comp) {
          dest[p] = 1u;
        }
      }
    }
  }
}

// =============================================================================
// Funções exportadas para R
// =============================================================================

//' Fast external erosion
//'
//' @param img A logical matrix representing the binary image.
//' @param raio Radius of the structuring element (integer >= 1).
//' @param forma Shape of the structuring element ("disc" or "square").
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

//' Standard morphological erosion
//'
//' @param img A logical matrix representing the binary image.
//' @param raio Radius of the structuring element (integer >= 1).
//' @param forma Shape of the structuring element ("disc" or "square").
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

//' Fast morphological dilation
//'
//' @param img A logical matrix representing the binary image.
//' @param raio Radius of the structuring element (integer >= 1).
//' @param forma Shape of the structuring element ("disc" or "square").
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

//' Fast internal hole filling (with optional maximum hole size and neck distance check)
//'
//' @param img A logical matrix representing the binary image.
//' @param max_size Maximum area (in pixels) of internal holes to fill. If <= 0 (default), fills all holes.
//' @param min_neck_dist Minimum foreground thickness (distance to background) surrounding a hole to allow filling. Default = 0.0 (fill all internal holes).
// [[Rcpp::export]]
LogicalMatrix fill_holes_cpp(LogicalMatrix img, double max_size = -1.0, double min_neck_dist = 0.0) {
  int nrow = img.nrow(), ncol = img.ncol();
  int N = nrow * ncol;
  Vec8 v = lm_to_vec(img);
  Vec8 dest(N);
  MemoBuffers memo(N);
  fill_holes_fast(v, nrow, ncol, dest, memo, max_size, min_neck_dist);
  return vec_to_lm(dest, nrow, ncol);
}

//' Fast conservative erosion pipeline
//'
//' @param img A logical matrix representing the binary image.
//' @param raio_erosao Radius for erosion step (integer >= 1).
//' @param raio_dilatacao Radius for dilation step (integer >= 1).
//' @param forma Shape of the structuring element ("disc" or "square").
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


// Filtro de mediana binária trabalhando diretamente em Vec8
static void median_filter_vec8(const Vec8& src, int nrow, int ncol, int radius, Vec8& dest) {
  int N = nrow * ncol;
  if (radius < 1) {
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; ++i) dest[i] = src[i];
    return;
  }
  
  #pragma omp parallel for schedule(static)
  for (int c = 0; c < ncol; ++c) {
    for (int r = 0; r < nrow; ++r) {
      int win_ones = 0;
      int win_total = 0;
      int r_start = std::max(0, r - radius);
      int r_end = std::min(nrow - 1, r + radius);
      int c_start = std::max(0, c - radius);
      int c_end = std::min(ncol - 1, c + radius);
      
      for (int wc = c_start; wc <= c_end; ++wc) {
        const uint8_t* col_ptr = src.data() + (size_t)wc * nrow;
        for (int wr = r_start; wr <= r_end; ++wr) {
          if (col_ptr[wr]) win_ones++;
          win_total++;
        }
      }
      dest[(size_t)c * nrow + r] = (win_ones >= (win_total + 1) / 2) ? 1u : 0u;
    }
  }
}

double help_otsu(SEXP img_sexp);

//' Grayscale Reflection Removal via Masked Median Inpainting
//'
//' Removes specular reflections from a grayscale index matrix before thresholding.
//' Algorithm:
//'   1. Auto-detects which class is background by sampling border pixels after Otsu.
//'   2. BFS from image borders marks all background pixels reachable from the exterior.
//'   3. Internal holes = background-class pixels NOT reached by BFS (enclosed by grain).
//'   4. Each internal hole pixel is replaced by the MEDIAN of its foreground (grain)
//'      neighbours within `radius`. This is orientation-agnostic and never touches
//'      the grain boundary or the external background.
//'
//' @param img Grayscale index matrix.
//' @param radius Neighbourhood radius for median sampling (default = 5).
//' @param has_white_bg Unused — background class is detected automatically.
//' @return Corrected grayscale matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix remove_reflection_cpp(NumericMatrix img, int radius = 5, bool has_white_bg = true) {
  int nrow = img.nrow(), ncol = img.ncol(), N = nrow * ncol;
  if (N <= 0 || radius < 1) return img;

  const double* p_in = REAL(img);

  // -----------------------------------------------------------------------
  // Step 1: Otsu threshold
  // -----------------------------------------------------------------------
  double thresh = help_otsu(img);
  {
    double mn = p_in[0], mx = p_in[0];
    for (int i = 1; i < N; i++) {
      if (p_in[i] < mn) mn = p_in[i];
      if (p_in[i] > mx) mx = p_in[i];
    }
    if (thresh <= mn + 1e-9 || thresh >= mx - 1e-9) thresh = (mn + mx) * 0.5;
  }

  // -----------------------------------------------------------------------
  // Step 2: Determine background class from border-pixel majority vote
  // -----------------------------------------------------------------------
  int n_above = 0, n_total = 0;
  for (int r = 0; r < nrow; r++) {
    if (p_in[0 * nrow + r]         >= thresh) n_above++;
    if (p_in[(ncol-1) * nrow + r]  >= thresh) n_above++;
    n_total += 2;
  }
  for (int c = 1; c < ncol - 1; c++) {
    if (p_in[c * nrow + 0]        >= thresh) n_above++;
    if (p_in[c * nrow + nrow - 1] >= thresh) n_above++;
    n_total += 2;
  }
  // bg_is_above: TRUE = background pixels have HIGH index values
  bool bg_is_above = (n_above * 2 >= n_total);  // majority vote

  // -----------------------------------------------------------------------
  // Step 3: BFS from image borders — mark all externally-connected bg pixels
  // -----------------------------------------------------------------------
  std::vector<uint8_t> ext_bg(N, 0u);
  std::vector<int> q;
  q.reserve(N / 4);

  auto try_seed = [&](int r, int c) {
    int idx = c * nrow + r;
    bool is_bg = bg_is_above ? (p_in[idx] >= thresh) : (p_in[idx] < thresh);
    if (is_bg && !ext_bg[idx]) { ext_bg[idx] = 1u; q.push_back(idx); }
  };
  for (int r = 0; r < nrow; r++) { try_seed(r, 0); try_seed(r, ncol - 1); }
  for (int c = 1; c < ncol - 1; c++) { try_seed(0, c); try_seed(nrow - 1, c); }

  for (size_t qi = 0; qi < q.size(); qi++) {
    int idx = q[qi];
    int r = idx % nrow, c = idx / nrow;
    auto try_push = [&](int nr, int nc) {
      int ni = nc * nrow + nr;
      bool is_bg = bg_is_above ? (p_in[ni] >= thresh) : (p_in[ni] < thresh);
      if (is_bg && !ext_bg[ni]) { ext_bg[ni] = 1u; q.push_back(ni); }
    };
    if (r > 0)        try_push(r - 1, c);
    if (r < nrow - 1) try_push(r + 1, c);
    if (c > 0)        try_push(r, c - 1);
    if (c < ncol - 1) try_push(r, c + 1);
  }

  // -----------------------------------------------------------------------
  // Step 4: For each internal hole pixel, replace with median of fg neighbours
  // -----------------------------------------------------------------------
  NumericMatrix out = clone(img);
  double* p_out = REAL(out);

  // (serial — nth_element is not trivially parallelisable, and holes are rare)
  for (int c = 0; c < ncol; c++) {
    for (int r = 0; r < nrow; r++) {
      int idx = c * nrow + r;
      bool is_bg = bg_is_above ? (p_in[idx] >= thresh) : (p_in[idx] < thresh);
      if (!is_bg || ext_bg[idx]) continue;  // skip grain pixels and external bg

      // Internal hole: collect foreground (grain) neighbours
      int r_lo = std::max(0, r - radius), r_hi = std::min(nrow - 1, r + radius);
      int c_lo = std::max(0, c - radius), c_hi = std::min(ncol - 1, c + radius);

      std::vector<double> fg_vals;
      fg_vals.reserve((r_hi - r_lo + 1) * (c_hi - c_lo + 1));
      for (int wc = c_lo; wc <= c_hi; wc++) {
        for (int wr = r_lo; wr <= r_hi; wr++) {
          int ni = wc * nrow + wr;
          bool n_is_fg = bg_is_above ? (p_in[ni] < thresh) : (p_in[ni] >= thresh);
          if (n_is_fg) fg_vals.push_back(p_in[ni]);
        }
      }

      if (!fg_vals.empty()) {
        size_t mid = fg_vals.size() / 2;
        std::nth_element(fg_vals.begin(), fg_vals.begin() + mid, fg_vals.end());
        p_out[idx] = fg_vals[mid];
      }
    }
  }

  return out;
}

// Pipeline morfológico sequencial unificado para binarização rápida
// [[Rcpp::export]]
LogicalMatrix help_binary_filters_cpp(SEXP img_sexp,
                                      int erode = 0,
                                      int dilate = 0,
                                      int opening = 0,
                                      int closing = 0,
                                      int filter = 0,
                                      bool fill_hull = false,
                                      CharacterVector filter_order = CharacterVector::create("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                                      double max_size = -1.0,
                                      double min_neck_dist = 0.0) {
  int nrow = Rf_nrows(img_sexp);
  int ncol = Rf_ncols(img_sexp);
  if (nrow == 0 || ncol == 0) {
    SEXP dims = Rf_getAttrib(img_sexp, R_DimSymbol);
    if (!Rf_isNull(dims) && Rf_length(dims) >= 2) {
      nrow = INTEGER(dims)[0];
      ncol = INTEGER(dims)[1];
    }
  }
  int N = nrow * ncol;
  if (N <= 0) return LogicalMatrix(0, 0);

  Vec8 v(N);
  if (TYPEOF(img_sexp) == LGLSXP || TYPEOF(img_sexp) == INTSXP) {
    const int* ptr = INTEGER(img_sexp);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) v[i] = ptr[i] ? 1u : 0u;
  } else if (TYPEOF(img_sexp) == RAWSXP) {
    const uint8_t* ptr = RAW(img_sexp);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) v[i] = ptr[i] ? 1u : 0u;
  } else if (TYPEOF(img_sexp) == REALSXP) {
    const double* ptr = REAL(img_sexp);
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < N; i++) v[i] = (ptr[i] != 0.0 && !std::isnan(ptr[i])) ? 1u : 0u;
  }

  if (erode <= 0 && dilate <= 0 && opening <= 0 && closing <= 0 && filter <= 1 && !fill_hull) {
    return vec_to_lm(v, nrow, ncol);
  }

  Vec8 temp(N);
  MemoBuffers memo(N);

  Vec8* p_src = &v;
  Vec8* p_dest = &temp;

  for (int idx = 0; idx < filter_order.size(); ++idx) {
    std::string op = Rcpp::as<std::string>(filter_order[idx]);
    if (op == "erode" && erode > 0) {
      ero_std_fast(*p_src, nrow, ncol, erode, true, *p_dest, memo);
      std::swap(p_src, p_dest);
    } else if (op == "dilate" && dilate > 0) {
      dil_fast(*p_src, nrow, ncol, dilate, true, *p_dest, memo.dist, memo.hdil);
      std::swap(p_src, p_dest);
    } else if (op == "opening" && opening > 0) {
      ero_std_fast(*p_src, nrow, ncol, opening, true, *p_dest, memo);
      std::swap(p_src, p_dest);
      dil_fast(*p_src, nrow, ncol, opening, true, *p_dest, memo.dist, memo.hdil);
      std::swap(p_src, p_dest);
    } else if (op == "closing" && closing > 0) {
      dil_fast(*p_src, nrow, ncol, closing, true, *p_dest, memo.dist, memo.hdil);
      std::swap(p_src, p_dest);
      ero_std_fast(*p_src, nrow, ncol, closing, true, *p_dest, memo);
      std::swap(p_src, p_dest);
    } else if (op == "filter" && filter > 1) {
      median_filter_vec8(*p_src, nrow, ncol, filter, *p_dest);
      std::swap(p_src, p_dest);
    } else if (op == "fill_hull" && fill_hull) {
      fill_holes_fast(*p_src, nrow, ncol, *p_dest, memo, max_size, min_neck_dist);
      std::swap(p_src, p_dest);
    }
  }

  return vec_to_lm(*p_src, nrow, ncol);
}


