// dense_zncc.h -- shared declarations for the dense stride-1 ZNCC exemplar matching
// used by PerSAM (`feat_res = 1024`).
//
// The kernel representation lives here so that both the ONNX runner (which builds the
// kernels from SAM mask outputs) and the FFT-based correlation backend can see it.

#ifndef PLIMAN_DENSE_ZNCC_H
#define PLIMAN_DENSE_ZNCC_H

#include <vector>
#include <cstddef>

// One sample of the exemplar template, at offset (dx, dy) from the object centre.
// `t_zero` is the mean-subtracted RGB template value, `w` the radial weight.
struct ZnccKernelEntry {
  int dx;
  int dy;
  float w;
  float t_zero[3];
};

// One maximal run of consecutive dx at a fixed dy. The kernel support is a disk
// intersected with the exemplar mask, so a row is usually one run and at worst a few.
struct ZnccSpan {
  int dy;
  int dx0;
  int len;
  size_t off;   // index of the first entry of this run in the packed arrays
};

struct DenseExemplarKernel {
  std::vector<ZnccKernelEntry> entries;
  float tmpl_std;
  float sum_w;
  float inv_sum_w;
  int radius;

  // Span (structure-of-arrays) form of `entries`, built by build_zncc_spans(). The
  // interior fast path walks these instead of the entry list: within a run the image
  // samples are contiguous in memory and the weights are contiguous in `w_arr` /
  // `wt_arr`, which turns a gather over a disk into vectorizable streaming loads.
  std::vector<ZnccSpan> spans;
  std::vector<float> w_arr;    // w, one per entry, in span order
  std::vector<float> wt_arr;   // w * t_zero[c], three per entry, interleaved
};

// Packs `k.entries` into runs of consecutive dx at equal dy.
void build_zncc_spans(DenseExemplarKernel& k);

// FFT-based dense ZNCC over the image interior.
//
// Computes, for every pixel of the interior region (those at least `radius` away from
// each border, where the direct path needs no clamping),
//
//     zncc_p(y,x) = dot_p(y,x) / (patch_std_p(y,x) * tmpl_std_p)
//
// for each exemplar kernel p, and writes `v_sim(y,x) * max_p zncc_p(y,x)` into
// `dense_out` wherever that maximum is positive. Pixels whose `vit` value is below
// `vit_thresh`, and all border pixels, are left untouched for the caller to handle.
//
// Each of the six terms the ZNCC needs is a correlation of an image plane (or its
// square) with a fixed kernel, so all of them are evaluated as products of spectra:
// the cost becomes independent of the kernel radius, where the direct evaluation grows
// with its square. The image spectra are cached across calls keyed on `img_key`, so
// several exemplars on one image pay for them once.
//
// Returns false when FFTW is unavailable or the geometry is unsupported, in which case
// the caller must fall back to direct evaluation.
bool dense_zncc_fft(const float* img_c0,
                    const float* img_c1,
                    const float* img_c2,
                    int W,
                    int H,
                    unsigned long long img_key,
                    const std::vector<DenseExemplarKernel>& kernels,
                    const std::vector<float>& vit,
                    float vit_thresh,
                    std::vector<float>& dense_out);

// Frees the cached image spectra (called when the ONNX session caches are cleared).
void dense_zncc_clear_cache();

// Diagnostics (see dense_zncc.cpp): attempts and successes of the FFT backend.
extern long long g_zncc_fft_attempts;
extern long long g_zncc_fft_success;

#endif  // PLIMAN_DENSE_ZNCC_H
