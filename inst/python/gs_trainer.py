#!/usr/bin/env python3
"""
GPU-Accelerated 3D Gaussian Splatting Training Engine for pliman
=================================================================

Faithful implementation of Kerbl et al. (2023) "3D Gaussian Splatting for
Real-Time Radiance Field Rendering":

* Differentiable EWA projection of anisotropic 3D Gaussians (scales +
  rotation quaternions) and spherical-harmonics colour evaluation are written
  in PyTorch, so their gradients come from autograd.
* Tile-based front-to-back alpha compositing and its exact analytical backward
  pass are custom CUDA kernels compiled at runtime with CuPy / NVRTC (no nvcc
  or MSVC toolchain needed).
* Loss: (1 - lambda) * L1 + lambda * (1 - SSIM).
* Adam with per-parameter learning rates and exponential decay of the position
  learning rate (scaled by the scene extent).
* Adaptive density control: clone / split driven by accumulated screen-space
  positional gradients, pruning of transparent / oversized Gaussians and
  periodic opacity reset.
* Progressive spherical-harmonics degree (one band every 1000 iterations).

Binary interchange format with R (all little-endian):
    int32 magic (0x32565347 = "GSV2"), int32 N, int32 K (SH coefficient
    columns), float32 positions[N,3], scales_log[N,3], quaternions[N,4],
    opacities_logit[N], sh[N,K], int32 semantic_ids[N], float32 scores[N]
The legacy format (no magic, K = 3) is also accepted on input.
"""

import sys
import os
import math
import time
import json
import argparse

import numpy as np

# Force UTF-8 stdout/stderr on Windows to prevent cp1252 charmap errors
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# Automatically resolve CUDA runtime headers / NVRTC shipped as pip wheels
try:
    import nvidia.cuda_runtime
    cuda_rt_dir = nvidia.cuda_runtime.__path__[0]
    if "CUDA_PATH" not in os.environ or not os.path.exists(os.environ["CUDA_PATH"]):
        os.environ["CUDA_PATH"] = cuda_rt_dir
except Exception:
    pass

import warnings
warnings.filterwarnings("ignore", message="CUDA path could not be detected")

import torch
import torch.nn.functional as F
import cupy
from PIL import Image, ImageOps

MAGIC_V2 = 0x32565347
TILE = 16
BLOCK = TILE * TILE

# -----------------------------------------------------------------------------
# CUDA kernels: tile rasterizer forward + exact backward
# -----------------------------------------------------------------------------
CUDA_SRC = r"""
#define TILE 16
#define BLOCK 256

extern "C" __global__ void render_fwd(
    const int*   __restrict__ ranges,      // [num_tiles, 2]
    const int*   __restrict__ point_list,  // [num_instances]
    const float* __restrict__ means2d,     // [P, 2]
    const float* __restrict__ conics,      // [P, 3]  (a, b, c)
    const float* __restrict__ opac,        // [P]
    const float* __restrict__ colors,      // [P, 3]
    const float* __restrict__ bg,          // [3]
    const int W, const int H, const int tiles_x,
    float* __restrict__ out_img,           // [H, W, 3]
    float* __restrict__ final_T,           // [H, W]
    int*   __restrict__ n_contrib)         // [H, W]
{
    const int tile = blockIdx.y * tiles_x + blockIdx.x;
    const int px = blockIdx.x * TILE + threadIdx.x;
    const int py = blockIdx.y * TILE + threadIdx.y;
    const int tid = threadIdx.y * TILE + threadIdx.x;
    const bool inside = (px < W) && (py < H);
    bool done = !inside;
    const float pfx = (float)px + 0.5f;
    const float pfy = (float)py + 0.5f;

    const int start = ranges[2 * tile];
    const int end   = ranges[2 * tile + 1];
    const int rounds = (end - start + BLOCK - 1) / BLOCK;
    int todo = end - start;

    __shared__ int    s_id[BLOCK];
    __shared__ float2 s_xy[BLOCK];
    __shared__ float4 s_co[BLOCK];

    float T = 1.0f;
    float C0 = 0.0f, C1 = 0.0f, C2 = 0.0f;
    int contributor = 0;
    int last_contributor = 0;

    for (int r = 0; r < rounds; r++, todo -= BLOCK) {
        int num_done = __syncthreads_count(done);
        if (num_done == BLOCK) break;
        int idx = start + r * BLOCK + tid;
        if (idx < end) {
            int g = point_list[idx];
            s_id[tid] = g;
            s_xy[tid] = make_float2(means2d[2 * g], means2d[2 * g + 1]);
            s_co[tid] = make_float4(conics[3 * g], conics[3 * g + 1], conics[3 * g + 2], opac[g]);
        }
        __syncthreads();
        int n = min(BLOCK, todo);
        for (int j = 0; !done && j < n; j++) {
            contributor++;
            float2 xy = s_xy[j];
            float4 co = s_co[j];
            float dx = xy.x - pfx;
            float dy = xy.y - pfy;
            float power = -0.5f * (co.x * dx * dx + co.z * dy * dy) - co.y * dx * dy;
            if (power > 0.0f) continue;
            float alpha = fminf(0.99f, co.w * __expf(power));
            if (alpha < (1.0f / 255.0f)) continue;
            float test_T = T * (1.0f - alpha);
            if (test_T < 0.0001f) { done = true; continue; }
            int g = s_id[j];
            float w = alpha * T;
            C0 += colors[3 * g + 0] * w;
            C1 += colors[3 * g + 1] * w;
            C2 += colors[3 * g + 2] * w;
            T = test_T;
            last_contributor = contributor;
        }
    }

    if (inside) {
        int pix = py * W + px;
        final_T[pix] = T;
        n_contrib[pix] = last_contributor;
        out_img[3 * pix + 0] = C0 + T * bg[0];
        out_img[3 * pix + 1] = C1 + T * bg[1];
        out_img[3 * pix + 2] = C2 + T * bg[2];
    }
}

extern "C" __global__ void render_bwd(
    const int*   __restrict__ ranges,
    const int*   __restrict__ point_list,
    const float* __restrict__ means2d,
    const float* __restrict__ conics,
    const float* __restrict__ opac,
    const float* __restrict__ colors,
    const float* __restrict__ bg,
    const int W, const int H, const int tiles_x,
    const float* __restrict__ final_T,
    const int*   __restrict__ n_contrib,
    const float* __restrict__ dL_dout,     // [H, W, 3]
    float* __restrict__ dL_dmeans2d,       // [P, 2]
    float* __restrict__ dL_dconics,        // [P, 3]
    float* __restrict__ dL_dopac,          // [P]
    float* __restrict__ dL_dcolors)        // [P, 3]
{
    const int tile = blockIdx.y * tiles_x + blockIdx.x;
    const int px = blockIdx.x * TILE + threadIdx.x;
    const int py = blockIdx.y * TILE + threadIdx.y;
    const int tid = threadIdx.y * TILE + threadIdx.x;
    const bool inside = (px < W) && (py < H);
    const float pfx = (float)px + 0.5f;
    const float pfy = (float)py + 0.5f;

    const int start = ranges[2 * tile];
    const int end   = ranges[2 * tile + 1];
    const int rounds = (end - start + BLOCK - 1) / BLOCK;
    int todo = end - start;

    __shared__ int    s_id[BLOCK];
    __shared__ float2 s_xy[BLOCK];
    __shared__ float4 s_co[BLOCK];
    __shared__ float  s_col[BLOCK * 3];

    const int pix = py * W + px;
    const float T_final = inside ? final_T[pix] : 0.0f;
    float T = T_final;
    int contributor = end - start;
    const int last_contributor = inside ? n_contrib[pix] : 0;

    float dL_dpix0 = 0.0f, dL_dpix1 = 0.0f, dL_dpix2 = 0.0f;
    if (inside) {
        dL_dpix0 = dL_dout[3 * pix + 0];
        dL_dpix1 = dL_dout[3 * pix + 1];
        dL_dpix2 = dL_dout[3 * pix + 2];
    }
    const float bg_dot = bg[0] * dL_dpix0 + bg[1] * dL_dpix1 + bg[2] * dL_dpix2;

    float acc0 = 0.0f, acc1 = 0.0f, acc2 = 0.0f;
    float last_alpha = 0.0f;
    float lc0 = 0.0f, lc1 = 0.0f, lc2 = 0.0f;

    for (int r = 0; r < rounds; r++, todo -= BLOCK) {
        __syncthreads();
        int progress = r * BLOCK + tid;
        if (start + progress < end) {
            int g = point_list[end - progress - 1];
            s_id[tid] = g;
            s_xy[tid] = make_float2(means2d[2 * g], means2d[2 * g + 1]);
            s_co[tid] = make_float4(conics[3 * g], conics[3 * g + 1], conics[3 * g + 2], opac[g]);
            s_col[3 * tid + 0] = colors[3 * g + 0];
            s_col[3 * tid + 1] = colors[3 * g + 1];
            s_col[3 * tid + 2] = colors[3 * g + 2];
        }
        __syncthreads();
        int n = min(BLOCK, todo);
        // Every thread of the block walks the same Gaussian j at the same time, so the
        // per-pixel gradient contributions are first reduced inside each warp with
        // shuffles and committed with one atomicAdd per warp (up to 32x fewer atomics).
        for (int j = 0; j < n; j++) {
            contributor--;
            float v0 = 0.f, v1 = 0.f, v2 = 0.f, v3 = 0.f, v4 = 0.f, v5 = 0.f, v6 = 0.f, v7 = 0.f, v8 = 0.f;
            bool contrib = false;
            if (inside && contributor < last_contributor) {
                float2 xy = s_xy[j];
                float4 co = s_co[j];
                float dx = xy.x - pfx;
                float dy = xy.y - pfy;
                float power = -0.5f * (co.x * dx * dx + co.z * dy * dy) - co.y * dx * dy;
                if (power <= 0.0f) {
                    float G = __expf(power);
                    float alpha = fminf(0.99f, co.w * G);
                    if (alpha >= (1.0f / 255.0f)) {
                        contrib = true;
                        T = T / (1.0f - alpha);
                        const float dchannel_dcolor = alpha * T;
                        float c0 = s_col[3 * j + 0], c1 = s_col[3 * j + 1], c2 = s_col[3 * j + 2];

                        acc0 = last_alpha * lc0 + (1.0f - last_alpha) * acc0;
                        acc1 = last_alpha * lc1 + (1.0f - last_alpha) * acc1;
                        acc2 = last_alpha * lc2 + (1.0f - last_alpha) * acc2;
                        lc0 = c0; lc1 = c1; lc2 = c2;

                        float dL_dalpha = (c0 - acc0) * dL_dpix0 + (c1 - acc1) * dL_dpix1 + (c2 - acc2) * dL_dpix2;
                        v6 = dchannel_dcolor * dL_dpix0;
                        v7 = dchannel_dcolor * dL_dpix1;
                        v8 = dchannel_dcolor * dL_dpix2;

                        dL_dalpha *= T;
                        last_alpha = alpha;
                        dL_dalpha += (-T_final / (1.0f - alpha)) * bg_dot;

                        const float dL_dG = co.w * dL_dalpha;
                        const float gdx = G * dx;
                        const float gdy = G * dy;
                        v0 = dL_dG * (-gdx * co.x - gdy * co.y);
                        v1 = dL_dG * (-gdy * co.z - gdx * co.y);
                        v2 = -0.5f * gdx * dx * dL_dG;
                        v3 = -gdx * dy * dL_dG;
                        v4 = -0.5f * gdy * dy * dL_dG;
                        v5 = G * dL_dalpha;
                    }
                }
            }
            if (__any_sync(0xffffffffu, contrib)) {
                #pragma unroll
                for (int off = 16; off > 0; off >>= 1) {
                    v0 += __shfl_down_sync(0xffffffffu, v0, off);
                    v1 += __shfl_down_sync(0xffffffffu, v1, off);
                    v2 += __shfl_down_sync(0xffffffffu, v2, off);
                    v3 += __shfl_down_sync(0xffffffffu, v3, off);
                    v4 += __shfl_down_sync(0xffffffffu, v4, off);
                    v5 += __shfl_down_sync(0xffffffffu, v5, off);
                    v6 += __shfl_down_sync(0xffffffffu, v6, off);
                    v7 += __shfl_down_sync(0xffffffffu, v7, off);
                    v8 += __shfl_down_sync(0xffffffffu, v8, off);
                }
                if ((tid & 31) == 0) {
                    const int g = s_id[j];
                    atomicAdd(&dL_dmeans2d[2 * g + 0], v0);
                    atomicAdd(&dL_dmeans2d[2 * g + 1], v1);
                    atomicAdd(&dL_dconics[3 * g + 0], v2);
                    atomicAdd(&dL_dconics[3 * g + 1], v3);
                    atomicAdd(&dL_dconics[3 * g + 2], v4);
                    atomicAdd(&dL_dopac[g], v5);
                    atomicAdd(&dL_dcolors[3 * g + 0], v6);
                    atomicAdd(&dL_dcolors[3 * g + 1], v7);
                    atomicAdd(&dL_dcolors[3 * g + 2], v8);
                }
            }
        }
    }
}
"""

_MODULE = None


def get_kernels():
    global _MODULE
    if _MODULE is None:
        _MODULE = cupy.RawModule(code=CUDA_SRC, options=("--use_fast_math",))
    return _MODULE.get_function("render_fwd"), _MODULE.get_function("render_bwd")


def _cp(t):
    """Zero-copy view of a contiguous torch CUDA tensor as a CuPy array."""
    return cupy.asarray(t.detach().contiguous())


class RasterizeGaussians(torch.autograd.Function):
    @staticmethod
    def forward(ctx, means2d, conics, opac, colors, point_list, ranges, W, H, bg):
        fwd, _ = get_kernels()
        tiles_x = (W + TILE - 1) // TILE
        tiles_y = (H + TILE - 1) // TILE
        dev = means2d.device
        means2d = means2d.contiguous()
        conics = conics.contiguous()
        opac = opac.contiguous()
        colors = colors.contiguous()
        out = torch.zeros((H, W, 3), dtype=torch.float32, device=dev)
        final_T = torch.zeros((H, W), dtype=torch.float32, device=dev)
        n_contrib = torch.zeros((H, W), dtype=torch.int32, device=dev)
        fwd((tiles_x, tiles_y), (TILE, TILE),
            (_cp(ranges), _cp(point_list), _cp(means2d), _cp(conics), _cp(opac), _cp(colors),
             _cp(bg), np.int32(W), np.int32(H), np.int32(tiles_x),
             _cp(out), _cp(final_T), _cp(n_contrib)))
        ctx.save_for_backward(means2d, conics, opac, colors, point_list, ranges, bg, final_T, n_contrib)
        ctx.W, ctx.H = W, H
        return out

    @staticmethod
    def backward(ctx, grad_out):
        _, bwd = get_kernels()
        means2d, conics, opac, colors, point_list, ranges, bg, final_T, n_contrib = ctx.saved_tensors
        W, H = ctx.W, ctx.H
        tiles_x = (W + TILE - 1) // TILE
        tiles_y = (H + TILE - 1) // TILE
        g_means = torch.zeros_like(means2d)
        g_conics = torch.zeros_like(conics)
        g_opac = torch.zeros_like(opac)
        g_colors = torch.zeros_like(colors)
        bwd((tiles_x, tiles_y), (TILE, TILE),
            (_cp(ranges), _cp(point_list), _cp(means2d), _cp(conics), _cp(opac), _cp(colors),
             _cp(bg), np.int32(W), np.int32(H), np.int32(tiles_x),
             _cp(final_T), _cp(n_contrib), _cp(grad_out.contiguous()),
             _cp(g_means), _cp(g_conics), _cp(g_opac), _cp(g_colors)))
        return g_means, g_conics, g_opac, g_colors, None, None, None, None, None


# -----------------------------------------------------------------------------
# Spherical harmonics (same conventions as Kerbl et al. / pliman C++ engine)
# -----------------------------------------------------------------------------
C0 = 0.28209479177387814
C1 = 0.4886025119029199
C2 = [1.0925484305920792, -1.0925484305920792, 0.31539156525252005,
      -1.0925484305920792, 0.5462742152960396]
C3 = [-0.5900435899266435, 2.890611442640554, -0.4570457994644658,
      0.3731763325901154, -0.4570457994644658, 1.445305721320277,
      -0.5900435899266435]


def eval_sh(deg, sh, dirs):
    """sh: [N, K, 3]; dirs: [N, 3] (unit). Returns [N, 3]."""
    result = C0 * sh[:, 0]
    if deg > 0:
        x, y, z = dirs[:, 0:1], dirs[:, 1:2], dirs[:, 2:3]
        result = result - C1 * y * sh[:, 1] + C1 * z * sh[:, 2] - C1 * x * sh[:, 3]
        if deg > 1:
            xx, yy, zz = x * x, y * y, z * z
            xy, yz, xz = x * y, y * z, x * z
            result = (result +
                      C2[0] * xy * sh[:, 4] +
                      C2[1] * yz * sh[:, 5] +
                      C2[2] * (2.0 * zz - xx - yy) * sh[:, 6] +
                      C2[3] * xz * sh[:, 7] +
                      C2[4] * (xx - yy) * sh[:, 8])
            if deg > 2:
                result = (result +
                          C3[0] * y * (3 * xx - yy) * sh[:, 9] +
                          C3[1] * xy * z * sh[:, 10] +
                          C3[2] * y * (4 * zz - xx - yy) * sh[:, 11] +
                          C3[3] * z * (2 * zz - 3 * xx - 3 * yy) * sh[:, 12] +
                          C3[4] * x * (4 * zz - xx - yy) * sh[:, 13] +
                          C3[5] * z * (xx - yy) * sh[:, 14] +
                          C3[6] * x * (xx - 3 * yy) * sh[:, 15])
    return result


def quat_to_rotmat(q):
    """q: [N, 4] (w, x, y, z), already normalised."""
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = torch.stack([
        1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
        2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
        2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)
    ], dim=-1).reshape(-1, 3, 3)
    return R


# -----------------------------------------------------------------------------
# Full differentiable render
# -----------------------------------------------------------------------------
def render(params, view, active_sh, bg, near=0.01):
    means = params["means"]
    dev = means.device
    W, H = view["width"], view["height"]
    fx, fy, cx, cy = view["fx"], view["fy"], view["cx"], view["cy"]
    Rc, tc = view["R"], view["t"]

    # ---- EWA projection, written element-wise (no tiny batched matmuls) ----
    p_cam = means @ Rc.T + tc
    x, y, zr = p_cam[:, 0], p_cam[:, 1], p_cam[:, 2]
    with torch.no_grad():
        in_front = zr > near
    z = torch.where(in_front, zr, torch.ones_like(zr))
    inv_z = 1.0 / z

    # Clamp the Jacobian evaluation point to an enlarged frustum (as in 3DGS)
    lim_x = 1.3 * (0.5 * W / fx)
    lim_y = 1.3 * (0.5 * H / fy)
    txz = torch.clamp(x * inv_z, -lim_x, lim_x)
    tyz = torch.clamp(y * inv_z, -lim_y, lim_y)
    j00 = fx * inv_z
    j02 = -fx * txz * inv_z
    j11 = fy * inv_z
    j12 = -fy * tyz * inv_z
    # T = J @ Rc  (rows t0, t1)
    t0 = j00.unsqueeze(-1) * Rc[0].unsqueeze(0) + j02.unsqueeze(-1) * Rc[2].unsqueeze(0)
    t1 = j11.unsqueeze(-1) * Rc[1].unsqueeze(0) + j12.unsqueeze(-1) * Rc[2].unsqueeze(0)

    q = F.normalize(params["quats"], dim=-1)
    qw, qx, qy, qz = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    Rq = torch.stack([
        1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qw * qz), 2 * (qx * qz + qw * qy),
        2 * (qx * qy + qw * qz), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qw * qx),
        2 * (qx * qz - qw * qy), 2 * (qy * qz + qw * qx), 1 - 2 * (qx * qx + qy * qy)
    ], dim=-1).reshape(-1, 3, 3)
    s = torch.exp(params["scales"])
    # cov2d = (T R S)(T R S)^T  ->  project the scaled rotation columns first
    A0 = (t0.unsqueeze(-1) * Rq).sum(1) * s     # [N, 3] = t0^T R diag(s)
    A1 = (t1.unsqueeze(-1) * Rq).sum(1) * s
    a = (A0 * A0).sum(-1) + 0.3
    b = (A0 * A1).sum(-1)
    c = (A1 * A1).sum(-1) + 0.3
    det = a * c - b * b

    u = fx * x * inv_z + cx
    v = fy * y * inv_z + cy

    with torch.no_grad():
        mid = 0.5 * (a + c)
        lam1 = mid + torch.sqrt(torch.clamp(mid * mid - det, min=0.1))
        radius = torch.ceil(3.0 * torch.sqrt(lam1))
        tiles_x = (W + TILE - 1) // TILE
        tiles_y = (H + TILE - 1) // TILE
        min_x = torch.clamp(torch.floor((u - radius) / TILE), 0, tiles_x).long()
        max_x = torch.clamp(torch.floor((u + radius + TILE - 1) / TILE), 0, tiles_x).long()
        min_y = torch.clamp(torch.floor((v - radius) / TILE), 0, tiles_y).long()
        max_y = torch.clamp(torch.floor((v + radius + TILE - 1) / TILE), 0, tiles_y).long()
        n_tiles = (max_x - min_x) * (max_y - min_y)
        ok = in_front & (det > 0) & (n_tiles > 0) & torch.isfinite(u) & torch.isfinite(v)

    sel = torch.nonzero(ok).squeeze(-1)
    if sel.numel() == 0:
        return None
    vis_idx = sel

    # single gather of all per-splat 2D quantities
    inv_det = 1.0 / det
    packed = torch.stack([u, v, c * inv_det, -b * inv_det, a * inv_det,
                          torch.sigmoid(params["opacities"])], dim=-1)[sel]
    means2d = packed[:, 0:2].contiguous()
    if means2d.requires_grad:
        means2d.retain_grad()
    conics = packed[:, 2:5].contiguous()
    opac = packed[:, 5].contiguous()

    K = (active_sh + 1) ** 2
    if active_sh > 0:
        sh = torch.cat([params["f_dc"], params["f_rest"][:, :K - 1]], dim=1)[sel]
        dirs = F.normalize(means[sel] - view["campos"], dim=-1)
    else:
        sh = params["f_dc"][sel]
        dirs = None
    colors = torch.clamp_min(eval_sh(active_sh, sh, dirs) + 0.5, 0.0)

    # ---- tile binning + depth sort (non-differentiable) ----
    with torch.no_grad():
        mnx, mxx = min_x[sel], max_x[sel]
        mny, mxy = min_y[sel], max_y[sel]
        wx = mxx - mnx
        per = n_tiles[sel]
        total = int(per.sum().item())
        P = sel.numel()
        gid = torch.repeat_interleave(torch.arange(P, device=dev), per, output_size=total)
        offs = torch.cumsum(per, 0) - per
        local = torch.arange(total, device=dev) - offs[gid]
        wxg = wx[gid]
        tx = mnx[gid] + local % wxg
        ty = mny[gid] + local // wxg
        tile_id = ty * tiles_x + tx
        depth_bits = zr[sel].detach().float().contiguous().view(torch.int32).long()
        keys = (tile_id << 32) | depth_bits[gid]
        keys_sorted, order = torch.sort(keys)
        point_list = gid[order].int().contiguous()
        tiles_sorted = keys_sorted >> 32
        counts = torch.bincount(tiles_sorted, minlength=tiles_x * tiles_y)
        ends = torch.cumsum(counts, 0)
        starts = ends - counts
        ranges = torch.stack([starts, ends], dim=-1).int().contiguous()

    img = RasterizeGaussians.apply(means2d, conics, opac, colors, point_list, ranges, W, H, bg)
    return {
        "image": img,
        "means2d": means2d,
        "vis_idx": vis_idx,
        "radii": radius[sel].detach(),
    }


# -----------------------------------------------------------------------------
# SSIM (separable 11x11 Gaussian window, sigma = 1.5)
# -----------------------------------------------------------------------------
def _gauss_window(size=11, sigma=1.5, channels=3, device="cuda"):
    coords = torch.arange(size, dtype=torch.float32, device=device) - size // 2
    g = torch.exp(-(coords ** 2) / (2 * sigma * sigma))
    g = g / g.sum()
    return g


def _blur(x, g):
    C = x.shape[1]
    k = g.numel()
    pad = k // 2
    x = F.conv2d(x, g.view(1, 1, 1, k).expand(C, 1, 1, k), padding=(0, pad), groups=C)
    x = F.conv2d(x, g.view(1, 1, k, 1).expand(C, 1, k, 1), padding=(pad, 0), groups=C)
    return x


def ssim(img1, img2, window):
    """img: [1, 3, H, W]; window: 1-D Gaussian from _gauss_window()."""
    stack = torch.cat([img1, img2, img1 * img1, img2 * img2, img1 * img2], dim=1)
    f = _blur(stack, window)
    mu1, mu2, e11, e22, e12 = torch.split(f, img1.shape[1], dim=1)
    mu1_sq, mu2_sq, mu12 = mu1 * mu1, mu2 * mu2, mu1 * mu2
    s1, s2, s12 = e11 - mu1_sq, e22 - mu2_sq, e12 - mu12
    Ck1, Ck2 = 0.01 ** 2, 0.03 ** 2
    m = ((2 * mu12 + Ck1) * (2 * s12 + Ck2)) / ((mu1_sq + mu2_sq + Ck1) * (s1 + s2 + Ck2))
    return m.mean()


# -----------------------------------------------------------------------------
# Binary I/O
# -----------------------------------------------------------------------------
def load_splat_binary(path):
    with open(path, "rb") as f:
        raw = f.read()
    head = np.frombuffer(raw, dtype=np.int32, count=1)[0]
    off = 0
    if head == MAGIC_V2:
        N, K = np.frombuffer(raw, dtype=np.int32, count=3)[1:3]
        off = 12
    else:
        N = head
        K = 3
        off = 4
    N, K = int(N), int(K)

    def take(dtype, count):
        nonlocal off
        arr = np.frombuffer(raw, dtype=dtype, count=count, offset=off)
        off += count * 4
        return arr

    pos = take(np.float32, N * 3).reshape(N, 3)
    scl = take(np.float32, N * 3).reshape(N, 3)
    qt = take(np.float32, N * 4).reshape(N, 4)
    op = take(np.float32, N)
    sh = take(np.float32, N * K).reshape(N, K)
    if off + N * 4 <= len(raw):
        sem = take(np.int32, N)
    else:
        sem = np.zeros(N, dtype=np.int32)
    if off + N * 4 <= len(raw):
        sc = take(np.float32, N)
    else:
        sc = np.zeros(N, dtype=np.float32)
    return dict(N=N, K=K, positions=pos, scales=scl, quats=qt, opac=op, sh=sh, sem=sem, scores=sc)


def save_splat_binary(path, pos, scl, qt, op, sh, sem, sc):
    N = pos.shape[0]
    K = sh.shape[1]
    with open(path, "wb") as f:
        np.array([MAGIC_V2, N, K], dtype=np.int32).tofile(f)
        pos.astype(np.float32).tofile(f)
        scl.astype(np.float32).tofile(f)
        qt.astype(np.float32).tofile(f)
        op.astype(np.float32).tofile(f)
        sh.astype(np.float32).tofile(f)
        sem.astype(np.int32).tofile(f)
        sc.astype(np.float32).tofile(f)


# -----------------------------------------------------------------------------
# Optimizer state surgery for densification
# -----------------------------------------------------------------------------
PARAM_NAMES = ["means", "f_dc", "f_rest", "opacities", "scales", "quats"]


def _replace_params(optimizer, params, fn):
    """fn(name, tensor, is_state) -> new tensor."""
    for group in optimizer.param_groups:
        name = group["name"]
        old = group["params"][0]
        state = optimizer.state.get(old, None)
        new_p = torch.nn.Parameter(fn(name, old.data, False).contiguous().requires_grad_(True))
        if state is not None:
            del optimizer.state[old]
            state["exp_avg"] = fn(name, state["exp_avg"], True).contiguous()
            state["exp_avg_sq"] = fn(name, state["exp_avg_sq"], True).contiguous()
            optimizer.state[new_p] = state
        group["params"][0] = new_p
        params[name] = new_p


def prune_points(optimizer, params, keep_mask, stats):
    _replace_params(optimizer, params, lambda n, t, s: t[keep_mask])
    for k in list(stats.keys()):
        stats[k] = stats[k][keep_mask]


def append_points(optimizer, params, new_vals, stats):
    def fn(name, t, is_state):
        add = new_vals[name]
        if is_state:
            add = torch.zeros_like(add)
        return torch.cat([t, add], dim=0)
    _replace_params(optimizer, params, fn)
    n_new = new_vals["means"].shape[0]
    for k in list(stats.keys()):
        stats[k] = torch.cat([stats[k], torch.zeros(n_new, dtype=stats[k].dtype, device=stats[k].device)], 0)


def densify_and_prune(optimizer, params, stats, args, extent, max_screen_size, max_gaussians):
    grads = stats["grad_accum"] / stats["denom"].clamp_min(1.0)
    grads[~torch.isfinite(grads)] = 0.0
    N0 = params["means"].shape[0]
    with torch.no_grad():
        scales = torch.exp(params["scales"])
        max_s = scales.max(dim=1).values
        high = grads >= args.grad_threshold
        clone_mask = high & (max_s <= args.percent_dense * extent)
        split_mask = high & (max_s > args.percent_dense * extent)

        budget = max(0, max_gaussians - N0)
        n_clone = int(clone_mask.sum())
        n_split = int(split_mask.sum())
        if n_clone + n_split > budget:
            # keep the highest-gradient candidates that fit the budget
            cand = torch.nonzero(clone_mask | split_mask).squeeze(-1)
            g = grads[cand]
            top = cand[torch.argsort(g, descending=True)[:budget]]
            keep = torch.zeros_like(clone_mask)
            keep[top] = True
            clone_mask &= keep
            split_mask &= keep

        new_vals = None
        if clone_mask.any():
            new_vals = {n: params[n].data[clone_mask] for n in PARAM_NAMES}

        split_vals = None
        if split_mask.any():
            Ns = 2
            sel = split_mask
            s = scales[sel].repeat(Ns, 1)
            samples = torch.randn_like(s) * s
            Rq = quat_to_rotmat(F.normalize(params["quats"].data[sel], dim=-1)).repeat(Ns, 1, 1)
            new_means = (Rq @ samples.unsqueeze(-1)).squeeze(-1) + params["means"].data[sel].repeat(Ns, 1)
            split_vals = {
                "means": new_means,
                "scales": torch.log(s / (0.8 * Ns)),
                "quats": params["quats"].data[sel].repeat(Ns, 1),
                "f_dc": params["f_dc"].data[sel].repeat(Ns, 1, 1),
                "f_rest": params["f_rest"].data[sel].repeat(Ns, 1, 1),
                "opacities": params["opacities"].data[sel].repeat(Ns),
            }

    if new_vals is not None:
        append_points(optimizer, params, new_vals, stats)
    if split_vals is not None:
        n_before = params["means"].shape[0]
        append_points(optimizer, params, split_vals, stats)
        # remove the split parents
        prune_parent = torch.zeros(params["means"].shape[0], dtype=torch.bool, device=params["means"].device)
        prune_parent[:N0][split_mask] = True
        prune_points(optimizer, params, ~prune_parent, stats)

    with torch.no_grad():
        opac = torch.sigmoid(params["opacities"])
        prune = opac < args.min_opacity
        if max_screen_size is not None:
            big_vs = stats["max_radii"] > max_screen_size
            big_ws = torch.exp(params["scales"]).max(dim=1).values > args.max_scale * extent
            prune = prune | big_vs | big_ws
        if params["means"].shape[0] - int(prune.sum()) < 16:
            prune = torch.zeros_like(prune)
    prune_points(optimizer, params, ~prune, stats)

    stats["grad_accum"].zero_()
    stats["denom"].zero_()
    stats["max_radii"].zero_()


def reset_opacity(optimizer, params):
    with torch.no_grad():
        new_op = torch.min(params["opacities"].data, torch.logit(torch.tensor(0.01, device=params["opacities"].device)))

    def fn(name, t, is_state):
        if name != "opacities":
            return t
        return torch.zeros_like(t) if is_state else new_op
    _replace_params(optimizer, params, fn)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def parse_int_safe(val):
    if val is None or val == "None" or val == "NA":
        return None
    return int(float(val))


def load_view(spec, downsample, device):
    img_path = spec["image_path"]
    pil = Image.open(img_path)
    spec_w = int(spec.get("width", pil.size[0]))
    spec_h = int(spec.get("height", pil.size[1]))
    # Honour EXIF orientation only when the calibrated size requires it
    if pil.size != (spec_w, spec_h):
        rot = ImageOps.exif_transpose(pil)
        if rot.size == (spec_w, spec_h):
            pil = rot
    pil = pil.convert("RGB")
    orig_w, orig_h = pil.size
    ds = max(1.0, float(downsample))
    tw = max(16, int(round(orig_w / ds)))
    th = max(16, int(round(orig_h / ds)))
    if (tw, th) != (orig_w, orig_h):
        pil = pil.resize((tw, th), Image.Resampling.LANCZOS)
    arr = np.asarray(pil, dtype=np.uint8).copy()
    sx = tw / float(spec_w)
    sy = th / float(spec_h)
    R = torch.tensor(np.asarray(spec["R"], dtype=np.float32).reshape(3, 3), device=device)
    t = torch.tensor(np.asarray(spec["t"], dtype=np.float32).reshape(3), device=device)
    campos = -(R.T @ t)
    return {
        "target_u8": torch.from_numpy(arr).to(device),   # [H, W, 3] uint8
        "width": tw, "height": th,
        "fx": float(spec["fx"]) * sx, "fy": float(spec["fy"]) * sy,
        "cx": float(spec["cx"]) * sx, "cy": float(spec["cy"]) * sy,
        "R": R, "t": t, "campos": campos,
        "label": spec.get("label", os.path.basename(img_path)),
    }


def main():
    p = argparse.ArgumentParser(description="pliman GPU 3D Gaussian Splatting trainer")
    p.add_argument("--splat_in", required=True)
    p.add_argument("--cameras_json", required=True)
    p.add_argument("--splat_out", required=True)
    p.add_argument("--iterations", type=parse_int_safe, default=7000)
    p.add_argument("--lr_pos", type=float, default=0.00016)
    p.add_argument("--lr_pos_final", type=float, default=0.0000016)
    p.add_argument("--lr_scale", type=float, default=0.005)
    p.add_argument("--lr_opac", type=float, default=0.05)
    p.add_argument("--lr_sh", type=float, default=0.0025)
    p.add_argument("--lr_rot", type=float, default=0.001)
    p.add_argument("--densify_interval", type=parse_int_safe, default=100)
    p.add_argument("--densify_from", type=parse_int_safe, default=None)
    p.add_argument("--densify_until", type=parse_int_safe, default=None)
    p.add_argument("--opacity_reset_interval", type=parse_int_safe, default=3000)
    p.add_argument("--grad_threshold", type=float, default=0.0002)
    p.add_argument("--percent_dense", "--scale_threshold", dest="percent_dense", type=float, default=0.01)
    p.add_argument("--min_opacity", type=float, default=0.005)
    p.add_argument("--max_scale", type=float, default=0.1)
    p.add_argument("--max_screen_size", type=float, default=20.0)
    p.add_argument("--max_gaussians", type=parse_int_safe, default=1000000)
    p.add_argument("--downsample", type=float, default=4)
    p.add_argument("--sh_degree", type=parse_int_safe, default=3)
    p.add_argument("--sh_interval", type=parse_int_safe, default=1000)
    p.add_argument("--lambda_dssim", type=float, default=0.2)
    p.add_argument("--white_background", type=parse_int_safe, default=0)
    p.add_argument("--random_background", type=parse_int_safe, default=0)
    p.add_argument("--seed", type=parse_int_safe, default=42)
    p.add_argument("--log_interval", type=parse_int_safe, default=None)
    p.add_argument("--stats_json", type=str, default=None)
    args = p.parse_args()

    if not torch.cuda.is_available():
        sys.stderr.write("CUDA is not available on this system. GPU training requires an NVIDIA GPU.\n")
        sys.exit(1)

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    device = torch.device("cuda:0")
    print(f"[INIT] GPU: {torch.cuda.get_device_name(device)} (CUDA {torch.version.cuda})", flush=True)

    try:
        get_kernels()
    except Exception as e:
        sys.stderr.write(f"Failed to compile CUDA rasterizer with CuPy/NVRTC: {e}\n")
        sys.exit(2)

    # ---- Load Gaussians ----
    data = load_splat_binary(args.splat_in)
    N = data["N"]
    max_deg = max(0, min(3, int(args.sh_degree)))
    K_full = (max_deg + 1) ** 2
    sh_in = data["sh"]
    K_in = sh_in.shape[1] // 3
    sh_r = np.zeros((N, K_full, 3), dtype=np.float32)
    k_copy = min(K_in, K_full)
    sh_r[:, :k_copy, :] = sh_in[:, :k_copy * 3].reshape(N, k_copy, 3)
    print(f"[INIT] Loaded {N:,} initial Gaussians (SH degree target: {max_deg}).", flush=True)

    def P(x):
        return torch.nn.Parameter(torch.tensor(x, dtype=torch.float32, device=device).contiguous())

    params = {
        "means": P(data["positions"]),
        "f_dc": P(sh_r[:, :1, :]),
        "f_rest": P(sh_r[:, 1:, :]),
        "opacities": P(data["opac"]),
        "scales": P(data["scales"]),
        "quats": P(data["quats"]),
    }
    semantic_ids = data["sem"]

    # ---- Load cameras / views ----
    with open(args.cameras_json, "r", encoding="utf-8") as f:
        cam_specs = json.load(f)
    if isinstance(cam_specs, dict):
        cam_specs = list(cam_specs.values())
    views = []
    for spec in cam_specs:
        if not os.path.exists(spec["image_path"]):
            continue
        views.append(load_view(spec, args.downsample, device))
    if not views:
        sys.stderr.write("Error: no training views could be loaded.\n")
        sys.exit(1)
    print(f"[INIT] Cached {len(views)} training views at {views[0]['width']}x{views[0]['height']} px "
          f"(downsample {args.downsample:g}x).", flush=True)

    # Scene extent (NeRF++ normalisation radius, as in 3DGS)
    centers = torch.stack([v["campos"] for v in views])
    c_mean = centers.mean(0)
    extent = float(torch.linalg.norm(centers - c_mean, dim=1).max().item() * 1.1)
    if not math.isfinite(extent) or extent <= 0:
        extent = 1.0
    print(f"[INIT] Scene extent (camera radius): {extent:.4f}", flush=True)

    iters = int(args.iterations)
    densify_from = args.densify_from if args.densify_from is not None else min(500, max(1, iters // 10))
    densify_until = args.densify_until if args.densify_until is not None else iters // 2
    opac_reset = args.opacity_reset_interval if args.opacity_reset_interval and args.opacity_reset_interval > 0 else 10 ** 9

    optimizer = torch.optim.Adam([
        {"params": [params["means"]], "lr": args.lr_pos * extent, "name": "means"},
        {"params": [params["f_dc"]], "lr": args.lr_sh, "name": "f_dc"},
        {"params": [params["f_rest"]], "lr": args.lr_sh / 20.0, "name": "f_rest"},
        {"params": [params["opacities"]], "lr": args.lr_opac, "name": "opacities"},
        {"params": [params["scales"]], "lr": args.lr_scale, "name": "scales"},
        {"params": [params["quats"]], "lr": args.lr_rot, "name": "quats"},
    ], lr=0.0, eps=1e-15)

    lr0 = args.lr_pos * extent
    lr1 = max(1e-12, args.lr_pos_final * extent)

    stats = {
        "grad_accum": torch.zeros(N, device=device),
        "denom": torch.zeros(N, device=device),
        "max_radii": torch.zeros(N, device=device),
    }

    bg_base = torch.tensor([1.0, 1.0, 1.0] if args.white_background else [0.0, 0.0, 0.0], device=device)
    window = _gauss_window(device=device)
    log_interval = args.log_interval or max(1, min(100, iters // 50))

    order = rng.permutation(len(views))
    ptr = 0
    recent_l, recent_psnr = [], []
    history = []
    t0 = time.time()
    print(f"[START] {iters} iterations | densify {densify_from}-{densify_until} every {args.densify_interval} "
          f"| opacity reset every {opac_reset if opac_reset < 10**9 else 'never'}", flush=True)

    for step in range(1, iters + 1):
        # position LR schedule (log-linear decay)
        tt = min(1.0, (step - 1) / max(1, iters - 1))
        lr_now = math.exp(math.log(lr0) * (1 - tt) + math.log(lr1) * tt)
        for g in optimizer.param_groups:
            if g["name"] == "means":
                g["lr"] = lr_now

        active_sh = min(max_deg, (step - 1) // max(1, args.sh_interval))

        if ptr >= len(order):
            order = rng.permutation(len(views))
            ptr = 0
        view = views[order[ptr]]
        ptr += 1

        bg = torch.rand(3, device=device) if args.random_background else bg_base
        out = render(params, view, active_sh, bg)
        if out is None:
            continue
        img = out["image"]
        gt = view["target_u8"].float() / 255.0

        l1 = torch.abs(img - gt).mean()
        if args.lambda_dssim > 0:
            ss = ssim(img.permute(2, 0, 1).unsqueeze(0), gt.permute(2, 0, 1).unsqueeze(0), window)
            loss = (1.0 - args.lambda_dssim) * l1 + args.lambda_dssim * (1.0 - ss)
        else:
            loss = l1
        loss.backward()

        with torch.no_grad():
            mse = torch.mean((img - gt) ** 2).item()
            psnr = 10.0 * math.log10(1.0 / max(mse, 1e-10))
            recent_l.append(loss.item())
            recent_psnr.append(psnr)
            if len(recent_l) > 50:
                recent_l.pop(0)
                recent_psnr.pop(0)

            # densification statistics
            if step < densify_until:
                vis = out["vis_idx"]
                g2 = out["means2d"].grad
                if g2 is not None:
                    gn = torch.sqrt((g2[:, 0] * view["width"] * 0.5) ** 2 + (g2[:, 1] * view["height"] * 0.5) ** 2)
                    stats["grad_accum"].index_add_(0, vis, gn)
                    stats["denom"].index_add_(0, vis, torch.ones_like(gn))
                    stats["max_radii"][vis] = torch.maximum(stats["max_radii"][vis], out["radii"])

        optimizer.step()
        optimizer.zero_grad(set_to_none=True)

        if step < densify_until:
            if step > densify_from and step % args.densify_interval == 0:
                n_before = params["means"].shape[0]
                size_thr = args.max_screen_size if step > opac_reset else None
                densify_and_prune(optimizer, params, stats, args, extent, size_thr, args.max_gaussians)
                n_after = params["means"].shape[0]
                if len(semantic_ids) != n_after:
                    semantic_ids = np.zeros(n_after, dtype=np.int32)
                print(f"* Step {step}/{iters} | Densify: {n_before:,} -> {n_after:,} Gaussians", flush=True)
            if step % opac_reset == 0 or (args.white_background and step == densify_from):
                reset_opacity(optimizer, params)
                print(f"* Step {step}/{iters} | Opacity reset", flush=True)

        if step % log_interval == 0 or step == iters:
            el = time.time() - t0
            ips = step / max(1e-3, el)
            rem = (iters - step) / max(1e-3, ips)
            eta = f"{rem / 60:.1f}m" if rem >= 60 else f"{int(rem)}s"
            ml = float(np.mean(recent_l)) if recent_l else 0.0
            mp = float(np.mean(recent_psnr)) if recent_psnr else 0.0
            history.append({"step": step, "loss": ml, "psnr": mp, "n": int(params["means"].shape[0])})
            print(f"* Step {step}/{iters} [{int(100 * step / iters)}%] | {ips:.1f} it/s | ETA: {eta} | "
                  f"Loss: {ml:.4f} | PSNR: {mp:.2f} dB | SH: {active_sh} | Gaussians: {params['means'].shape[0]:,}",
                  flush=True)

    total = time.time() - t0
    n_final = params["means"].shape[0]
    print(f"* GPU training finished in {total:.1f}s | Final Gaussians: {n_final:,} | "
          f"PSNR: {float(np.mean(recent_psnr)) if recent_psnr else 0:.2f} dB", flush=True)

    with torch.no_grad():
        q = F.normalize(params["quats"].data, dim=-1)
        sh_out = torch.cat([params["f_dc"].data, params["f_rest"].data], dim=1)  # [N, K, 3]
        sh_out = sh_out.reshape(n_final, -1).cpu().numpy()
        if len(semantic_ids) != n_final:
            semantic_ids = np.zeros(n_final, dtype=np.int32)
        save_splat_binary(
            args.splat_out,
            params["means"].data.cpu().numpy(),
            params["scales"].data.cpu().numpy(),
            q.cpu().numpy(),
            params["opacities"].data.cpu().numpy(),
            sh_out,
            semantic_ids,
            np.zeros(n_final, dtype=np.float32),
        )

    if args.stats_json:
        with open(args.stats_json, "w", encoding="utf-8") as f:
            json.dump({"extent": extent, "seconds": total, "history": history,
                       "sh_degree": max_deg, "num_gaussians": int(n_final)}, f)


if __name__ == "__main__":
    main()
