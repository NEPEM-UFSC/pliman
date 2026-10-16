#!/usr/bin/env python3
"""
GPU-Accelerated 3D Gaussian Splatting Training Engine for pliman
Utilizes PyTorch CUDA and CuPy NVRTC for high-throughput differentiable rasterization.
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

# Automatically resolve CUDA headers from nvidia-cuda-runtime if needed
try:
    import nvidia.cuda_runtime
    cuda_rt_dir = nvidia.cuda_runtime.__path__[0]
    if "CUDA_PATH" not in os.environ or not os.path.exists(os.environ["CUDA_PATH"]):
        os.environ["CUDA_PATH"] = cuda_rt_dir
except Exception:
    pass

import torch
import cupy
from PIL import Image

# -----------------------------------------------------------------------------
# CUDA Kernel for Tile-Based Forward Rasterization & Photometric Backward Pass
# -----------------------------------------------------------------------------
CUDA_RASTER_SRC = r"""
extern "C" {

__device__ __forceinline__ float clamp_val(float v, float lo, float hi) {
    return fminf(fmaxf(v, lo), hi);
}

__global__ void rasterize_and_backward_kernel(
    const int* __restrict__ tile_ranges,       // [num_tiles + 1]
    const int* __restrict__ item_to_splat_idx, // [total_items]
    const float* __restrict__ splat_u,         // [num_vis]
    const float* __restrict__ splat_v,
    const float* __restrict__ splat_inv_xx,
    const float* __restrict__ splat_inv_xy,
    const float* __restrict__ splat_inv_yy,
    const float* __restrict__ splat_opacity,
    const float* __restrict__ splat_r,
    const float* __restrict__ splat_g,
    const float* __restrict__ splat_b,
    const int*   __restrict__ splat_g_idx,
    const float* __restrict__ splat_cam_x,
    const float* __restrict__ splat_cam_y,
    const float* __restrict__ splat_cam_z,
    const float* __restrict__ gt_img,          // [3, H, W]
    float*       __restrict__ grad_pos,        // [N, 3]
    float*       __restrict__ grad_scale,      // [N, 3]
    float*       __restrict__ grad_opac,       // [N]
    float*       __restrict__ grad_sh,         // [N, 3]
    float*       __restrict__ grad_2d_norm,    // [N]
    float*       __restrict__ tile_loss,       // [num_tiles]
    int tiles_x, int tiles_y,
    int width, int height,
    float fx, float fy,
    const float* __restrict__ cam_R_T          // [3, 3] world to cam transpose
) {
    int tx = blockIdx.x;
    int ty = blockIdx.y;
    int tile_idx = ty * tiles_x + tx;

    int px = tx * 16 + threadIdx.x;
    int py = ty * 16 + threadIdx.y;

    int start_item = tile_ranges[tile_idx];
    int end_item = tile_ranges[tile_idx + 1];

    float loss_thread = 0.0f;
    float inv_pix = 1.0f / (float)(width * height);

    if (px < width && py < height && start_item < end_item) {
        float T = 1.0f;
        float acc_r = 0.0f, acc_g = 0.0f, acc_b = 0.0f;

        const int MAX_C = 48;
        int contrib_sidx[MAX_C];
        float contrib_alpha[MAX_C];
        float contrib_weight[MAX_C];
        float contrib_dx[MAX_C];
        float contrib_dy[MAX_C];
        int num_c = 0;

        for (int k = start_item; k < end_item; ++k) {
            int s = item_to_splat_idx[k];
            float u = splat_u[s];
            float v = splat_v[s];
            float dx = (float)px - u;
            float dy = (float)py - v;

            float power = -0.5f * (dx * dx * splat_inv_xx[s] + 
                                  2.0f * dx * dy * splat_inv_xy[s] + 
                                  dy * dy * splat_inv_yy[s]);
            if (power > 0.0f) continue;

            float alpha = fminf(0.99f, splat_opacity[s] * __expf(power));
            if (alpha < (1.0f / 255.0f)) continue;

            float weight = alpha * T;
            acc_r += splat_r[s] * weight;
            acc_g += splat_g[s] * weight;
            acc_b += splat_b[s] * weight;

            if (num_c < MAX_C) {
                contrib_sidx[num_c] = s;
                contrib_alpha[num_c] = alpha;
                contrib_weight[num_c] = weight;
                contrib_dx[num_c] = dx;
                contrib_dy[num_c] = dy;
                num_c++;
            }

            T *= (1.0f - alpha);
            if (T < 0.001f) break;
        }

        int pix_idx = py * width + px;
        float gt_r = gt_img[pix_idx];
        float gt_g = gt_img[width * height + pix_idx];
        float gt_b = gt_img[2 * width * height + pix_idx];

        float err_r = acc_r - gt_r;
        float err_g = acc_g - gt_g;
        float err_b = acc_b - gt_b;

        loss_thread = (fabsf(err_r) + fabsf(err_g) + fabsf(err_b)) * inv_pix;

        float dL_r = clamp_val(err_r, -1.0f, 1.0f) * inv_pix;
        float dL_g = clamp_val(err_g, -1.0f, 1.0f) * inv_pix;
        float dL_b = clamp_val(err_b, -1.0f, 1.0f) * inv_pix;

        const float SH_C0 = 0.28209479177387814f;

        float accum_color_grad = 0.0f;
        for (int i = num_c - 1; i >= 0; --i) {
            int s = contrib_sidx[i];
            int g = splat_g_idx[s];
            float alpha = contrib_alpha[i];
            float weight = contrib_weight[i];
            float dx = contrib_dx[i];
            float dy = contrib_dy[i];

            float r = splat_r[s];
            float g_col = splat_g[s];
            float b = splat_b[s];

            float d_sh_r = dL_r * weight * SH_C0;
            float d_sh_g = dL_g * weight * SH_C0;
            float d_sh_b = dL_b * weight * SH_C0;

            atomicAdd(&grad_sh[g * 3 + 0], d_sh_r);
            atomicAdd(&grad_sh[g * 3 + 1], d_sh_g);
            atomicAdd(&grad_sh[g * 3 + 2], d_sh_b);

            float dot_grad = r * dL_r + g_col * dL_g + b * dL_b;
            float d_alpha = (weight / fmaxf(1e-5f, alpha)) * dot_grad - accum_color_grad;
            accum_color_grad += weight * dot_grad;

            float d_u = d_alpha * alpha * (dx * splat_inv_xx[s] + dy * splat_inv_xy[s]);
            float d_v = d_alpha * alpha * (dx * splat_inv_xy[s] + dy * splat_inv_yy[s]);

            float tz = splat_cam_z[s];
            float tz2 = tz * tz;
            float dt_x = d_u * (fx / tz);
            float dt_y = d_v * (fy / tz);
            float dt_z = d_u * (-(fx * splat_cam_x[s]) / tz2) + d_v * (-(fy * splat_cam_y[s]) / tz2);

            float dP_x = cam_R_T[0] * dt_x + cam_R_T[1] * dt_y + cam_R_T[2] * dt_z;
            float dP_y = cam_R_T[3] * dt_x + cam_R_T[4] * dt_y + cam_R_T[5] * dt_z;
            float dP_z = cam_R_T[6] * dt_x + cam_R_T[7] * dt_y + cam_R_T[8] * dt_z;

            atomicAdd(&grad_pos[g * 3 + 0], dP_x);
            atomicAdd(&grad_pos[g * 3 + 1], dP_y);
            atomicAdd(&grad_pos[g * 3 + 2], dP_z);

            float g2d = sqrtf(dP_x * dP_x + dP_y * dP_y + dP_z * dP_z);
            atomicAdd(&grad_2d_norm[g], g2d);

            float d_opac = d_alpha * (1.0f - splat_opacity[s]);
            atomicAdd(&grad_opac[g], d_opac);

            float d_s = d_alpha * 0.001f;
            atomicAdd(&grad_scale[g * 3 + 0], d_s);
            atomicAdd(&grad_scale[g * 3 + 1], d_s);
            atomicAdd(&grad_scale[g * 3 + 2], d_s);
        }
    }
    atomicAdd(&tile_loss[tile_idx], loss_thread);
}

}
"""

def load_splat_binary(file_path, device):
    """Load splat attributes from binary format."""
    with open(file_path, "rb") as f:
        N = np.fromfile(f, dtype=np.int32, count=1)[0]
        pos = np.fromfile(f, dtype=np.float32, count=N * 3).reshape(N, 3)
        scales_log = np.fromfile(f, dtype=np.float32, count=N * 3).reshape(N, 3)
        quaternions = np.fromfile(f, dtype=np.float32, count=N * 4).reshape(N, 4)
        opacities_logit = np.fromfile(f, dtype=np.float32, count=N)
        sh_coeffs = np.fromfile(f, dtype=np.float32, count=N * 3).reshape(N, 3)
        semantic_ids = np.fromfile(f, dtype=np.int32, count=N) if f.peek(4) else np.zeros(N, dtype=np.int32)
        semantic_scores = np.fromfile(f, dtype=np.float32, count=N) if f.peek(4) else np.zeros(N, dtype=np.float32)

    return {
        "num_gaussians": int(N),
        "positions": torch.tensor(pos, dtype=torch.float32, device=device),
        "scales_log": torch.tensor(scales_log, dtype=torch.float32, device=device),
        "quaternions": torch.tensor(quaternions, dtype=torch.float32, device=device),
        "opacities_logit": torch.tensor(opacities_logit, dtype=torch.float32, device=device),
        "sh_coeffs": torch.tensor(sh_coeffs, dtype=torch.float32, device=device),
        "semantic_ids": torch.tensor(semantic_ids, dtype=torch.int32, device=device),
        "semantic_scores": torch.tensor(semantic_scores, dtype=torch.float32, device=device),
    }

def save_splat_binary(splat, file_path):
    """Save splat attributes to binary format."""
    N = int(splat["num_gaussians"])
    with open(file_path, "wb") as f:
        np.int32(N).tofile(f)
        splat["positions"].detach().cpu().numpy().astype(np.float32).tofile(f)
        splat["scales_log"].detach().cpu().numpy().astype(np.float32).tofile(f)
        splat["quaternions"].detach().cpu().numpy().astype(np.float32).tofile(f)
        splat["opacities_logit"].detach().cpu().numpy().astype(np.float32).tofile(f)
        splat["sh_coeffs"].detach().cpu().numpy().astype(np.float32).tofile(f)
        splat["semantic_ids"].detach().cpu().numpy().astype(np.int32).tofile(f)
        splat["semantic_scores"].detach().cpu().numpy().astype(np.float32).tofile(f)

def quat_to_rot_matrix(q):
    """q: [N, 4] (w, x, y, z)"""
    q = F.normalize(q, dim=-1)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = torch.stack([
        1.0 - 2.0 * (y*y + z*z), 2.0 * (x*y - w*z), 2.0 * (x*z + w*y),
        2.0 * (x*y + w*z), 1.0 - 2.0 * (x*x + z*z), 2.0 * (y*z - w*x),
        2.0 * (x*z - w*y), 2.0 * (y*z + w*x), 1.0 - 2.0 * (x*x + y*y)
    ], dim=-1).reshape(-1, 3, 3)
    return R

def parse_int_safe(val):
    if val is None or val == "None":
        return None
    return int(float(val))

def main():
    parser = argparse.ArgumentParser(description="GPU 3D Gaussian Splatting Trainer")
    parser.add_argument("--splat_in", type=str, required=True, help="Input splat binary file")
    parser.add_argument("--cameras_json", type=str, required=True, help="Cameras JSON specification")
    parser.add_argument("--splat_out", type=str, required=True, help="Output splat binary file")
    parser.add_argument("--iterations", type=parse_int_safe, default=500, help="Training iterations")
    parser.add_argument("--lr_pos", type=float, default=0.00016, help="Position learning rate")
    parser.add_argument("--lr_scale", type=float, default=0.005, help="Scale learning rate")
    parser.add_argument("--lr_opac", type=float, default=0.05, help="Opacity learning rate")
    parser.add_argument("--lr_sh", type=float, default=0.0025, help="SH color learning rate")
    parser.add_argument("--densify_interval", type=parse_int_safe, default=100, help="Densification interval")
    parser.add_argument("--densify_until", type=parse_int_safe, default=None, help="Densify until step")
    parser.add_argument("--grad_threshold", type=float, default=0.0002, help="Densification grad threshold")
    parser.add_argument("--scale_threshold", type=float, default=0.01, help="Scale split threshold")
    parser.add_argument("--min_opacity", type=float, default=0.05, help="Pruning minimum opacity")
    parser.add_argument("--max_scale", type=float, default=0.5, help="Pruning maximum scale")
    parser.add_argument("--max_gaussians", type=parse_int_safe, default=200000, help="Maximum number of Gaussians")
    parser.add_argument("--downsample", type=parse_int_safe, default=4, help="Image downsampling factor")
    args = parser.parse_args()

    if not torch.cuda.is_available():
        sys.stderr.write("CUDA is not available on this system. GPU training requires an NVIDIA GPU.\n")
        sys.exit(1)

    device = torch.device("cuda:0")
    gpu_name = torch.cuda.get_device_name(device)
    print(f"[INIT] Using GPU: {gpu_name} (CUDA {torch.version.cuda})")

    # Load splat
    splat = load_splat_binary(args.splat_in, device)
    N = splat["num_gaussians"]
    print(f"[INIT] Loaded {N:,} initial 3D Gaussians.")

    # Load cameras
    with open(args.cameras_json, "r", encoding="utf-8") as f:
        cam_specs = json.load(f)

    if len(cam_specs) == 0:
        sys.stderr.write("Error: No camera views found in cameras JSON.\n")
        sys.exit(1)

    # Preload and downsample training views onto GPU memory / system cache
    print(f"[INIT] Loading and caching {len(cam_specs)} training views (downsample: {args.downsample}x)...")
    training_views = []
    for spec in cam_specs:
        img_path = spec["image_path"]
        if not os.path.exists(img_path):
            continue

        pil_img = Image.open(img_path).convert("RGB")
        orig_w, orig_h = pil_img.size
        
        target_w = max(32, orig_w // args.downsample)
        target_h = max(32, orig_h // args.downsample)
        
        # Ensure dimensions are multiples of 16 for clean tiling
        target_w = (target_w // 16) * 16
        target_h = (target_h // 16) * 16

        pil_resized = pil_img.resize((target_w, target_h), Image.Resampling.BILINEAR)
        img_arr = np.array(pil_resized, dtype=np.float32) / 255.0
        # Transpose to [3, H, W]
        img_tensor = torch.tensor(img_arr.transpose(2, 0, 1), dtype=torch.float32, device=device)

        scale_w = float(target_w) / float(spec.get("width", orig_w))
        scale_h = float(target_h) / float(spec.get("height", orig_h))

        cam_R = torch.tensor(spec["R"], dtype=torch.float32, device=device)
        cam_t = torch.tensor(spec["t"], dtype=torch.float32, device=device)

        training_views.append({
            "target": img_tensor,
            "width": target_w,
            "height": target_h,
            "fx": spec["fx"] * scale_w,
            "fy": spec["fy"] * scale_h,
            "cx": spec["cx"] * scale_w,
            "cy": spec["cy"] * scale_h,
            "R": cam_R,
            "t": cam_t,
            "R_T": cam_R.T.contiguous()
        })

    num_views = len(training_views)
    print(f"[INIT] Prepared {num_views} training views.")

    # Compile NVRTC CUDA Kernel
    raw_kernel = cupy.RawKernel(CUDA_RASTER_SRC, "rasterize_and_backward_kernel")

    # Parameters & Optimizer Moments
    pos = splat["positions"]
    scales_log = splat["scales_log"]
    quaternions = splat["quaternions"]
    opacities_logit = splat["opacities_logit"]
    sh_coeffs = splat["sh_coeffs"]

    m_pos = torch.zeros_like(pos); v_pos = torch.zeros_like(pos)
    m_scale = torch.zeros_like(scales_log); v_scale = torch.zeros_like(scales_log)
    m_opac = torch.zeros_like(opacities_logit); v_opac = torch.zeros_like(opacities_logit)
    m_sh = torch.zeros_like(sh_coeffs); v_sh = torch.zeros_like(sh_coeffs)

    accum_grad_2d = torch.zeros(N, device=device)
    accum_count = 0
    recent_losses = []

    densify_until = args.densify_until if args.densify_until is not None else int(args.iterations * 0.8)
    TILE_SIZE = 16
    start_time = time.time()

    beta1 = 0.9
    beta2 = 0.999
    eps = 1e-8

    print(f"[START] Beginning {args.iterations} GPU training steps...")

    for step in range(1, args.iterations + 1):
        v_idx = (step - 1) % num_views
        view = training_views[v_idx]

        w, h = view["width"], view["height"]
        fx, fy = view["fx"], view["fy"]
        cx, cy = view["cx"], view["cy"]
        cam_R = view["R"]
        cam_t = view["t"]
        cam_R_T = view["R_T"]
        gt_img = view["target"]

        tiles_x = (w + TILE_SIZE - 1) // TILE_SIZE
        tiles_y = (h + TILE_SIZE - 1) // TILE_SIZE
        total_tiles = tiles_x * tiles_y

        # 1. 3D-to-2D View Projection on GPU
        p_cam = torch.matmul(pos, cam_R.T) + cam_t
        tz = p_cam[:, 2]
        vis = (tz > 0.1)

        u = fx * (p_cam[:, 0] / tz) + cx
        v = fy * (p_cam[:, 1] / tz) + cy
        vis = vis & (u >= -32) & (u <= w + 32) & (v >= -32) & (v <= h + 32)

        vis_indices = torch.nonzero(vis).squeeze(-1)
        num_vis = vis_indices.shape[0]

        if num_vis > 0:
            u_vis = u[vis_indices]
            v_vis = v[vis_indices]
            tz_vis = tz[vis_indices]
            px_vis = p_cam[vis_indices, 0]
            py_vis = p_cam[vis_indices, 1]

            s_vis = torch.exp(scales_log[vis_indices])
            r_max = torch.max(s_vis, dim=-1)[0]
            proj_rad = (fx * r_max / tz_vis) * 3.0
            radius_vis = torch.clamp(torch.ceil(proj_rad), min=1.0, max=40.0)

            var = torch.clamp((fx * r_max / tz_vis)**2 + 0.3, min=0.1)
            inv_var = 1.0 / var
            inv_xx = inv_var
            inv_xy = torch.zeros_like(inv_var)
            inv_yy = inv_var

            opac_vis = torch.sigmoid(opacities_logit[vis_indices])
            sh_vis = sh_coeffs[vis_indices]
            SH_C0 = 0.28209479177387814
            r_vis = torch.clamp(0.5 + SH_C0 * sh_vis[:, 0], 0.0, 1.0)
            g_vis = torch.clamp(0.5 + SH_C0 * sh_vis[:, 1], 0.0, 1.0)
            b_vis = torch.clamp(0.5 + SH_C0 * sh_vis[:, 2], 0.0, 1.0)

            # 2. Tile binning & depth sorting
            min_tx = torch.clamp(((u_vis - radius_vis) / TILE_SIZE).long(), 0, tiles_x - 1)
            max_tx = torch.clamp(((u_vis + radius_vis) / TILE_SIZE).long(), 0, tiles_x - 1)
            min_ty = torch.clamp(((v_vis - radius_vis) / TILE_SIZE).long(), 0, tiles_y - 1)
            max_ty = torch.clamp(((v_vis + radius_vis) / TILE_SIZE).long(), 0, tiles_y - 1)

            rep_x = (max_tx - min_tx + 1)
            rep_y = (max_ty - min_ty + 1)
            num_tiles_per_s = rep_x * rep_y

            depth_order = torch.argsort(tz_vis)
            s_idx_expanded = torch.repeat_interleave(depth_order, num_tiles_per_s[depth_order])

            tile_ranges = torch.zeros(total_tiles + 1, dtype=torch.int32, device=device)
            tile_ranges[-1] = s_idx_expanded.shape[0]
            items_per_tile = s_idx_expanded.shape[0] // total_tiles
            tile_ranges[1:-1] = (torch.arange(1, total_tiles, device=device, dtype=torch.int32) * items_per_tile)

            # Gradients
            grad_pos = torch.zeros_like(pos)
            grad_scale = torch.zeros_like(scales_log)
            grad_opac = torch.zeros_like(opacities_logit)
            grad_sh = torch.zeros_like(sh_coeffs)
            grad_2d_norm = torch.zeros(pos.shape[0], device=device)
            tile_loss = torch.zeros(total_tiles, device=device)

            # Zero-copy DLPack to CuPy
            cp_tile_ranges = cupy.from_dlpack(tile_ranges)
            cp_items = cupy.from_dlpack(s_idx_expanded.contiguous().to(torch.int32))
            cp_u = cupy.from_dlpack(u_vis.contiguous())
            cp_v = cupy.from_dlpack(v_vis.contiguous())
            cp_inv_xx = cupy.from_dlpack(inv_xx.contiguous())
            cp_inv_xy = cupy.from_dlpack(inv_xy.contiguous())
            cp_inv_yy = cupy.from_dlpack(inv_yy.contiguous())
            cp_opac = cupy.from_dlpack(opac_vis.contiguous())
            cp_r = cupy.from_dlpack(r_vis.contiguous())
            cp_g = cupy.from_dlpack(g_vis.contiguous())
            cp_b = cupy.from_dlpack(b_vis.contiguous())
            cp_g_idx = cupy.from_dlpack(vis_indices.contiguous().to(torch.int32))
            cp_cam_x = cupy.from_dlpack(px_vis.contiguous())
            cp_cam_y = cupy.from_dlpack(py_vis.contiguous())
            cp_cam_z = cupy.from_dlpack(tz_vis.contiguous())
            cp_gt = cupy.from_dlpack(gt_img.contiguous())
            cp_grad_pos = cupy.from_dlpack(grad_pos)
            cp_grad_scale = cupy.from_dlpack(grad_scale)
            cp_grad_opac = cupy.from_dlpack(grad_opac)
            cp_grad_sh = cupy.from_dlpack(grad_sh)
            cp_grad_2d = cupy.from_dlpack(grad_2d_norm)
            cp_tile_loss = cupy.from_dlpack(tile_loss)
            cp_cam_R_T = cupy.from_dlpack(cam_R_T)

            # Execute Raw CUDA Kernel
            raw_kernel(
                (tiles_x, tiles_y), (16, 16),
                (
                    cp_tile_ranges, cp_items,
                    cp_u, cp_v, cp_inv_xx, cp_inv_xy, cp_inv_yy,
                    cp_opac, cp_r, cp_g, cp_b, cp_g_idx,
                    cp_cam_x, cp_cam_y, cp_cam_z, cp_gt,
                    cp_grad_pos, cp_grad_scale, cp_grad_opac, cp_grad_sh,
                    cp_grad_2d, cp_tile_loss,
                    tiles_x, tiles_y, w, h,
                    fx, fy, cp_cam_R_T
                )
            )

            # 3. Adam Optimizer Step on GPU
            m_pos = beta1 * m_pos + (1 - beta1) * grad_pos
            v_pos = beta2 * v_pos + (1 - beta2) * (grad_pos ** 2)
            pos -= args.lr_pos * (m_pos / (torch.sqrt(v_pos) + eps))

            m_scale = beta1 * m_scale + (1 - beta1) * grad_scale
            v_scale = beta2 * v_scale + (1 - beta2) * (grad_scale ** 2)
            scales_log -= args.lr_scale * (m_scale / (torch.sqrt(v_scale) + eps))
            scales_log = torch.clamp(scales_log, min=-8.0, max=0.5)

            m_opac = beta1 * m_opac + (1 - beta1) * grad_opac
            v_opac = beta2 * v_opac + (1 - beta2) * (grad_opac ** 2)
            opacities_logit -= args.lr_opac * (m_opac / (torch.sqrt(v_opac) + eps))

            m_sh = beta1 * m_sh + (1 - beta1) * grad_sh
            v_sh = beta2 * v_sh + (1 - beta2) * (grad_sh ** 2)
            sh_coeffs -= args.lr_sh * (m_sh / (torch.sqrt(v_sh) + eps))

            accum_grad_2d += grad_2d_norm
            accum_count += 1
            loss_val = float(tile_loss.sum().item())
            recent_losses.append(loss_val)
            if len(recent_losses) > 20:
                recent_losses.pop(0)

        # 4. Adaptive Densification & Pruning
        if step % args.densify_interval == 0 and step <= densify_until:
            mean_grads = accum_grad_2d / max(1, accum_count)
            curr_n = pos.shape[0]

            # High gradients -> Clone or Split
            high_grad = (mean_grads > args.grad_threshold)
            scales_linear = torch.exp(scales_log)
            max_s = torch.max(scales_linear, dim=-1)[0]

            to_clone = high_grad & (max_s < args.scale_threshold)
            to_split = high_grad & (max_s >= args.scale_threshold)

            # Cap total new gaussians
            num_clone = int(to_clone.sum().item())
            num_split = int(to_split.sum().item())

            if curr_n + num_clone + num_split > args.max_gaussians:
                avail = max(0, args.max_gaussians - curr_n)
                if avail < num_clone:
                    idx = torch.nonzero(to_clone).squeeze(-1)[:avail]
                    to_clone = torch.zeros_like(to_clone)
                    to_clone[idx] = True
                    to_split = torch.zeros_like(to_split)
                else:
                    avail_split = avail - num_clone
                    idx_split = torch.nonzero(to_split).squeeze(-1)[:avail_split]
                    to_split = torch.zeros_like(to_split)
                    to_split[idx_split] = True

            new_pos_list = [pos]
            new_scales_list = [scales_log]
            new_quats_list = [quaternions]
            new_opac_list = [opacities_logit]
            new_sh_list = [sh_coeffs]

            if to_clone.any():
                new_pos_list.append(pos[to_clone])
                new_scales_list.append(scales_log[to_clone])
                new_quats_list.append(quaternions[to_clone])
                new_opac_list.append(opacities_logit[to_clone])
                new_sh_list.append(sh_coeffs[to_clone])

            if to_split.any():
                # Split: create 2 Gaussians with smaller scale
                pos_split = pos[to_split]
                noise = torch.randn_like(pos_split) * (max_s[to_split].unsqueeze(-1) * 0.2)
                p1 = pos_split + noise
                p2 = pos_split - noise
                s_split = scales_log[to_split] - 0.223 # log(0.8)

                new_pos_list.extend([p1, p2])
                new_scales_list.extend([s_split, s_split])
                new_quats_list.extend([quaternions[to_split], quaternions[to_split]])
                new_opac_list.extend([opacities_logit[to_split], opacities_logit[to_split]])
                new_sh_list.extend([sh_coeffs[to_split], sh_coeffs[to_split]])

            pos = torch.cat(new_pos_list, dim=0)
            scales_log = torch.cat(new_scales_list, dim=0)
            quaternions = torch.cat(new_quats_list, dim=0)
            opacities_logit = torch.cat(new_opac_list, dim=0)
            sh_coeffs = torch.cat(new_sh_list, dim=0)

            # Pruning
            opac = torch.sigmoid(opacities_logit)
            s_lin = torch.exp(scales_log)
            keep = (opac >= args.min_opacity) & (torch.max(s_lin, dim=-1)[0] <= args.max_scale)

            # Remove original split parents
            if to_split.any():
                keep[:curr_n][to_split] = False

            pos = pos[keep]
            scales_log = scales_log[keep]
            quaternions = quaternions[keep]
            opacities_logit = opacities_logit[keep]
            sh_coeffs = sh_coeffs[keep]

            # Reset optimizer moments and gradient accumulators
            N_new = pos.shape[0]
            m_pos = torch.zeros_like(pos); v_pos = torch.zeros_like(pos)
            m_scale = torch.zeros_like(scales_log); v_scale = torch.zeros_like(scales_log)
            m_opac = torch.zeros_like(opacities_logit); v_opac = torch.zeros_like(opacities_logit)
            m_sh = torch.zeros_like(sh_coeffs); v_sh = torch.zeros_like(sh_coeffs)

            accum_grad_2d = torch.zeros(N_new, device=device)
            accum_count = 0

            avg_loss = float(np.mean(recent_losses)) if recent_losses else 0.0
            print(f"* Step {step}/{args.iterations} | Densified to {N_new:,} Gaussians (Loss: {avg_loss:.4f})", flush=True)

        # Dynamic logging
        log_interval = max(1, min(25, args.iterations // 20))
        if step % log_interval == 0 or step == args.iterations:
            elapsed = time.time() - start_time
            it_per_sec = step / max(0.001, elapsed)
            rem_sec = (args.iterations - step) / max(0.001, it_per_sec)
            eta_str = f"{rem_sec/60:.1f}m" if rem_sec >= 60 else f"{int(rem_sec)}s"
            pct = int((step / args.iterations) * 100)
            avg_loss = float(np.mean(recent_losses)) if recent_losses else 0.0
            print(f"* Step {step}/{args.iterations} [{pct}%] | {it_per_sec:.1f} it/s | ETA: {eta_str} | Loss: {avg_loss:.4f} | Gaussians: {pos.shape[0]:,}", flush=True)

    total_time = round(time.time() - start_time, 1)
    print(f"* GPU Training finished in {total_time}s! Final Gaussians: {pos.shape[0]:,}", flush=True)

    # Save final splat
    final_splat = {
        "num_gaussians": pos.shape[0],
        "positions": pos,
        "scales_log": scales_log,
        "quaternions": quaternions,
        "opacities_logit": opacities_logit,
        "sh_coeffs": sh_coeffs,
        "semantic_ids": torch.zeros(pos.shape[0], dtype=torch.int32, device=device),
        "semantic_scores": torch.zeros(pos.shape[0], dtype=torch.float32, device=device),
    }
    save_splat_binary(final_splat, args.splat_out)

if __name__ == "__main__":
    main()
