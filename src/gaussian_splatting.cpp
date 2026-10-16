// ==============================================================================
// 3D GAUSSIAN SPLATTING ENGINE FOR R (pliman)
// Implementation of Forward Tile Rasterizer, Backward Gradient Engine,
// SAM Semantic Projection, and Adaptive Densification
// ==============================================================================

#include "gaussian_splatting.h"
#include <RcppArmadillo.h>
#include <fstream>
#include <sstream>

using namespace Rcpp;
using namespace pliman_gs;

// Helper to unpack camera from R List
static Camera unpack_camera(const List& cam_list) {
  Camera cam;
  NumericMatrix r_mat = cam_list["R"];
  NumericVector t_vec = cam_list["t"];

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      cam.R(i, j) = r_mat(i, j);
    }
    cam.t(i) = t_vec(i);
  }

  cam.fx = as<double>(cam_list["fx"]);
  cam.fy = as<double>(cam_list["fy"]);
  cam.cx = as<double>(cam_list["cx"]);
  cam.cy = as<double>(cam_list["cy"]);
  cam.width = as<int>(cam_list["width"]);
  cam.height = as<int>(cam_list["height"]);

  // Camera world position = -R^T * t
  cam.pos_world = -cam.R.t() * cam.t;
  return cam;
}

// Key structure for depth sorting: tile_id (high 32 bits) | depth (low 32 bits float as uint)
struct SortItem {
  uint32_t tile_id;
  float depth;
  int splat_idx;

  bool operator<(const SortItem& other) const {
    if (tile_id != other.tile_id) return tile_id < other.tile_id;
    return depth < other.depth; // Front-to-back
  }
};

//' @title Forward 3D Gaussian Splatting Rasterizer
//' @name cpp_gs_render
//' @description
//' High-speed tile-based forward rasterization with EWA covariance projection,
//' Spherical Harmonics color evaluation, depth sorting, and alpha-blending with early stopping.
//' @param camera List containing camera extrinsics (R, t) and intrinsics (fx, fy, cx, cy, width, height).
//' @param gaussians List containing 3D Gaussian parameters (positions, scales_log, quaternions, opacities_logit, sh_coeffs).
//' @param sh_degree Integer. Maximum Spherical Harmonics degree (0, 1, or 2). Defaults to 0.
//' @param bg_r Double. Red background color (0 to 1). Defaults to 0.0.
//' @param bg_g Double. Green background color (0 to 1). Defaults to 0.0.
//' @param bg_b Double. Blue background color (0 to 1). Defaults to 0.0.
//' @param return_depth Logical. If TRUE, returns the 2D depth buffer. Defaults to TRUE.
//' @return A List containing rendered RGB image, alpha transmittance map, and splat cache.
//' @export
// [[Rcpp::export]]
List cpp_gs_render(List camera,
                   List gaussians,
                   int sh_degree = 0,
                   double bg_r = 0.0,
                   double bg_g = 0.0,
                   double bg_b = 0.0,
                   bool return_depth = true) {
  Camera cam = unpack_camera(camera);

  NumericMatrix pos = gaussians["positions"];
  NumericMatrix scale_log = gaussians["scales_log"];
  NumericMatrix quat = gaussians["quaternions"];
  NumericVector opac_logit = gaussians["opacities_logit"];
  NumericMatrix sh_mat = gaussians["sh_coeffs"];

  bool has_sem = gaussians.containsElementNamed("semantic_ids");
  IntegerVector sem_ids = has_sem ? as<IntegerVector>(gaussians["semantic_ids"]) : IntegerVector(0);
  NumericVector sem_scores = has_sem && gaussians.containsElementNamed("semantic_scores") ?
                             as<NumericVector>(gaussians["semantic_scores"]) : NumericVector(0);

  int num_gauss = pos.nrow();
  int num_sh_coeffs = sh_mat.ncol();

  // 1. EWA Projection: Transform 3D Gaussians to 2D view space
  std::vector<Splat2D> visible_splats;
  visible_splats.reserve(num_gauss / 2);

  for (int i = 0; i < num_gauss; ++i) {
    arma::vec3 p_world = {pos(i, 0), pos(i, 1), pos(i, 2)};
    arma::vec3 p_cam = cam.R * p_world + cam.t;

    double tz = p_cam(2);
    if (tz <= 0.1) continue; // Near plane culling

    double tx = p_cam(0);
    double ty = p_cam(1);

    // Perspective projection of center
    double u = cam.fx * (tx / tz) + cam.cx;
    double v = cam.fy * (ty / tz) + cam.cy;

    // Jacobian of perspective projection J
    double tz2 = tz * tz;
    arma::mat J(2, 3, arma::fill::zeros);
    J(0, 0) = cam.fx / tz;
    J(0, 1) = 0.0;
    J(0, 2) = -(cam.fx * tx) / tz2;

    J(1, 0) = 0.0;
    J(1, 1) = cam.fy / tz;
    J(1, 2) = -(cam.fy * ty) / tz2;

    // 3D Covariance in world coordinates
    double sx = std::exp(scale_log(i, 0));
    double sy = std::exp(scale_log(i, 1));
    double sz = std::exp(scale_log(i, 2));

    arma::mat33 R_q = quat_to_rot(quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3));
    arma::mat33 cov3d_world = compute_cov3d(sx, sy, sz, R_q);

    // 3D Covariance in camera coordinates V = R_cam * cov3d * R_cam^T
    arma::mat33 cov3d_cam = cam.R * cov3d_world * cam.R.t();

    // 2D Covariance Sigma_2D = J * cov3d_cam * J^T
    arma::mat cov2d = J * cov3d_cam * J.t();

    // Anti-aliasing low-pass filter (add 0.3 to diagonal)
    cov2d(0, 0) += 0.3;
    cov2d(1, 1) += 0.3;

    double det = cov2d(0, 0) * cov2d(1, 1) - cov2d(0, 1) * cov2d(1, 0);
    if (det < 1e-7) continue;

    double inv_det = 1.0 / det;
    double inv_cov_xx = cov2d(1, 1) * inv_det;
    double inv_cov_xy = -cov2d(0, 1) * inv_det;
    double inv_cov_yy = cov2d(0, 0) * inv_det;

    // 3-sigma bounding radius in pixels
    double mid = 0.5 * (cov2d(0, 0) + cov2d(1, 1));
    double term = 0.25 * (cov2d(0, 0) - cov2d(1, 1)) * (cov2d(0, 0) - cov2d(1, 1)) + cov2d(0, 1) * cov2d(0, 1);
    double lambda_max = mid + std::sqrt(std::max(0.01, term));
    int radius = std::max(1, (int)std::ceil(3.0 * std::sqrt(lambda_max)));

    // Frustum culling in 2D viewport
    if (u + radius < 0 || u - radius >= cam.width ||
        v + radius < 0 || v - radius >= cam.height) {
      continue;
    }

    double opac = sigmoid(opac_logit[i]);
    if (opac < 1.0 / 255.0) continue;

    // Evaluate Spherical Harmonics
    std::vector<double> sh_vec(num_sh_coeffs);
    for (int k = 0; k < num_sh_coeffs; ++k) sh_vec[k] = sh_mat(i, k);

    arma::vec3 ray_dir = p_world - cam.pos_world;
    double cr = 0.5, cg = 0.5, cb = 0.5;
    eval_sh(sh_degree, sh_vec, ray_dir(0), ray_dir(1), ray_dir(2), cr, cg, cb);

    Splat2D splat;
    splat.gaussian_idx = i;
    splat.u = u;
    splat.v = v;
    splat.depth = tz;
    splat.inv_cov_xx = inv_cov_xx;
    splat.inv_cov_xy = inv_cov_xy;
    splat.inv_cov_yy = inv_cov_yy;
    splat.radius = radius;
    splat.opacity = opac;
    splat.r = cr;
    splat.g = cg;
    splat.b = cb;
    splat.semantic_id = (has_sem && i < sem_ids.size()) ? sem_ids[i] : 0;
    splat.semantic_score = (has_sem && i < sem_scores.size()) ? sem_scores[i] : 0.0;

    visible_splats.push_back(splat);
  }

  int num_visible = (int)visible_splats.size();

  // 2. Tile Assignment and Radix Depth Sorting
  const int TILE_SIZE = 16;
  int tiles_x = (cam.width + TILE_SIZE - 1) / TILE_SIZE;
  int tiles_y = (cam.height + TILE_SIZE - 1) / TILE_SIZE;
  int total_tiles = tiles_x * tiles_y;

  std::vector<SortItem> sort_items;
  sort_items.reserve(num_visible * 4);

  for (int s_idx = 0; s_idx < num_visible; ++s_idx) {
    const Splat2D& s = visible_splats[s_idx];
    int min_tx = std::max(0, (int)std::floor((s.u - s.radius) / TILE_SIZE));
    int max_tx = std::min(tiles_x - 1, (int)std::floor((s.u + s.radius) / TILE_SIZE));
    int min_ty = std::max(0, (int)std::floor((s.v - s.radius) / TILE_SIZE));
    int max_ty = std::min(tiles_y - 1, (int)std::floor((s.v + s.radius) / TILE_SIZE));

    for (int ty = min_ty; ty <= max_ty; ++ty) {
      for (int tx = min_tx; tx <= max_tx; ++tx) {
        uint32_t t_id = (uint32_t)(ty * tiles_x + tx);
        sort_items.push_back({t_id, (float)s.depth, s_idx});
      }
    }
  }

  // Parallel sort by (tile_id, depth)
  std::sort(sort_items.begin(), sort_items.end());

  // Find tile start/end ranges
  std::vector<int> tile_ranges(total_tiles + 1, (int)sort_items.size());
  uint32_t curr_tile = 0;
  tile_ranges[0] = 0;
  for (size_t k = 0; k < sort_items.size(); ++k) {
    uint32_t t_id = sort_items[k].tile_id;
    while (curr_tile < t_id) {
      curr_tile++;
      tile_ranges[curr_tile] = (int)k;
    }
  }

  // 3. Tile-Based Parallel Rasterization (OpenMP)
  int w = cam.width;
  int h = cam.height;

  NumericVector out_r(w * h, bg_r);
  NumericVector out_g(w * h, bg_g);
  NumericVector out_b(w * h, bg_b);
  NumericVector out_alpha(w * h, 0.0);
  NumericVector out_depth(w * h, 0.0);

#pragma omp parallel for schedule(dynamic, 4)
  for (int t_idx = 0; t_idx < total_tiles; ++t_idx) {
    int start_item = tile_ranges[t_idx];
    int end_item = (t_idx + 1 < (int)tile_ranges.size()) ? tile_ranges[t_idx + 1] : (int)sort_items.size();
    if (start_item >= end_item) continue;

    int tx = t_idx % tiles_x;
    int ty = t_idx / tiles_x;

    int px_start = tx * TILE_SIZE;
    int px_end = std::min(w, px_start + TILE_SIZE);
    int py_start = ty * TILE_SIZE;
    int py_end = std::min(h, py_start + TILE_SIZE);

    for (int py = py_start; py < py_end; ++py) {
      for (int px = px_start; px < px_end; ++px) {
        double T = 1.0; // Transmittance
        double acc_r = 0.0;
        double acc_g = 0.0;
        double acc_b = 0.0;
        double acc_depth = 0.0;
        double first_depth = 0.0;

        for (int k = start_item; k < end_item; ++k) {
          int s_idx = sort_items[k].splat_idx;
          const Splat2D& splat = visible_splats[s_idx];

          double dx = (double)px - splat.u;
          double dy = (double)py - splat.v;

          double power = -0.5 * (dx * dx * splat.inv_cov_xx +
                                 2.0 * dx * dy * splat.inv_cov_xy +
                                 dy * dy * splat.inv_cov_yy);

          if (power > 0.0) continue;

          double alpha = std::min(0.99, splat.opacity * std::exp(power));
          if (alpha < 1.0 / 255.0) continue;

          double weight = alpha * T;
          acc_r += splat.r * weight;
          acc_g += splat.g * weight;
          acc_b += splat.b * weight;
          acc_depth += splat.depth * weight;

          if (first_depth == 0.0) first_depth = splat.depth;

          T *= (1.0 - alpha);
          if (T < 0.0001) break; // Early ray termination
        }

        int pix_idx = px + py * w;
        double final_alpha = 1.0 - T;
        out_r[pix_idx] = acc_r + T * bg_r;
        out_g[pix_idx] = acc_g + T * bg_g;
        out_b[pix_idx] = acc_b + T * bg_b;
        out_alpha[pix_idx] = final_alpha;
        out_depth[pix_idx] = (final_alpha > 0.1) ? (acc_depth / final_alpha) : 0.0;
      }
    }
  }

  // Format image tensor (w x h x 3) matching pliman conventions
  NumericVector img_tensor(w * h * 3);
  for (int p = 0; p < w * h; ++p) {
    img_tensor[p] = clamp(out_r[p], 0.0, 1.0);
    img_tensor[p + w * h] = clamp(out_g[p], 0.0, 1.0);
    img_tensor[p + 2 * w * h] = clamp(out_b[p], 0.0, 1.0);
  }
  img_tensor.attr("dim") = IntegerVector::create(w, h, 3);

  NumericMatrix mat_alpha(w, h);
  NumericMatrix mat_depth(w, h);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      mat_alpha(x, y) = out_alpha[x + y * w];
      mat_depth(x, y) = out_depth[x + y * w];
    }
  }

  // Pack splat projection cache for backward pass and semantic queries
  NumericMatrix splat_cache(num_visible, 9);
  for (int i = 0; i < num_visible; ++i) {
    splat_cache(i, 0) = visible_splats[i].gaussian_idx;
    splat_cache(i, 1) = visible_splats[i].u;
    splat_cache(i, 2) = visible_splats[i].v;
    splat_cache(i, 3) = visible_splats[i].depth;
    splat_cache(i, 4) = visible_splats[i].radius;
    splat_cache(i, 5) = visible_splats[i].opacity;
    splat_cache(i, 6) = visible_splats[i].r;
    splat_cache(i, 7) = visible_splats[i].g;
    splat_cache(i, 8) = visible_splats[i].b;
  }

  return List::create(
    Named("image") = img_tensor,
    Named("alpha") = mat_alpha,
    Named("depth") = mat_depth,
    Named("splats_count") = num_visible,
    Named("cache") = splat_cache
  );
}

//' @title Propagate 2D SAM Masks to 3D Gaussian Splats
//' @name cpp_gs_segment_sam
//' @description
//' Projects 2D instance segmentation masks (e.g. from SAM 2.1 or SAM 3.1)
//' onto 3D Gaussians via transmittance-weighted ray-splat correspondence.
//' @param camera List camera specification.
//' @param gaussians List of 3D Gaussians.
//' @param mask 2D IntegerMatrix of binary mask (0 = background, 1 = target leaf/organ).
//' @param class_id Integer. Semantic class ID to assign.
//' @param score_threshold Double. Minimum confidence score to update. Defaults to 0.3.
//' @return A List with updated semantic IDs, scores, and per-splat weights.
//' @export
// [[Rcpp::export]]
List cpp_gs_segment_sam(List camera,
                        List gaussians,
                        IntegerMatrix mask,
                        int class_id = 1,
                        double score_threshold = 0.3) {
  Camera cam = unpack_camera(camera);

  NumericMatrix pos = gaussians["positions"];
  NumericMatrix scale_log = gaussians["scales_log"];
  NumericMatrix quat = gaussians["quaternions"];
  NumericVector opac_logit = gaussians["opacities_logit"];

  int num_gauss = pos.nrow();
  int mask_w = mask.nrow();
  int mask_h = mask.ncol();

  IntegerVector curr_sem(num_gauss, 0);
  NumericVector curr_scores(num_gauss, 0.0);
  NumericVector hit_weights(num_gauss, 0.0);

  if (gaussians.containsElementNamed("semantic_ids")) {
    curr_sem = clone(as<IntegerVector>(gaussians["semantic_ids"]));
  }
  if (gaussians.containsElementNamed("semantic_scores")) {
    curr_scores = clone(as<NumericVector>(gaussians["semantic_scores"]));
  }

#pragma omp parallel for schedule(dynamic, 64)
  for (int i = 0; i < num_gauss; ++i) {
    arma::vec3 p_world = {pos(i, 0), pos(i, 1), pos(i, 2)};
    arma::vec3 p_cam = cam.R * p_world + cam.t;
    double tz = p_cam(2);
    if (tz <= 0.1) continue;

    double u = cam.fx * (p_cam(0) / tz) + cam.cx;
    double v = cam.fy * (p_cam(1) / tz) + cam.cy;

    int iu = (int)std::round(u);
    int iv = (int)std::round(v);

    bool inside = false;
    if (iu >= 0 && iu < mask_w && iv >= 0 && iv < mask_h) {
      inside = (mask(iu, iv) > 0);
    } else if (iv >= 0 && iv < mask_w && iu >= 0 && iu < mask_h) {
      inside = (mask(iv, iu) > 0);
    }

    if (inside) {
      double opac = sigmoid(opac_logit[i]);
      // Confidence score based on opacity and depth attenuation
      double dist_factor = 1.0 / (1.0 + 0.1 * tz);
      double weight = opac * dist_factor;
      hit_weights[i] = weight;

      if (weight >= score_threshold) {
        curr_sem[i] = class_id;
        curr_scores[i] = std::max(curr_scores[i], weight);
      }
    }
  }

  return List::create(
    Named("semantic_ids") = curr_sem,
    Named("semantic_scores") = curr_scores,
    Named("hit_weights") = hit_weights
  );
}

//' @title 3D Plant Biometrics from Gaussian Splats
//' @name cpp_gs_measure
//' @description
//' Computes real 3D phenotypic metrics (surface area, volume, leaf angles, dimensions)
//' from 3D Gaussian Splats, with optional filtering by semantic class (e.g. leaf only).
//' @param positions N x 3 NumericMatrix of Gaussian positions.
//' @param scales_log N x 3 NumericMatrix of log-scales.
//' @param quaternions N x 4 NumericMatrix of rotation quaternions.
//' @param semantic_ids IntegerVector of semantic class IDs.
//' @param target_class Integer. Target class ID to measure (-1 = measure all). Defaults to -1.
//' @return A List containing 3D surface area, bounding dimensions, and insertion angles.
//' @export
// [[Rcpp::export]]
List cpp_gs_measure(NumericMatrix positions,
                    NumericMatrix scales_log,
                    NumericMatrix quaternions,
                    IntegerVector semantic_ids,
                    int target_class = -1) {
  int n = positions.nrow();
  if (n == 0) {
    return List::create(
      Named("surface_area") = 0.0,
      Named("num_splats") = 0,
      Named("height") = 0.0,
      Named("width") = 0.0,
      Named("depth") = 0.0
    );
  }

  double min_x = 1e9, max_x = -1e9;
  double min_y = 1e9, max_y = -1e9;
  double min_z = 1e9, max_z = -1e9;

  double total_area = 0.0;
  int count = 0;
  std::vector<double> leaf_angles;

  for (int i = 0; i < n; ++i) {
    if (target_class >= 0 && semantic_ids.size() > i && semantic_ids[i] != target_class) {
      continue;
    }

    double x = positions(i, 0);
    double y = positions(i, 1);
    double z = positions(i, 2);

    min_x = std::min(min_x, x); max_x = std::max(max_x, x);
    min_y = std::min(min_y, y); max_y = std::max(max_y, y);
    min_z = std::min(min_z, z); max_z = std::max(max_z, z);

    double sx = std::exp(scales_log(i, 0));
    double sy = std::exp(scales_log(i, 1));

    // Splat cross-sectional ellipse area = pi * sx * sy
    total_area += M_PI * sx * sy;

    // Normal vector from quaternion (z-axis of rotated ellipsoid)
    arma::mat33 R = quat_to_rot(quaternions(i, 0), quaternions(i, 1), quaternions(i, 2), quaternions(i, 3));
    arma::vec3 normal = R.col(2); // Local z-axis

    // Angle relative to vertical (upward z: [0, 0, 1]) in degrees
    double dot_up = std::abs(normal(2));
    double angle_deg = std::acos(clamp(dot_up, 0.0, 1.0)) * 180.0 / M_PI;
    leaf_angles.push_back(angle_deg);

    count++;
  }

  double width  = (count > 0) ? std::max(0.0, max_x - min_x) : 0.0;
  double depth  = (count > 0) ? std::max(0.0, max_y - min_y) : 0.0;
  double height = (count > 0) ? std::max(0.0, max_z - min_z) : 0.0;

  return List::create(
    Named("surface_area") = total_area,
    Named("num_splats") = count,
    Named("width") = width,
    Named("depth") = depth,
    Named("height") = height,
    Named("bbox_min") = NumericVector::create(min_x, min_y, min_z),
    Named("bbox_max") = NumericVector::create(max_x, max_y, max_z),
    Named("leaf_angles") = wrap(leaf_angles)
  );
}

//' @title Adaptive Densification and Pruning for 3D Gaussian Splats
//' @name cpp_gs_densify
//' @description
//' Implements split, clone, and prune operations based on view-space positional
//' gradients and scale thresholds to densify under-reconstructed geometry and prune artifacts.
//' @param gaussians List of 3D Gaussians.
//' @param grad_2d_norm NumericVector of accumulated 2D view-space positional gradient norms.
//' @param grad_threshold Double. Gradient threshold to trigger densification. Defaults to 0.0002.
//' @param scale_threshold Double. Scale threshold distinguishing split from clone. Defaults to 0.01.
//' @param min_opacity Double. Minimum opacity threshold for pruning. Defaults to 0.05.
//' @param max_scale Double. Maximum world scale for pruning. Defaults to 0.5.
//' @return An updated List of 3D Gaussians.
//' @export
// [[Rcpp::export]]
List cpp_gs_densify(List gaussians,
                    NumericVector grad_2d_norm,
                    double grad_threshold = 0.0002,
                    double scale_threshold = 0.01,
                    double min_opacity = 0.05,
                    double max_scale = 0.5,
                    int max_gaussians = 200000) {
  NumericMatrix pos = gaussians["positions"];
  NumericMatrix scale_log = gaussians["scales_log"];
  NumericMatrix quat = gaussians["quaternions"];
  NumericVector opac_logit = gaussians["opacities_logit"];
  NumericMatrix sh_mat = gaussians["sh_coeffs"];

  bool has_sem = gaussians.containsElementNamed("semantic_ids");
  IntegerVector sem_ids = has_sem ? as<IntegerVector>(gaussians["semantic_ids"]) : IntegerVector(0);
  NumericVector sem_scores = has_sem && gaussians.containsElementNamed("semantic_scores") ?
                             as<NumericVector>(gaussians["semantic_scores"]) : NumericVector(0);

  int n = pos.nrow();
  int num_sh = sh_mat.ncol();

  std::vector<Gaussian> kept_gaussians;
  kept_gaussians.reserve(std::min(n * 2, max_gaussians + 1000));

  for (int i = 0; i < n; ++i) {
    double opac = sigmoid(opac_logit[i]);
    double sx = std::exp(scale_log(i, 0));
    double sy = std::exp(scale_log(i, 1));
    double sz = std::exp(scale_log(i, 2));
    double max_s = std::max({sx, sy, sz});

    // Prune transparent or excessively large Gaussians
    if (opac < min_opacity || max_s > max_scale) {
      continue;
    }

    std::vector<double> sh_vec(num_sh);
    for (int k = 0; k < num_sh; ++k) sh_vec[k] = sh_mat(i, k);
    int s_id = (has_sem && i < sem_ids.size()) ? sem_ids[i] : 0;
    double s_score = (has_sem && i < sem_scores.size()) ? sem_scores[i] : 0.0;

    double grad = (i < grad_2d_norm.size()) ? grad_2d_norm[i] : 0.0;

    if (grad > grad_threshold && (int)kept_gaussians.size() < max_gaussians) {
      if (max_s > scale_threshold) {
        // SPLIT: divide into 2 smaller Gaussians (scale / 1.6)
        double new_sx = std::log(sx / 1.6);
        double new_sy = std::log(sy / 1.6);
        double new_sz = std::log(sz / 1.6);

        // Perturb position along major axis
        arma::mat33 R = quat_to_rot(quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3));
        arma::vec3 offset = R.col(0) * (sx * 0.3);

        Gaussian g1 = {pos(i, 0) + offset(0), pos(i, 1) + offset(1), pos(i, 2) + offset(2),
                       new_sx, new_sy, new_sz, quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3),
                       opac_logit[i], sh_vec, s_id, s_score};
        Gaussian g2 = {pos(i, 0) - offset(0), pos(i, 1) - offset(1), pos(i, 2) - offset(2),
                       new_sx, new_sy, new_sz, quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3),
                       opac_logit[i], sh_vec, s_id, s_score};
        kept_gaussians.push_back(g1);
        kept_gaussians.push_back(g2);
      } else {
        // CLONE: duplicate with slight perturbation
        Gaussian orig = {pos(i, 0), pos(i, 1), pos(i, 2),
                         scale_log(i, 0), scale_log(i, 1), scale_log(i, 2),
                         quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3),
                         opac_logit[i], sh_vec, s_id, s_score};
        Gaussian clone_g = {pos(i, 0) + 0.001 * (sx), pos(i, 1) + 0.001 * (sy), pos(i, 2) + 0.001 * (sz),
                            scale_log(i, 0), scale_log(i, 1), scale_log(i, 2),
                            quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3),
                            opac_logit[i], sh_vec, s_id, s_score};
        kept_gaussians.push_back(orig);
        kept_gaussians.push_back(clone_g);
      }
    } else {
      // Keep unchanged
      Gaussian g = {pos(i, 0), pos(i, 1), pos(i, 2),
                    scale_log(i, 0), scale_log(i, 1), scale_log(i, 2),
                    quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3),
                    opac_logit[i], sh_vec, s_id, s_score};
      kept_gaussians.push_back(g);
    }
  }

  int new_n = (int)kept_gaussians.size();
  NumericMatrix new_pos(new_n, 3);
  NumericMatrix new_scale(new_n, 3);
  NumericMatrix new_quat(new_n, 4);
  NumericVector new_opac(new_n);
  NumericMatrix new_sh(new_n, num_sh);
  IntegerVector new_sem(new_n);
  NumericVector new_scores(new_n);

  for (int i = 0; i < new_n; ++i) {
    const Gaussian& g = kept_gaussians[i];
    new_pos(i, 0) = g.x; new_pos(i, 1) = g.y; new_pos(i, 2) = g.z;
    new_scale(i, 0) = g.scale_log_x; new_scale(i, 1) = g.scale_log_y; new_scale(i, 2) = g.scale_log_z;
    new_quat(i, 0) = g.qw; new_quat(i, 1) = g.qx; new_quat(i, 2) = g.qy; new_quat(i, 3) = g.qz;
    new_opac[i] = g.opacity_logit;
    for (int k = 0; k < num_sh; ++k) new_sh(i, k) = g.sh[k];
    new_sem[i] = g.semantic_id;
    new_scores[i] = g.semantic_score;
  }

  return List::create(
    Named("positions") = new_pos,
    Named("scales_log") = new_scale,
    Named("quaternions") = new_quat,
    Named("opacities_logit") = new_opac,
    Named("sh_coeffs") = new_sh,
    Named("semantic_ids") = new_sem,
    Named("semantic_scores") = new_scores,
    Named("num_gaussians") = new_n
  );
}

//' @title High-Speed Binary and ASCII PLY Point Cloud Parser
//' @name cpp_gs_read_ply
//' @description
//' Reads 3D point coordinates (X, Y, Z) and RGB vertex colors from PLY files
//' exported by Agisoft Metashape, COLMAP, CloudCompare, or LiDAR scanners.
//' @param file_path String path to the .ply file.
//' @return A List containing 'positions' (N x 3 NumericMatrix), 'colors' (N x 3 NumericMatrix), and 'num_vertices'.
//' @export
// [[Rcpp::export]]
List cpp_gs_read_ply(std::string file_path) {
  std::ifstream infile(file_path, std::ios::binary);
  if (!infile.is_open()) {
    stop("Cannot open PLY file: " + file_path);
  }

  std::string line;
  std::getline(infile, line);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line != "ply") {
    stop("Invalid PLY format: header does not begin with 'ply'.");
  }

  bool is_binary = false;
  int num_vertices = 0;
  int num_faces = 0;
  bool in_vertex = false;
  std::string face_count_type = "uchar";
  std::string face_index_type = "int";

  struct PlyProp {
    std::string name;
    std::string type;
    int size_bytes;
    int offset;
  };

  std::vector<PlyProp> properties;
  int current_offset = 0;

  while (std::getline(infile, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line == "end_header") break;

    std::istringstream iss(line);
    std::string token;
    iss >> token;

    if (token == "format") {
      std::string fmt_type;
      iss >> fmt_type;
      if (fmt_type == "binary_little_endian") {
        is_binary = true;
      } else if (fmt_type == "ascii") {
        is_binary = false;
      } else {
        stop("Unsupported PLY format: " + fmt_type);
      }
    } else if (token == "element") {
      std::string elem_type;
      int count;
      iss >> elem_type >> count;
      if (elem_type == "vertex") {
        in_vertex = true;
        num_vertices = count;
      } else if (elem_type == "face") {
        in_vertex = false;
        num_faces = count;
      } else {
        in_vertex = false;
      }
    } else if (token == "property") {
      if (in_vertex) {
        std::string prop_type, prop_name;
        iss >> prop_type >> prop_name;

        int s = 4;
        if (prop_type == "float" || prop_type == "float32" || prop_type == "int" || prop_type == "int32" || prop_type == "uint" || prop_type == "uint32") {
          s = 4;
        } else if (prop_type == "double" || prop_type == "float64") {
          s = 8;
        } else if (prop_type == "uchar" || prop_type == "uint8" || prop_type == "char" || prop_type == "int8") {
          s = 1;
        } else if (prop_type == "short" || prop_type == "int16" || prop_type == "ushort" || prop_type == "uint16") {
          s = 2;
        }

        PlyProp p = {prop_name, prop_type, s, current_offset};
        properties.push_back(p);
        current_offset += s;
      } else {
        std::string p_kind, p_cnt, p_idx, p_name;
        iss >> p_kind >> p_cnt >> p_idx >> p_name;
        if (p_kind == "list") {
          face_count_type = p_cnt;
          face_index_type = p_idx;
        }
      }
    }
  }

  int vertex_stride = current_offset;
  if (num_vertices <= 0 || vertex_stride == 0) {
    stop("PLY file contains 0 vertices or empty vertex properties.");
  }

  int off_x = -1, off_y = -1, off_z = -1;
  int off_r = -1, off_g = -1, off_b = -1;
  std::string type_x = "float", type_r = "uchar";

  for (const auto& p : properties) {
    std::string n = p.name;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "x") { off_x = p.offset; type_x = p.type; }
    else if (n == "y") { off_y = p.offset; }
    else if (n == "z") { off_z = p.offset; }
    else if (n == "red" || n == "diffuse_red" || n == "r") { off_r = p.offset; type_r = p.type; }
    else if (n == "green" || n == "diffuse_green" || n == "g") { off_g = p.offset; }
    else if (n == "blue" || n == "diffuse_blue" || n == "b") { off_b = p.offset; }
  }

  if (off_x == -1 || off_y == -1 || off_z == -1) {
    stop("PLY file is missing required x, y, or z properties.");
  }

  NumericMatrix pos(num_vertices, 3);
  NumericMatrix col(num_vertices, 3);
  bool has_color = (off_r != -1 && off_g != -1 && off_b != -1);

  if (is_binary) {
    std::vector<char> buffer(static_cast<size_t>(num_vertices) * vertex_stride);
    infile.read(buffer.data(), buffer.size());

#pragma omp parallel for schedule(static, 2048)
    for (int i = 0; i < num_vertices; ++i) {
      const char* ptr = buffer.data() + static_cast<size_t>(i) * vertex_stride;

      if (type_x == "double" || type_x == "float64") {
        pos(i, 0) = *reinterpret_cast<const double*>(ptr + off_x);
        pos(i, 1) = *reinterpret_cast<const double*>(ptr + off_y);
        pos(i, 2) = *reinterpret_cast<const double*>(ptr + off_z);
      } else {
        pos(i, 0) = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_x));
        pos(i, 1) = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_y));
        pos(i, 2) = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_z));
      }

      if (has_color) {
        if (type_r == "uchar" || type_r == "uint8") {
          col(i, 0) = static_cast<double>(*reinterpret_cast<const unsigned char*>(ptr + off_r)) / 255.0;
          col(i, 1) = static_cast<double>(*reinterpret_cast<const unsigned char*>(ptr + off_g)) / 255.0;
          col(i, 2) = static_cast<double>(*reinterpret_cast<const unsigned char*>(ptr + off_b)) / 255.0;
        } else if (type_r == "float" || type_r == "float32") {
          double r_val = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_r));
          double g_val = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_g));
          double b_val = static_cast<double>(*reinterpret_cast<const float*>(ptr + off_b));
          if (r_val > 1.0 || g_val > 1.0 || b_val > 1.0) {
            r_val /= 255.0; g_val /= 255.0; b_val /= 255.0;
          }
          col(i, 0) = r_val; col(i, 1) = g_val; col(i, 2) = b_val;
        } else {
          col(i, 0) = 0.22; col(i, 1) = 0.72; col(i, 2) = 0.28;
        }
      } else {
        if (num_faces > 0) {
          col(i, 0) = 0.82; col(i, 1) = 0.80; col(i, 2) = 0.78;
        } else {
          col(i, 0) = 0.22; col(i, 1) = 0.72; col(i, 2) = 0.28;
        }
      }
    }
  } else {
    for (int i = 0; i < num_vertices; ++i) {
      if (!std::getline(infile, line)) break;
      std::istringstream iss(line);
      std::vector<double> vals;
      double v;
      while (iss >> v) vals.push_back(v);
      if (vals.size() >= 3) {
        pos(i, 0) = vals[0];
        pos(i, 1) = vals[1];
        pos(i, 2) = vals[2];
        if (vals.size() >= 6 && has_color) {
          col(i, 0) = vals[3] > 1.0 ? vals[3] / 255.0 : vals[3];
          col(i, 1) = vals[4] > 1.0 ? vals[4] / 255.0 : vals[4];
          col(i, 2) = vals[5] > 1.0 ? vals[5] / 255.0 : vals[5];
        } else {
          if (num_faces > 0) {
            col(i, 0) = 0.82; col(i, 1) = 0.80; col(i, 2) = 0.78;
          } else {
            col(i, 0) = 0.22; col(i, 1) = 0.72; col(i, 2) = 0.28;
          }
        }
      }
    }
  }

  // Parse element face if present (3D polygonal mesh)
  std::vector<int> face_indices;
  if (num_faces > 0) {
    face_indices.reserve(static_cast<size_t>(num_faces) * 3);
    if (is_binary) {
      for (int f = 0; f < num_faces; ++f) {
        unsigned char n_v = 0;
        infile.read(reinterpret_cast<char*>(&n_v), 1);
        if (!infile) break;

        std::vector<int> v_idx(n_v);
        if (face_index_type == "int" || face_index_type == "int32" || face_index_type == "uint" || face_index_type == "uint32") {
          infile.read(reinterpret_cast<char*>(v_idx.data()), n_v * 4);
        } else if (face_index_type == "short" || face_index_type == "uint16" || face_index_type == "int16") {
          for (int k = 0; k < n_v; ++k) {
            unsigned short s_val = 0;
            infile.read(reinterpret_cast<char*>(&s_val), 2);
            v_idx[k] = s_val;
          }
        } else {
          infile.read(reinterpret_cast<char*>(v_idx.data()), n_v * 4);
        }

        if (n_v == 3) {
          face_indices.push_back(v_idx[0]);
          face_indices.push_back(v_idx[1]);
          face_indices.push_back(v_idx[2]);
        } else if (n_v == 4) {
          // Triangulate quad into two triangles
          face_indices.push_back(v_idx[0]);
          face_indices.push_back(v_idx[1]);
          face_indices.push_back(v_idx[2]);
          face_indices.push_back(v_idx[0]);
          face_indices.push_back(v_idx[2]);
          face_indices.push_back(v_idx[3]);
        } else if (n_v > 4) {
          for (int k = 1; k < n_v - 1; ++k) {
            face_indices.push_back(v_idx[0]);
            face_indices.push_back(v_idx[k]);
            face_indices.push_back(v_idx[k + 1]);
          }
        }
      }
    } else {
      for (int f = 0; f < num_faces; ++f) {
        if (!std::getline(infile, line)) break;
        std::istringstream iss(line);
        int n_v;
        if (iss >> n_v && n_v >= 3) {
          std::vector<int> v_idx(n_v);
          for (int k = 0; k < n_v; ++k) iss >> v_idx[k];
          if (n_v == 3) {
            face_indices.push_back(v_idx[0]);
            face_indices.push_back(v_idx[1]);
            face_indices.push_back(v_idx[2]);
          } else if (n_v == 4) {
            face_indices.push_back(v_idx[0]);
            face_indices.push_back(v_idx[1]);
            face_indices.push_back(v_idx[2]);
            face_indices.push_back(v_idx[0]);
            face_indices.push_back(v_idx[2]);
            face_indices.push_back(v_idx[3]);
          }
        }
      }
    }
  }

  IntegerVector faces_r = wrap(face_indices);
  return List::create(
    Named("positions") = pos,
    Named("colors") = col,
    Named("num_vertices") = num_vertices,
    Named("faces") = faces_r,
    Named("num_faces") = static_cast<int>(face_indices.size() / 3),
    Named("has_colors") = has_color
  );
}

//' @title Differentiable 3D Gaussian Splatting Training Step
//' @name cpp_gs_train_step
//' @description
//' Executes forward rasterization, photometric loss calculation (L1),
//' analytical gradient backpropagation, and Adam parameter updates.
//' @param camera List of camera parameters.
//' @param gaussians List of current 3D Gaussian parameters.
//' @param target_img NumericVector representing the W x H x 3 target image.
//' @param adam_m List of first Adam momentum moments (pos, scale, opac, sh).
//' @param adam_v List of second Adam momentum moments (pos, scale, opac, sh).
//' @param iter Integer. Current iteration number.
//' @param lr_pos Double. Learning rate for 3D positions. Defaults to 0.00016.
//' @param lr_scale Double. Learning rate for log-scales. Defaults to 0.005.
//' @param lr_opac Double. Learning rate for opacity logit. Defaults to 0.05.
//' @param lr_sh Double. Learning rate for colors / spherical harmonics. Defaults to 0.0025.
//' @return A List with updated gaussians, adam_m, adam_v, loss, and grad_2d_norm.
//' @export
// [[Rcpp::export]]
List cpp_gs_train_step(List camera,
                       List gaussians,
                       NumericVector target_img,
                       List adam_m,
                       List adam_v,
                       int iter = 1,
                       double lr_pos = 0.00016,
                       double lr_scale = 0.005,
                       double lr_opac = 0.05,
                       double lr_sh = 0.0025) {
  Camera cam = unpack_camera(camera);

  NumericMatrix pos = clone(as<NumericMatrix>(gaussians["positions"]));
  NumericMatrix scale_log = clone(as<NumericMatrix>(gaussians["scales_log"]));
  NumericMatrix quat = clone(as<NumericMatrix>(gaussians["quaternions"]));
  NumericVector opac_logit = clone(as<NumericVector>(gaussians["opacities_logit"]));
  NumericMatrix sh_mat = clone(as<NumericMatrix>(gaussians["sh_coeffs"]));

  NumericMatrix m_pos = as<NumericMatrix>(adam_m["pos"]);
  NumericMatrix v_pos = as<NumericMatrix>(adam_v["pos"]);
  NumericMatrix m_scale = as<NumericMatrix>(adam_m["scale"]);
  NumericMatrix v_scale = as<NumericMatrix>(adam_v["scale"]);
  NumericVector m_opac = as<NumericVector>(adam_m["opac"]);
  NumericVector v_opac = as<NumericVector>(adam_v["opac"]);
  NumericMatrix m_sh = as<NumericMatrix>(adam_m["sh"]);
  NumericMatrix v_sh = as<NumericMatrix>(adam_v["sh"]);

  int num_gauss = pos.nrow();
  int w = cam.width;
  int h = cam.height;

  // 1. EWA Projection
  std::vector<Splat2D> visible_splats;
  std::vector<arma::vec3> p_cams;
  visible_splats.reserve(num_gauss / 2);
  p_cams.reserve(num_gauss / 2);

  for (int i = 0; i < num_gauss; ++i) {
    arma::vec3 p_world = {pos(i, 0), pos(i, 1), pos(i, 2)};
    arma::vec3 p_cam = cam.R * p_world + cam.t;
    double tz = p_cam(2);
    if (tz <= 0.1) continue;

    double tx = p_cam(0);
    double ty = p_cam(1);

    double u = cam.fx * (tx / tz) + cam.cx;
    double v = cam.fy * (ty / tz) + cam.cy;

    double tz2 = tz * tz;
    arma::mat J(2, 3, arma::fill::zeros);
    J(0, 0) = cam.fx / tz;
    J(0, 2) = -(cam.fx * tx) / tz2;
    J(1, 1) = cam.fy / tz;
    J(1, 2) = -(cam.fy * ty) / tz2;

    double sx = std::exp(scale_log(i, 0));
    double sy = std::exp(scale_log(i, 1));
    double sz = std::exp(scale_log(i, 2));

    arma::mat33 R_q = quat_to_rot(quat(i, 0), quat(i, 1), quat(i, 2), quat(i, 3));
    arma::mat33 cov3d_world = compute_cov3d(sx, sy, sz, R_q);
    arma::mat33 cov3d_cam = cam.R * cov3d_world * cam.R.t();
    arma::mat cov2d = J * cov3d_cam * J.t();

    cov2d(0, 0) += 0.3;
    cov2d(1, 1) += 0.3;

    double det = cov2d(0, 0) * cov2d(1, 1) - cov2d(0, 1) * cov2d(1, 0);
    if (det < 1e-7) continue;

    double inv_det = 1.0 / det;
    double inv_cov_xx = cov2d(1, 1) * inv_det;
    double inv_cov_xy = -cov2d(0, 1) * inv_det;
    double inv_cov_yy = cov2d(0, 0) * inv_det;

    double mid = 0.5 * (cov2d(0, 0) + cov2d(1, 1));
    double term = 0.25 * (cov2d(0, 0) - cov2d(1, 1)) * (cov2d(0, 0) - cov2d(1, 1)) + cov2d(0, 1) * cov2d(0, 1);
    double lambda_max = mid + std::sqrt(std::max(0.01, term));
    int radius = std::max(1, (int)std::ceil(3.0 * std::sqrt(lambda_max)));

    if (u + radius < 0 || u - radius >= w || v + radius < 0 || v - radius >= h) continue;

    double opac = sigmoid(opac_logit[i]);
    if (opac < 0.005) continue;

    double cr = 0.5 + SH_C0 * sh_mat(i, 0);
    double cg = 0.5 + SH_C0 * sh_mat(i, 1);
    double cb = 0.5 + SH_C0 * sh_mat(i, 2);

    Splat2D splat;
    splat.gaussian_idx = i;
    splat.u = u; splat.v = v;
    splat.depth = tz;
    splat.inv_cov_xx = inv_cov_xx;
    splat.inv_cov_xy = inv_cov_xy;
    splat.inv_cov_yy = inv_cov_yy;
    splat.radius = radius;
    splat.opacity = opac;
    splat.r = clamp(cr, 0.0, 1.0);
    splat.g = clamp(cg, 0.0, 1.0);
    splat.b = clamp(cb, 0.0, 1.0);

    visible_splats.push_back(splat);
    p_cams.push_back(p_cam);
  }

  int num_vis = (int)visible_splats.size();
  if (num_vis == 0) {
    return List::create(
      Named("gaussians") = gaussians,
      Named("adam_m") = adam_m,
      Named("adam_v") = adam_v,
      Named("loss") = 0.0,
      Named("grad_2d_norm") = NumericVector(num_gauss, 0.0)
    );
  }

  // 2. Tile binning & sorting
  const int TILE_SIZE = 16;
  int tiles_x = (w + TILE_SIZE - 1) / TILE_SIZE;
  int tiles_y = (h + TILE_SIZE - 1) / TILE_SIZE;
  int total_tiles = tiles_x * tiles_y;

  std::vector<SortItem> sort_items;
  sort_items.reserve(num_vis * 4);

  for (int s_idx = 0; s_idx < num_vis; ++s_idx) {
    const Splat2D& s = visible_splats[s_idx];
    int min_tx = std::max(0, (int)std::floor((s.u - s.radius) / TILE_SIZE));
    int max_tx = std::min(tiles_x - 1, (int)std::floor((s.u + s.radius) / TILE_SIZE));
    int min_ty = std::max(0, (int)std::floor((s.v - s.radius) / TILE_SIZE));
    int max_ty = std::min(tiles_y - 1, (int)std::floor((s.v + s.radius) / TILE_SIZE));

    for (int ty = min_ty; ty <= max_ty; ++ty) {
      for (int tx = min_tx; tx <= max_tx; ++tx) {
        uint32_t t_id = (uint32_t)(ty * tiles_x + tx);
        sort_items.push_back({t_id, (float)s.depth, s_idx});
      }
    }
  }

  std::sort(sort_items.begin(), sort_items.end());

  std::vector<int> tile_ranges(total_tiles + 1, (int)sort_items.size());
  uint32_t curr_tile = 0;
  tile_ranges[0] = 0;
  for (size_t k = 0; k < sort_items.size(); ++k) {
    uint32_t t_id = sort_items[k].tile_id;
    while (curr_tile < t_id) {
      curr_tile++;
      tile_ranges[curr_tile] = (int)k;
    }
  }

  // Gradient accumulators
  std::vector<double> grad_pos_x(num_gauss, 0.0);
  std::vector<double> grad_pos_y(num_gauss, 0.0);
  std::vector<double> grad_pos_z(num_gauss, 0.0);
  std::vector<double> grad_scale_x(num_gauss, 0.0);
  std::vector<double> grad_scale_y(num_gauss, 0.0);
  std::vector<double> grad_scale_z(num_gauss, 0.0);
  std::vector<double> grad_opac(num_gauss, 0.0);
  std::vector<double> grad_sh0(num_gauss, 0.0);
  std::vector<double> grad_sh1(num_gauss, 0.0);
  std::vector<double> grad_sh2(num_gauss, 0.0);
  NumericVector grad_2d_norm(num_gauss, 0.0);

  double total_loss = 0.0;
  double inv_pixel_count = 1.0 / (double)(w * h);

  // 3. Tile Rasterization + Photometric Gradient Backpropagation (OpenMP)
#pragma omp parallel for schedule(dynamic, 4) reduction(+:total_loss)
  for (int t_idx = 0; t_idx < total_tiles; ++t_idx) {
    int start_item = tile_ranges[t_idx];
    int end_item = (t_idx + 1 < (int)tile_ranges.size()) ? tile_ranges[t_idx + 1] : (int)sort_items.size();
    if (start_item >= end_item) continue;

    int tx = t_idx % tiles_x;
    int ty = t_idx / tiles_x;

    int px_start = tx * TILE_SIZE;
    int px_end = std::min(w, px_start + TILE_SIZE);
    int py_start = ty * TILE_SIZE;
    int py_end = std::min(h, py_start + TILE_SIZE);

    int num_tile_items = end_item - start_item;
    struct TileGrad {
      float dP_x = 0.0f;
      float dP_y = 0.0f;
      float dP_z = 0.0f;
      float d_sh_r = 0.0f;
      float d_sh_g = 0.0f;
      float d_sh_b = 0.0f;
      float d_opac = 0.0f;
      float d_scale = 0.0f;
    };
    std::vector<TileGrad> tile_grads(num_tile_items);

    for (int py = py_start; py < py_end; ++py) {
      for (int px = px_start; px < px_end; ++px) {
        // Forward pass per pixel
        double T = 1.0;
        double acc_r = 0.0, acc_g = 0.0, acc_b = 0.0;

        // Stack-allocated cache for contributing splats at this pixel (zero heap malloc)
        struct SplatContrib {
          int local_idx;
          float alpha;
          float weight;
          float dx;
          float dy;
        };
        const int MAX_CONTRIBS = 64;
        SplatContrib contribs[MAX_CONTRIBS];
        int num_contribs = 0;

        for (int k = start_item; k < end_item; ++k) {
          int s_idx = sort_items[k].splat_idx;
          const Splat2D& splat = visible_splats[s_idx];

          double dx = (double)px - splat.u;
          double dy = (double)py - splat.v;

          double power = -0.5 * (dx * dx * splat.inv_cov_xx +
                                 2.0 * dx * dy * splat.inv_cov_xy +
                                 dy * dy * splat.inv_cov_yy);

          if (power > 0.0) continue;

          double alpha = std::min(0.99, splat.opacity * std::exp(power));
          if (alpha < 1.0 / 255.0) continue;

          double weight = alpha * T;
          acc_r += splat.r * weight;
          acc_g += splat.g * weight;
          acc_b += splat.b * weight;

          if (num_contribs < MAX_CONTRIBS) {
            contribs[num_contribs].local_idx = k - start_item;
            contribs[num_contribs].alpha = (float)alpha;
            contribs[num_contribs].weight = (float)weight;
            contribs[num_contribs].dx = (float)dx;
            contribs[num_contribs].dy = (float)dy;
            num_contribs++;
          }

          T *= (1.0 - alpha);
          if (T < 0.001) break;
        }

        // Ground-truth target pixel: target_img layout (w * h * 3) or (h * w * 3)
        int pix_idx = px + py * w;
        double gt_r = (pix_idx < target_img.size()) ? target_img[pix_idx] : 0.0;
        double gt_g = (pix_idx + w * h < target_img.size()) ? target_img[pix_idx + w * h] : 0.0;
        double gt_b = (pix_idx + 2 * w * h < target_img.size()) ? target_img[pix_idx + 2 * w * h] : 0.0;
        if (gt_r > 1.0 || gt_g > 1.0 || gt_b > 1.0) {
          gt_r /= 255.0;
          gt_g /= 255.0;
          gt_b /= 255.0;
        }

        double err_r = acc_r - gt_r;
        double err_g = acc_g - gt_g;
        double err_b = acc_b - gt_b;

        total_loss += (std::abs(err_r) + std::abs(err_g) + std::abs(err_b)) * inv_pixel_count;

        // Smooth photometric gradient dL/dI (Huber / Smooth-L1 derivative)
        double dL_r = clamp(err_r, -1.0, 1.0) * inv_pixel_count;
        double dL_g = clamp(err_g, -1.0, 1.0) * inv_pixel_count;
        double dL_b = clamp(err_b, -1.0, 1.0) * inv_pixel_count;

        // Backward pass through contributing splats into fast tile-local buffer
        double accum_color_grad = 0.0;
        for (int c_idx = num_contribs - 1; c_idx >= 0; --c_idx) {
          const auto& c = contribs[c_idx];
          int l_idx = c.local_idx;
          int s_idx = sort_items[start_item + l_idx].splat_idx;
          const Splat2D& splat = visible_splats[s_idx];

          // dL / dColor
          double d_sh_r = dL_r * c.weight * SH_C0;
          double d_sh_g = dL_g * c.weight * SH_C0;
          double d_sh_b = dL_b * c.weight * SH_C0;

          // dL / dAlpha
          double splat_dot_grad = (splat.r * dL_r + splat.g * dL_g + splat.b * dL_b);
          double d_alpha = (c.weight / std::max(1e-5f, c.alpha)) * splat_dot_grad - accum_color_grad;
          accum_color_grad += c.weight * splat_dot_grad;

          // dAlpha / d(Position 2D)
          double d_u = d_alpha * c.alpha * (c.dx * splat.inv_cov_xx + c.dy * splat.inv_cov_xy);
          double d_v = d_alpha * c.alpha * (c.dx * splat.inv_cov_xy + c.dy * splat.inv_cov_yy);

          // Backpropagate to 3D position using camera projection Jacobian
          const arma::vec3& p_cam = p_cams[s_idx];
          double tz = p_cam(2);
          double tz2 = tz * tz;
          double dt_x = d_u * (cam.fx / tz);
          double dt_y = d_v * (cam.fy / tz);
          double dt_z = d_u * (-(cam.fx * p_cam(0)) / tz2) + d_v * (-(cam.fy * p_cam(1)) / tz2);

          // World 3D position gradient = R^T * dt
          arma::vec3 dt_cam = {dt_x, dt_y, dt_z};
          arma::vec3 dP_world = cam.R.t() * dt_cam;

          // dAlpha / d(Opacity logit): sigma'(o) = sigma(o) * (1 - sigma(o))
          double d_opac = d_alpha * (1.0 - splat.opacity);

          tile_grads[l_idx].dP_x += (float)dP_world(0);
          tile_grads[l_idx].dP_y += (float)dP_world(1);
          tile_grads[l_idx].dP_z += (float)dP_world(2);
          tile_grads[l_idx].d_sh_r += (float)d_sh_r;
          tile_grads[l_idx].d_sh_g += (float)d_sh_g;
          tile_grads[l_idx].d_sh_b += (float)d_sh_b;
          tile_grads[l_idx].d_opac += (float)d_opac;
          tile_grads[l_idx].d_scale += (float)(d_alpha * 0.001);
        }
      }
    }

    // High-efficiency atomic accumulation once per tile (replaces millions of per-pixel atomic locks)
    for (int l_idx = 0; l_idx < num_tile_items; ++l_idx) {
      const auto& tg = tile_grads[l_idx];
      if (tg.dP_x == 0.0f && tg.dP_y == 0.0f && tg.dP_z == 0.0f && tg.d_opac == 0.0f) continue;
      int s_idx = sort_items[start_item + l_idx].splat_idx;
      int g_idx = visible_splats[s_idx].gaussian_idx;

#pragma omp atomic
      grad_pos_x[g_idx] += tg.dP_x;
#pragma omp atomic
      grad_pos_y[g_idx] += tg.dP_y;
#pragma omp atomic
      grad_pos_z[g_idx] += tg.dP_z;
#pragma omp atomic
      grad_sh0[g_idx] += tg.d_sh_r;
#pragma omp atomic
      grad_sh1[g_idx] += tg.d_sh_g;
#pragma omp atomic
      grad_sh2[g_idx] += tg.d_sh_b;
#pragma omp atomic
      grad_opac[g_idx] += tg.d_opac;
#pragma omp atomic
      grad_scale_x[g_idx] += tg.d_scale;
#pragma omp atomic
      grad_scale_y[g_idx] += tg.d_scale;
#pragma omp atomic
      grad_scale_z[g_idx] += tg.d_scale;
    }
  }

  // 4. Adam Optimizer Parameter Updates
  const double beta1 = 0.9;
  const double beta2 = 0.999;
  const double eps = 1e-8;
  double beta1_t = std::pow(beta1, (double)iter);
  double beta2_t = std::pow(beta2, (double)iter);
  double corr1 = 1.0 - beta1_t;

#pragma omp parallel for schedule(static, 1024)
  for (int i = 0; i < num_gauss; ++i) {
    // 2D view-space gradient norm for adaptive densification
    double gx = grad_pos_x[i], gy = grad_pos_y[i], gz = grad_pos_z[i];
    grad_2d_norm[i] = std::sqrt(gx * gx + gy * gy + gz * gz);

    // Update 3D Positions
    double g_p[3] = {gx, gy, gz};
    for (int d = 0; d < 3; ++d) {
      m_pos(i, d) = beta1 * m_pos(i, d) + (1.0 - beta1) * g_p[d];
      v_pos(i, d) = beta2 * v_pos(i, d) + (1.0 - beta2) * g_p[d] * g_p[d];
      double m_hat = m_pos(i, d) / corr1;
      double v_hat = v_pos(i, d) / (1.0 - beta2_t);
      pos(i, d) -= lr_pos * (m_hat / (std::sqrt(v_hat) + eps));
    }

    // Update Scales
    double g_s[3] = {grad_scale_x[i], grad_scale_y[i], grad_scale_z[i]};
    for (int d = 0; d < 3; ++d) {
      m_scale(i, d) = beta1 * m_scale(i, d) + (1.0 - beta1) * g_s[d];
      v_scale(i, d) = beta2 * v_scale(i, d) + (1.0 - beta2) * g_s[d] * g_s[d];
      double m_hat = m_scale(i, d) / corr1;
      double v_hat = v_scale(i, d) / (1.0 - beta2_t);
      scale_log(i, d) -= lr_scale * (m_hat / (std::sqrt(v_hat) + eps));
      // Clamp scale within sane botanical bounds [-8.0, 1.0]
      scale_log(i, d) = clamp(scale_log(i, d), -8.0, 0.5);
    }

    // Update Opacity Logit
    double g_o = grad_opac[i];
    m_opac[i] = beta1 * m_opac[i] + (1.0 - beta1) * g_o;
    v_opac[i] = beta2 * v_opac[i] + (1.0 - beta2) * g_o * g_o;
    double m_hat_o = m_opac[i] / corr1;
    double v_hat_o = v_opac[i] / (1.0 - beta2_t);
    opac_logit[i] -= lr_opac * (m_hat_o / (std::sqrt(v_hat_o) + eps));
    opac_logit[i] = clamp(opac_logit[i], -5.0, 5.0);

    // Update Colors (SH0)
    double g_sh[3] = {grad_sh0[i], grad_sh1[i], grad_sh2[i]};
    for (int d = 0; d < 3; ++d) {
      m_sh(i, d) = beta1 * m_sh(i, d) + (1.0 - beta1) * g_sh[d];
      v_sh(i, d) = beta2 * v_sh(i, d) + (1.0 - beta2) * g_sh[d] * g_sh[d];
      double m_hat = m_sh(i, d) / corr1;
      double v_hat = v_sh(i, d) / (1.0 - beta2_t);
      sh_mat(i, d) -= lr_sh * (m_hat / (std::sqrt(v_hat) + eps));
    }
  }

  // Update gaussians list
  gaussians["positions"] = pos;
  gaussians["scales_log"] = scale_log;
  gaussians["quaternions"] = quat;
  gaussians["opacities_logit"] = opac_logit;
  gaussians["sh_coeffs"] = sh_mat;

  return List::create(
    Named("gaussians") = gaussians,
    Named("adam_m") = adam_m,
    Named("adam_v") = adam_v,
    Named("loss") = total_loss,
    Named("grad_2d_norm") = grad_2d_norm
  );
}

// ==============================================================================
// High-Performance 3D Kd-Tree and Statistical Outlier Removal (SOR / ROR)
// 100% Native C++ (OpenMP accelerated, zero external package dependencies)
// ==============================================================================

struct KdNode3D {
  float x, y, z;
  int point_idx;
  int axis;
  int left;
  int right;
};

class KdTree3D {
public:
  struct PointEntry {
    float x, y, z;
    int idx;
  };

  std::vector<PointEntry> pts;
  std::vector<KdNode3D> nodes;

  void build(const NumericMatrix& pos) {
    int n = pos.nrow();
    pts.resize(n);
    for (int i = 0; i < n; ++i) {
      pts[i] = {static_cast<float>(pos(i, 0)),
                static_cast<float>(pos(i, 1)),
                static_cast<float>(pos(i, 2)),
                i};
    }
    nodes.resize(n);
    int node_counter = 0;
    if (n > 0) {
      build_recursive(0, n, 0, node_counter);
    }
  }

  int build_recursive(int start, int end, int depth, int& node_counter) {
    if (start >= end) return -1;
    int axis = depth % 3;
    int mid = start + (end - start) / 2;

    auto comp = [axis](const PointEntry& a, const PointEntry& b) {
      if (axis == 0) return a.x < b.x;
      if (axis == 1) return a.y < b.y;
      return a.z < b.z;
    };

    std::nth_element(pts.begin() + start, pts.begin() + mid, pts.begin() + end, comp);

    int cur_idx = node_counter++;
    nodes[cur_idx].x = pts[mid].x;
    nodes[cur_idx].y = pts[mid].y;
    nodes[cur_idx].z = pts[mid].z;
    nodes[cur_idx].point_idx = pts[mid].idx;
    nodes[cur_idx].axis = axis;

    nodes[cur_idx].left = build_recursive(start, mid, depth + 1, node_counter);
    nodes[cur_idx].right = build_recursive(mid + 1, end, depth + 1, node_counter);

    return cur_idx;
  }

  // Find k nearest neighbors using a max-heap of squared distances
  void knn_search(float qx, float qy, float qz, int k, bool exclude_self, int query_idx,
                  std::priority_queue<std::pair<float, int>>& heap, int node_idx = 0) const {
    if (node_idx < 0 || node_idx >= static_cast<int>(nodes.size())) return;

    const auto& node = nodes[node_idx];

    float dx = qx - node.x;
    float dy = qy - node.y;
    float dz = qz - node.z;
    float d2 = dx * dx + dy * dy + dz * dz;

    if (!exclude_self || node.point_idx != query_idx) {
      if (static_cast<int>(heap.size()) < k) {
        heap.push({d2, node.point_idx});
      } else if (d2 < heap.top().first) {
        heap.pop();
        heap.push({d2, node.point_idx});
      }
    }

    float diff = (node.axis == 0) ? dx : ((node.axis == 1) ? dy : dz);
    int first = (diff <= 0.0f) ? node.left : node.right;
    int second = (diff <= 0.0f) ? node.right : node.left;

    knn_search(qx, qy, qz, k, exclude_self, query_idx, heap, first);

    float diff2 = diff * diff;
    if (static_cast<int>(heap.size()) < k || diff2 < heap.top().first) {
      knn_search(qx, qy, qz, k, exclude_self, query_idx, heap, second);
    }
  }

  // Count neighbors within radius r
  int radius_count(float qx, float qy, float qz, float r2, bool exclude_self, int query_idx, int node_idx = 0) const {
    if (node_idx < 0 || node_idx >= static_cast<int>(nodes.size())) return 0;

    const auto& node = nodes[node_idx];
    float dx = qx - node.x;
    float dy = qy - node.y;
    float dz = qz - node.z;
    float d2 = dx * dx + dy * dy + dz * dz;

    int count = 0;
    if ((!exclude_self || node.point_idx != query_idx) && d2 <= r2) {
      count = 1;
    }

    float diff = (node.axis == 0) ? dx : ((node.axis == 1) ? dy : dz);
    int first = (diff <= 0.0f) ? node.left : node.right;
    int second = (diff <= 0.0f) ? node.right : node.left;

    count += radius_count(qx, qy, qz, r2, exclude_self, query_idx, first);

    if (diff * diff <= r2) {
      count += radius_count(qx, qy, qz, r2, exclude_self, query_idx, second);
    }
    return count;
  }
};

//' @title Statistical and Radius Outlier Removal (SOR / ROR) for 3D Point Clouds
//' @name cpp_gs_filter_outliers
//' @description
//' High-speed native C++ Spatial Kd-Tree implementation of Statistical Outlier Removal (SOR)
//' and Radius Outlier Removal (ROR). Filters out sparse, floating airborne noise and background
//' artifacts from 3D point clouds without any external library dependencies.
//' @param points N x 3 NumericMatrix of 3D coordinates.
//' @param k Integer. Number of nearest neighbors to evaluate (default 20).
//' @param std_ratio Double. Standard deviation multiplier threshold (default 1.0).
//' @param radius Double. Optional radius for spherical density filtering (0 = disabled).
//' @param min_neighbors Integer. Minimum neighbors inside radius (0 = disabled).
//' @return A List containing 'keep' (LogicalVector), 'mean_dists' (NumericVector),
//'   'cutoff' (double), 'global_mean' (double), 'global_sd' (double), 'num_kept' (int), and 'num_removed' (int).
//' @export
// [[Rcpp::export]]
List cpp_gs_filter_outliers(NumericMatrix points,
                            int k = 20,
                            double std_ratio = 1.0,
                            double radius = 0.0,
                            int min_neighbors = 0) {
  int n = points.nrow();
  if (n == 0) {
    return List::create(
      Named("keep") = LogicalVector(0),
      Named("mean_dists") = NumericVector(0),
      Named("cutoff") = 0.0,
      Named("global_mean") = 0.0,
      Named("global_sd") = 0.0,
      Named("num_kept") = 0,
      Named("num_removed") = 0
    );
  }

  // 1. Build Spatial 3D Kd-Tree
  KdTree3D tree;
  tree.build(points);

  int k_eval = std::min(std::max(1, k), n - 1);
  NumericVector mean_dists(n);
  LogicalVector keep(n, true);

  // 2. Parallel k-NN query across points
#pragma omp parallel for schedule(dynamic, 512)
  for (int i = 0; i < n; ++i) {
    float qx = static_cast<float>(points(i, 0));
    float qy = static_cast<float>(points(i, 1));
    float qz = static_cast<float>(points(i, 2));

    std::priority_queue<std::pair<float, int>> heap;
    tree.knn_search(qx, qy, qz, k_eval, true, i, heap);

    double sum_d = 0.0;
    int count = 0;
    while (!heap.empty()) {
      sum_d += std::sqrt(static_cast<double>(heap.top().first));
      heap.pop();
      count++;
    }
    mean_dists[i] = (count > 0) ? (sum_d / count) : 0.0;
  }

  // 3. Compute global mean and standard deviation of neighbor distances
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += mean_dists[i];
  }
  double global_mean = sum / n;

  double sq_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    double diff = mean_dists[i] - global_mean;
    sq_sum += diff * diff;
  }
  double global_sd = (n > 1) ? std::sqrt(sq_sum / (n - 1)) : 0.0;
  double cutoff = global_mean + std_ratio * global_sd;

  // 4. Thresholding & Optional Radius check
  float r2 = static_cast<float>(radius * radius);
  int kept_count = 0;

  for (int i = 0; i < n; ++i) {
    bool is_inlier = (mean_dists[i] <= cutoff);

    if (is_inlier && radius > 0.0 && min_neighbors > 0) {
      float qx = static_cast<float>(points(i, 0));
      float qy = static_cast<float>(points(i, 1));
      float qz = static_cast<float>(points(i, 2));
      int neighbors = tree.radius_count(qx, qy, qz, r2, true, i);
      if (neighbors < min_neighbors) {
        is_inlier = false;
      }
    }

    keep[i] = is_inlier;
    if (is_inlier) kept_count++;
  }

  return List::create(
    Named("keep") = keep,
    Named("mean_dists") = mean_dists,
    Named("cutoff") = cutoff,
    Named("global_mean") = global_mean,
    Named("global_sd") = global_sd,
    Named("num_kept") = kept_count,
    Named("num_removed") = (n - kept_count)
  );
}

//' @title Fast 3D k-Nearest Neighbors Query
//' @name cpp_gs_knn
//' @description
//' High-speed native Kd-Tree k-nearest neighbors search in 3D Euclidean space.
//' Completely replaces external packages like RANN with 100% native C++ code.
//' @param points N x 3 NumericMatrix of reference point coordinates.
//' @param queries M x 3 NumericMatrix of query coordinates (or NULL to query 'points' against itself).
//' @param k Integer. Number of nearest neighbors.
//' @return A List containing 'indices' (M x k IntegerMatrix, 1-based) and 'dists' (M x k NumericMatrix).
//' @export
// [[Rcpp::export]]
List cpp_gs_knn(NumericMatrix points, Nullable<NumericMatrix> queries = R_NilValue, int k = 10) {
  int n = points.nrow();
  if (n == 0) {
    return List::create(
      Named("indices") = IntegerMatrix(0, 0),
      Named("dists") = NumericMatrix(0, 0)
    );
  }

  KdTree3D tree;
  tree.build(points);

  bool self_query = queries.isNull();
  NumericMatrix q_mat = self_query ? points : as<NumericMatrix>(queries);
  int m = q_mat.nrow();
  int k_eval = std::min(std::max(1, k), self_query ? (n - 1) : n);

  IntegerMatrix indices(m, k_eval);
  NumericMatrix dists(m, k_eval);

#pragma omp parallel for schedule(dynamic, 512)
  for (int i = 0; i < m; ++i) {
    float qx = static_cast<float>(q_mat(i, 0));
    float qy = static_cast<float>(q_mat(i, 1));
    float qz = static_cast<float>(q_mat(i, 2));

    std::priority_queue<std::pair<float, int>> heap;
    tree.knn_search(qx, qy, qz, k_eval, self_query, i, heap);

    int pos = k_eval - 1;
    while (!heap.empty()) {
      if (pos >= 0) {
        indices(i, pos) = heap.top().second + 1; // 1-based R index
        dists(i, pos) = std::sqrt(static_cast<double>(heap.top().first));
      }
      heap.pop();
      pos--;
    }
  }

  return List::create(
    Named("indices") = indices,
    Named("dists") = dists
  );
}
