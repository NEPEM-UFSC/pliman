// ==============================================================================
// 3D GAUSSIAN SPLATTING ENGINE FOR R (pliman)
// High-Performance Differentiable Point-Based 3D Reconstruction
// Authors: Tiago Olivoto & pliman team
// ==============================================================================

#ifndef PLIMAN_GAUSSIAN_SPLATTING_H
#define PLIMAN_GAUSSIAN_SPLATTING_H

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace pliman_gs {

// Constant Spherical Harmonics normalization coefficients
constexpr double SH_C0 = 0.28209479177387814;
constexpr double SH_C1 = 0.4886025119029199;
constexpr double SH_C2_0 = 1.0925484305920792;
constexpr double SH_C2_1 = -1.0925484305920792;
constexpr double SH_C2_2 = 0.31539156525252005;
constexpr double SH_C2_3 = -1.0925484305920792;
constexpr double SH_C2_4 = 0.5462742152960396;

inline double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

inline double sigmoid_inv(double y) {
  y = std::max(1e-6, std::min(1.0 - 1e-6, y));
  return std::log(y / (1.0 - y));
}

inline double clamp(double val, double min_v, double max_v) {
  return std::max(min_v, std::min(max_v, val));
}

// 3D Gaussian Primitive
struct Gaussian {
  double x, y, z;               // Position (mean mu)
  double scale_log_x;           // Log-scale x
  double scale_log_y;           // Log-scale y
  double scale_log_z;           // Log-scale z
  double qw, qx, qy, qz;        // Orientation quaternion
  double opacity_logit;         // Opacity logit
  std::vector<double> sh;       // Spherical harmonics (size 3, 12, 27, or 48)
  int semantic_id;              // Semantic class (0 = default, 1 = leaf, 2 = stem, etc.)
  double semantic_score;        // Probability / confidence
};

// Projected 2D Splat (view-space cache for rasterization)
struct Splat2D {
  int gaussian_idx;
  double u, v;                  // 2D screen coordinate
  double depth;                 // Camera z depth
  double inv_cov_xx;            // 2D inverse covariance xx
  double inv_cov_xy;            // 2D inverse covariance xy
  double inv_cov_yy;            // 2D inverse covariance yy
  int radius;                   // 3-sigma bounding radius in pixels
  double opacity;               // Sigmoid opacity
  double r, g, b;               // Evaluated RGB color in [0, 1]
  int semantic_id;              // Semantic ID
  double semantic_score;        // Semantic score
};

// Camera specification
struct Camera {
  arma::mat33 R;                // World-to-Camera rotation (3x3)
  arma::vec3 t;                 // World-to-Camera translation (3x1)
  double fx, fy;                // Focal lengths
  double cx, cy;                // Principal point
  int width, height;            // Viewport dimensions
  arma::vec3 pos_world;         // Camera position in world coordinates (-R^T * t)
};

// Convert quaternion to 3x3 rotation matrix
inline arma::mat33 quat_to_rot(double qw, double qx, double qy, double qz) {
  double norm = std::sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
  if (norm < 1e-9) norm = 1.0;
  double w = qw / norm;
  double x = qx / norm;
  double y = qy / norm;
  double z = qz / norm;

  arma::mat33 R;
  R(0, 0) = 1.0 - 2.0 * (y * y + z * z);
  R(0, 1) = 2.0 * (x * y - w * z);
  R(0, 2) = 2.0 * (x * z + w * y);

  R(1, 0) = 2.0 * (x * y + w * z);
  R(1, 1) = 1.0 - 2.0 * (x * x + z * z);
  R(1, 2) = 2.0 * (y * z - w * x);

  R(2, 0) = 2.0 * (x * z - w * y);
  R(2, 1) = 2.0 * (y * z + w * x);
  R(2, 2) = 1.0 - 2.0 * (x * x + y * y);
  return R;
}

// Compute 3D Covariance Sigma = R * S * S^T * R^T
inline arma::mat33 compute_cov3d(double sx, double sy, double sz, const arma::mat33& R) {
  arma::mat33 S;
  S.zeros();
  S(0, 0) = sx;
  S(1, 1) = sy;
  S(2, 2) = sz;
  arma::mat33 M = R * S;
  return M * M.t();
}

// Evaluate Spherical Harmonics for direction vector
inline void eval_sh(int degree, const std::vector<double>& sh, double dir_x, double dir_y, double dir_z,
                    double& r, double& g, double& b) {
  double norm = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
  if (norm > 1e-8) {
    dir_x /= norm;
    dir_y /= norm;
    dir_z /= norm;
  }

  // Degree 0 (ambient color)
  r = 0.5 + SH_C0 * sh[0];
  g = 0.5 + SH_C0 * sh[1];
  b = 0.5 + SH_C0 * sh[2];

  if (degree >= 1 && sh.size() >= 12) {
    double y1 = -SH_C1 * dir_y;
    double y2 =  SH_C1 * dir_z;
    double y3 = -SH_C1 * dir_x;

    r += y1 * sh[3]  + y2 * sh[6]  + y3 * sh[9];
    g += y1 * sh[4]  + y2 * sh[7]  + y3 * sh[10];
    b += y1 * sh[5]  + y2 * sh[8]  + y3 * sh[11];
  }

  if (degree >= 2 && sh.size() >= 27) {
    double xx = dir_x * dir_x, yy = dir_y * dir_y, zz = dir_z * dir_z;
    double xy = dir_x * dir_y, yz = dir_y * dir_z, xz = dir_x * dir_z;

    double y4 = SH_C2_0 * xy;
    double y5 = SH_C2_1 * yz;
    double y6 = SH_C2_2 * (2.0 * zz - xx - yy);
    double y7 = SH_C2_3 * xz;
    double y8 = SH_C2_4 * (xx - yy);

    r += y4 * sh[12] + y5 * sh[15] + y6 * sh[18] + y7 * sh[21] + y8 * sh[24];
    g += y4 * sh[13] + y5 * sh[16] + y6 * sh[19] + y7 * sh[22] + y8 * sh[25];
    b += y4 * sh[14] + y5 * sh[17] + y6 * sh[20] + y7 * sh[23] + y8 * sh[26];
  }

  r = clamp(r, 0.0, 1.0);
  g = clamp(g, 0.0, 1.0);
  b = clamp(b, 0.0, 1.0);
}

} // namespace pliman_gs

#endif // PLIMAN_GAUSSIAN_SPLATTING_H
