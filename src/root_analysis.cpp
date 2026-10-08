// root_analysis.cpp — Advanced Root System Phenotyping Engine for pliman
// Exact Euclidean Distance Transform, Guo-Hall Topological Skeletonization,
// Graph Topology Extraction, Primary vs Lateral Classification, Branching Angles,
// Depth Layering, and RhizoVision-Surpassing Phenotypic Traits.
// Zero Python, Zero OpenCV, Zero External Dependencies.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// 1. Exact Euclidean Distance Transform (Felzenszwalb-Huttenlocher Algorithm)
// -----------------------------------------------------------------------------
static const double INF_DIST = 1e8;

// 1D squared distance transform along an array of size n
static void edt_1d(const double* f, double* d, int* v, double* z, int n) {
  int k = 0;
  v[0] = 0;
  z[0] = -INF_DIST;
  z[1] = +INF_DIST;

  for (int q = 1; q < n; ++q) {
    double s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * (double)q - 2.0 * (double)v[k]);
    while (s <= z[k]) {
      k--;
      s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * (double)q - 2.0 * (double)v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k + 1] = +INF_DIST;
  }

  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[k + 1] < (double)q) {
      k++;
    }
    double diff = (double)q - (double)v[k];
    d[q] = diff * diff + f[v[k]];
  }
}

// 2D exact Euclidean distance transform on binary image (1 = foreground, 0 = background)
// Returns distance to nearest background pixel for each foreground pixel.
static std::vector<double> compute_edt_2d(const uint8_t* bin, int w, int h) {
  std::vector<double> f(w * h);
  std::vector<double> d(w * h);

  // Initialize: 0 for background, infinity for foreground
  for (int i = 0; i < w * h; ++i) {
    f[i] = (bin[i] == 0) ? 0.0 : INF_DIST;
  }

  // Transform along columns (x fixed, vary y)
  // Note: in pliman, index is x + y * w
  #pragma omp parallel
  {
    std::vector<double> f_col(h);
    std::vector<double> d_col(h);
    std::vector<int> v_col(h);
    std::vector<double> z_col(h + 1);

    #pragma omp for schedule(static)
    for (int x = 0; x < w; ++x) {
      for (int y = 0; y < h; ++y) {
        f_col[y] = f[x + y * w];
      }
      edt_1d(f_col.data(), d_col.data(), v_col.data(), z_col.data(), h);
      for (int y = 0; y < h; ++y) {
        d[x + y * w] = d_col[y];
      }
    }
  }

  // Transform along rows (y fixed, vary x)
  #pragma omp parallel
  {
    std::vector<double> d_row(w);
    std::vector<double> out_row(w);
    std::vector<int> v_row(w);
    std::vector<double> z_row(w + 1);

    #pragma omp for schedule(static)
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        d_row[x] = d[x + y * w];
      }
      edt_1d(d_row.data(), out_row.data(), v_row.data(), z_row.data(), w);
      for (int x = 0; x < w; ++x) {
        f[x + y * w] = out_row[x];
      }
    }
  }

  // Take square root for Euclidean distance
  std::vector<double> dist(w * h);
  for (int i = 0; i < w * h; ++i) {
    if (bin[i] == 0) {
      dist[i] = 0.0;
    } else {
      dist[i] = std::sqrt(f[i]);
    }
  }
  return dist;
}

// -----------------------------------------------------------------------------
// 2. Fast Guo-Hall Skeletonization
// -----------------------------------------------------------------------------
static std::vector<uint8_t> guo_hall_skeleton(const uint8_t* bin, int w, int h) {
  std::vector<uint8_t> grid(w * h);
  for (int i = 0; i < w * h; ++i) grid[i] = (bin[i] > 0) ? 1 : 0;

  bool changed = true;
  while (changed) {
    changed = false;
    for (int sub = 0; sub < 2; ++sub) {
      bool even = (sub == 1);
      std::vector<int> to_remove;

      for (int y = 1; y < h - 1; ++y) {
        for (int x = 1; x < w - 1; ++x) {
          int k = x + y * w;
          if (grid[k] == 0) continue;

          int p2 = grid[(x) + (y - 1) * w];
          int p3 = grid[(x + 1) + (y - 1) * w];
          int p4 = grid[(x + 1) + y * w];
          int p5 = grid[(x + 1) + (y + 1) * w];
          int p6 = grid[x + (y + 1) * w];
          int p7 = grid[(x - 1) + (y + 1) * w];
          int p8 = grid[(x - 1) + y * w];
          int p9 = grid[(x - 1) + (y - 1) * w];

          int C = ((!p2) & (p3 | p4)) + ((!p4) & (p5 | p6)) + ((!p6) & (p7 | p8)) + ((!p8) & (p9 | p2));
          if (C != 1) continue;

          int N1 = (p9 | p2) + (p3 | p4) + (p5 | p6) + (p7 | p8);
          int N2 = (p2 | p3) + (p4 | p5) + (p6 | p7) + (p8 | p9);
          int N = (N1 < N2) ? N1 : N2;
          if (N < 2 || N > 3) continue;

          int m = even ? ((p6 | p7 | (!p9)) & p8) : ((p2 | p3 | (!p5)) & p4);
          if (m == 0) {
            to_remove.push_back(k);
          }
        }
      }

      if (!to_remove.empty()) {
        changed = true;
        for (int idx : to_remove) grid[idx] = 0;
      }
    }
  }
  return grid;
}

// -----------------------------------------------------------------------------
// 3. 8-Neighbor Connectivity & Skeleton Degree
// -----------------------------------------------------------------------------
static inline int get_skel_degree(const uint8_t* skel, int w, int h, int x, int y) {
  int deg = 0;
  for (int dy = -1; dy <= 1; ++dy) {
    for (int dx = -1; dx <= 1; ++dx) {
      if (dx == 0 && dy == 0) continue;
      int nx = x + dx;
      int ny = y + dy;
      if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
        if (skel[nx + ny * w] > 0) deg++;
      }
    }
  }
  return deg;
}

// -----------------------------------------------------------------------------
// 4. Radius-Aware Adaptive Topological Spur Pruning
// Prunes terminal branches whose length <= max(min_length_px, parent_radius * radius_factor).
// Completely eliminates bark hairs and surface roughness without shortening genuine lateral roots.
// -----------------------------------------------------------------------------
static std::vector<uint8_t> adaptive_prune_skel(const uint8_t* skel_in,
                                                const double* dist_map,
                                                int w, int h,
                                                int min_length_px,
                                                double radius_factor = 1.2) {
  std::vector<uint8_t> skel(w * h);
  std::copy(skel_in, skel_in + w * h, skel.begin());

  int min_y = 1e9;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (skel[x + y * w] > 0) {
        if (y < min_y) min_y = y;
        break;
      }
    }
  }

  bool changed = true;
  int iter = 0;

  while (changed && iter < 3) {
    changed = false;
    iter++;
    std::vector<int> spurs_to_remove;

    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        int idx = x + y * w;
        if (skel[idx] == 0) continue;

        int deg = get_skel_degree(skel.data(), w, h, x, y);
        if (deg == 1) {
          // Trace along branch
          std::vector<int> branch_pixels;
          branch_pixels.push_back(idx);

          int curr = idx;
          int prev = -1;
          bool hit_junction = false;
          int junc_pix = -1;

          for (int step = 0; step < 250; ++step) {
            int cx = curr % w;
            int cy = curr / w;

            int next_pix = -1;
            for (int dy = -1; dy <= 1; ++dy) {
              for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;
                int nx = cx + dx;
                int ny = cy + dy;
                if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
                  int nidx = nx + ny * w;
                  if (skel[nidx] > 0 && nidx != prev) {
                    next_pix = nidx;
                    break;
                  }
                }
              }
              if (next_pix != -1) break;
            }

            if (next_pix == -1) break;

            int next_deg = get_skel_degree(skel.data(), w, h, next_pix % w, next_pix / w);
            if (next_deg >= 3) {
              hit_junction = true;
              junc_pix = next_pix;
              break;
            }

            branch_pixels.push_back(next_pix);
            prev = curr;
            curr = next_pix;
          }

          if (hit_junction && junc_pix >= 0) {
            // A genuine spur has radius <= 2.0 px (surface hair/noise on the bark).
            // Any branch with thick pixels (radius > 2.2 px) or near the crown collar is part of the root system!
            double max_b_r = 0.0;
            double sum_b_r = 0.0;
            for (int bpix : branch_pixels) {
              double r = (dist_map != nullptr) ? dist_map[bpix] : 1.0;
              max_b_r = std::max(max_b_r, r);
              sum_b_r += r;
            }
            double avg_b_r = sum_b_r / (double)branch_pixels.size();

            if (max_b_r > 2.2 || avg_b_r > 1.8) continue; // Thick root trunk/branch, not a spur!

            int tip_y = idx / w;
            if (tip_y <= min_y + 15) continue; // Crown collar tip, never prune!

            double junc_r = (dist_map != nullptr) ? dist_map[junc_pix] : 0.0;
            double max_len = std::max((double)min_length_px, junc_r * radius_factor);
            if ((double)branch_pixels.size() <= max_len) {
              for (int bpix : branch_pixels) {
                spurs_to_remove.push_back(bpix);
              }
            }
          }
        }
      }
    }

    if (!spurs_to_remove.empty()) {
      changed = true;
      for (int sidx : spurs_to_remove) {
        skel[sidx] = 0;
      }
    }
  }

  return skel;
}

// -----------------------------------------------------------------------------
// 5. Convex Hull Helper (Monotone Chain Algorithm)
// -----------------------------------------------------------------------------
struct Point2D {
  double x, y;
};

static bool compare_pts(const Point2D& a, const Point2D& b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

static double cross_product(const Point2D& o, const Point2D& a, const Point2D& b) {
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

static double compute_convex_hull_area(std::vector<Point2D>& pts) {
  int n = (int)pts.size();
  if (n < 3) return 0.0;

  std::sort(pts.begin(), pts.end(), compare_pts);
  std::vector<Point2D> h(2 * n);
  int k = 0;

  // Lower hull
  for (int i = 0; i < n; ++i) {
    while (k >= 2 && cross_product(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
    h[k++] = pts[i];
  }

  // Upper hull
  for (int i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t && cross_product(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
    h[k++] = pts[i];
  }
  h.resize(k - 1);

  // Shoelace formula
  double area = 0.0;
  int m = (int)h.size();
  for (int i = 0; i < m; ++i) {
    int j = (i + 1) % m;
    area += h[i].x * h[j].y - h[j].x * h[i].y;
  }
  return std::abs(area) * 0.5;
}

// -----------------------------------------------------------------------------
// 6. Graph Topology Extraction & Root System Analysis
// -----------------------------------------------------------------------------
struct Node {
  int id;
  double x, y;
  int type; // 1 = Tip, 2 = Junction, 3 = Crown
  std::vector<int> edge_ids;
};

struct Edge {
  int id;
  int node1, node2;
  std::vector<Point2D> coords;
  double length_px;
  double length_mm;
  double avg_diameter_mm;
  double median_diameter_mm;
  double max_diameter_mm;
  double min_diameter_mm;
  double surface_area_mm2;
  double volume_mm3;
  int order; // 1 = Primary root, 2 = Lateral root, 3 = Tertiary
  double branching_angle; // degrees relative to parent root
  double insertion_depth;
};

/// -----------------------------------------------------------------------------
// 4b. Binary Solidification Helpers
// -----------------------------------------------------------------------------
static std::vector<uint8_t> fill_small_holes(const std::vector<uint8_t>& bin, int w, int h, int max_size = 200) {
  std::vector<uint8_t> filled = bin;
  std::vector<int> bg_label(w * h, -1);
  int label = 0;

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int idx = x + y * w;
      if (bin[idx] == 0 && bg_label[idx] == -1) {
        std::vector<int> comp;
        std::queue<int> q;
        bool touches_border = false;
        q.push(idx);
        bg_label[idx] = label;

        while (!q.empty()) {
          int curr = q.front(); q.pop();
          comp.push_back(curr);
          int cx = curr % w, cy = curr / w;
          if (cx == 0 || cx == w - 1 || cy == 0 || cy == h - 1) touches_border = true;

          const int dx[] = {1, -1, 0, 0};
          const int dy[] = {0, 0, 1, -1};
          for (int k = 0; k < 4; ++k) {
            int nx = cx + dx[k], ny = cy + dy[k];
            if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
              int nidx = nx + ny * w;
              if (bin[nidx] == 0 && bg_label[nidx] == -1) {
                bg_label[nidx] = label;
                q.push(nidx);
              }
            }
          }
        }

        if (!touches_border && (int)comp.size() <= max_size) {
          for (int p : comp) filled[p] = 1;
        }
        label++;
      }
    }
  }
  return filled;
}

static std::vector<uint8_t> morph_closing(const std::vector<uint8_t>& bin, int w, int h, int radius) {
  if (radius <= 0) return bin;
  std::vector<uint8_t> inv(w * h);
  for (int i = 0; i < w * h; ++i) inv[i] = (bin[i] > 0) ? 0 : 1;
  std::vector<double> dt_fg = compute_edt_2d(inv.data(), w, h);

  std::vector<uint8_t> dilated(w * h, 0);
  for (int i = 0; i < w * h; ++i) {
    if (bin[i] > 0 || dt_fg[i] <= (double)radius) dilated[i] = 1;
  }

  std::vector<double> dt_bg = compute_edt_2d(dilated.data(), w, h);
  std::vector<uint8_t> closed(w * h, 0);
  for (int i = 0; i < w * h; ++i) {
    if (dt_bg[i] > (double)radius) closed[i] = 1;
  }
  return closed;
}

static void draw_line_8conn(Rcpp::NumericMatrix& mat, int x0, int y0, int x1, int y1) {
  int w = mat.nrow();
  int h = mat.ncol();
  int dx = std::abs(x1 - x0);
  int dy = std::abs(y1 - y0);
  int sx = (x0 < x1) ? 1 : -1;
  int sy = (y0 < y1) ? 1 : -1;
  int err = dx - dy;

  while (true) {
    if (x0 >= 0 && x0 < w && y0 >= 0 && y0 < h) {
      mat(x0, y0) = 1.0;
    }
    if (x0 == x1 && y0 == y1) break;
    int e2 = 2 * err;
    if (e2 > -dy) {
      err -= dy;
      x0 += sx;
    }
    if (e2 < dx) {
      err += dx;
      y0 += sy;
    }
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix skeleton_mat_cpp(Rcpp::NumericMatrix binary_mat,
                                    int closing_rad = 1,
                                    int fill_size = 200,
                                    int min_len = 10,
                                    double rad_fac = 1.2,
                                    bool dissolve_cycles = true,
                                    bool smooth = true) {
  int w = binary_mat.nrow();
  int h = binary_mat.ncol();
  int npix = w * h;

  std::vector<uint8_t> bin(npix, 0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (binary_mat(x, y) > 0.5) bin[x + y * w] = 1;
    }
  }

  std::vector<uint8_t> filled = (fill_size > 0) ? fill_small_holes(bin, w, h, fill_size) : bin;
  std::vector<uint8_t> closed = (closing_rad > 0) ? morph_closing(filled, w, h, closing_rad) : filled;
  std::vector<double> edt = compute_edt_2d(closed.data(), w, h);

  int pw = w + 4, ph = h + 4;
  std::vector<uint8_t> pbin(pw * ph, 0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (closed[x + y * w] > 0) pbin[(x + 2) + (y + 2) * pw] = 1;
    }
  }
  std::vector<uint8_t> raw_pskel = guo_hall_skeleton(pbin.data(), pw, ph);
  std::vector<uint8_t> raw_skel(npix, 0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      raw_skel[x + y * w] = raw_pskel[(x + 2) + (y + 2) * pw];
    }
  }

  std::vector<uint8_t> skel = adaptive_prune_skel(raw_skel.data(), edt.data(), w, h, min_len, rad_fac);

  std::vector<int> deg(npix, 0);
  std::vector<int> skel_pixels;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int idx = x + y * w;
      if (skel[idx]) {
        deg[idx] = get_skel_degree(skel.data(), w, h, x, y);
        skel_pixels.push_back(idx);
      }
    }
  }

  std::vector<int> junction_cluster(npix, -1);
  std::map<int, int> pixel_to_node;
  int next_node_id = 0;

  for (int idx : skel_pixels) {
    if (deg[idx] >= 3 && junction_cluster[idx] == -1) {
      std::queue<int> q;
      q.push(idx);
      junction_cluster[idx] = next_node_id;

      while (!q.empty()) {
        int curr = q.front(); q.pop();
        pixel_to_node[curr] = next_node_id;
        int cx = curr % w, cy = curr / w;
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            if (!dx && !dy) continue;
            int nx = cx + dx, ny = cy + dy;
            if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
              int nidx = nx + ny * w;
              if (skel[nidx] && deg[nidx] >= 3 && junction_cluster[nidx] == -1) {
                junction_cluster[nidx] = next_node_id;
                q.push(nidx);
              }
            }
          }
        }
      }
      next_node_id++;
    }
  }

  for (int idx : skel_pixels) {
    if (deg[idx] == 1) {
      pixel_to_node[idx] = next_node_id++;
    }
  }

  int total_nodes = next_node_id;

  struct MatEdge {
    int u, v;
    std::vector<Point2D> coords;
    double weight;
    double length_px;
    double avg_radius;
  };

  std::vector<MatEdge> edges;
  std::vector<bool> visited_edge_pix(npix, false);
  std::set<std::pair<int, int>> visited_transitions;

  for (int nid = 0; nid < total_nodes; ++nid) {
    std::vector<int> seeds;
    for (const auto& p2n : pixel_to_node) {
      if (p2n.second == nid) seeds.push_back(p2n.first);
    }

    for (int sidx : seeds) {
      int sx = sidx % w, sy = sidx / w;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (!dx && !dy) continue;
          int nx = sx + dx, ny = sy + dy;
          if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
          int nidx = nx + ny * w;
          if (!skel[nidx]) continue;
          if (pixel_to_node.count(nidx) && pixel_to_node[nidx] == nid) continue;
          if (visited_edge_pix[nidx] && deg[nidx] == 2) continue;

          std::vector<Point2D> path;
          path.push_back({(double)sx, (double)sy});

          int curr_x = nx, curr_y = ny, curr_idx = nidx, prev_idx = sidx;
          int target_node_id = -1;

          while (true) {
            path.push_back({(double)curr_x, (double)curr_y});
            visited_edge_pix[curr_idx] = true;

            if (pixel_to_node.count(curr_idx)) {
              target_node_id = pixel_to_node[curr_idx];
              break;
            }

            int next_x = -1, next_y = -1, next_idx = -1;
            for (int ddy = -1; ddy <= 1; ++ddy) {
              for (int ddx = -1; ddx <= 1; ++ddx) {
                if (!ddx && !ddy) continue;
                int tx = curr_x + ddx, ty = curr_y + ddy;
                if (tx < 0 || tx >= w || ty < 0 || ty >= h) continue;
                int tidx = tx + ty * w;
                if (skel[tidx] && tidx != prev_idx) {
                  next_x = tx; next_y = ty; next_idx = tidx;
                  break;
                }
              }
              if (next_idx != -1) break;
            }

            if (next_idx == -1) break;

            prev_idx = curr_idx;
            curr_x = next_x; curr_y = next_y; curr_idx = next_idx;
          }

          if (target_node_id != -1 && target_node_id != nid) {
            int na = std::min(nid, target_node_id);
            int nb = std::max(nid, target_node_id);
            if (!visited_transitions.count({na, nb})) {
              visited_transitions.insert({na, nb});

              double sum_r = 0.0, len_px = 0.0;
              for (size_t k = 0; k < path.size(); ++k) {
                int px = (int)std::round(path[k].x);
                int py = (int)std::round(path[k].y);
                if (px >= 0 && px < w && py >= 0 && py < h) sum_r += edt[px + py * w];
                if (k > 0) {
                  double ddx = path[k].x - path[k - 1].x;
                  double ddy = path[k].y - path[k - 1].y;
                  len_px += std::sqrt(ddx * ddx + ddy * ddy);
                }
              }

              double avg_r = sum_r / (double)path.size();
              double dx_e = std::abs(path.back().x - path.front().x);
              double dy_e = std::abs(path.back().y - path.front().y);
              double vert_ratio = (dy_e + 2.0) / (dx_e + 2.0);
              // Prioritize thickness and vertical downward continuity; penalize horizontal ladder rungs and long detour loops
              double wt = std::pow(std::max(0.1, avg_r), 3.5) * std::pow(vert_ratio, 1.5) / std::sqrt(std::max(1.0, len_px));

              MatEdge me;
              me.u = na;
              me.v = nb;
              me.coords = path;
              me.length_px = len_px;
              me.avg_radius = avg_r;
              me.weight = wt;
              edges.push_back(me);
            }
          }
        }
      }
    }
  }

  std::vector<bool> keep_edge(edges.size(), true);
  if (dissolve_cycles && !edges.empty()) {
    std::vector<int> edge_order(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) edge_order[i] = (int)i;

    std::sort(edge_order.begin(), edge_order.end(), [&](int a, int b) {
      return edges[a].weight > edges[b].weight;
    });

    std::vector<int> dsu(total_nodes);
    for (int i = 0; i < total_nodes; ++i) dsu[i] = i;
    std::function<int(int)> find_root = [&](int i) -> int {
      if (dsu[i] == i) return i;
      return dsu[i] = find_root(dsu[i]);
    };

    for (int eid : edge_order) {
      int ru = find_root(edges[eid].u);
      int rv = find_root(edges[eid].v);
      if (ru != rv) {
        dsu[ru] = rv;
      } else {
        keep_edge[eid] = false;
      }
    }
  }

  // Post-MST pruning
  // Post-MST pruning: remove only genuine thin hair spurs created by cycle dissolution
  bool pruned = true;
  int p_iter = 0;
  while (pruned && p_iter < 2) {
    pruned = false;
    p_iter++;
    std::vector<int> node_deg(total_nodes, 0);
    std::vector<std::vector<int>> node_inc(total_nodes);
    for (size_t i = 0; i < edges.size(); ++i) {
      if (keep_edge[i]) {
        node_deg[edges[i].u]++;
        node_deg[edges[i].v]++;
        node_inc[edges[i].u].push_back((int)i);
        node_inc[edges[i].v].push_back((int)i);
      }
    }

    for (size_t i = 0; i < edges.size(); ++i) {
      if (!keep_edge[i]) continue;
      if (edges[i].avg_radius > 1.8) continue; // Thick root segment, never a spur!

      int u = edges[i].u, v = edges[i].v;
      bool u_tip = (node_deg[u] == 1);
      bool v_tip = (node_deg[v] == 1);

      if (u_tip && !v_tip) {
        double max_parent_r = 0.0;
        for (int e_other : node_inc[v]) {
          if (e_other != (int)i && keep_edge[e_other]) max_parent_r = std::max(max_parent_r, edges[e_other].avg_radius);
        }
        if (edges[i].length_px <= std::max((double)min_len, max_parent_r * rad_fac)) {
          keep_edge[i] = false;
          pruned = true;
        }
      } else if (v_tip && !u_tip) {
        double max_parent_r = 0.0;
        for (int e_other : node_inc[u]) {
          if (e_other != (int)i && keep_edge[e_other]) max_parent_r = std::max(max_parent_r, edges[e_other].avg_radius);
        }
        if (edges[i].length_px <= std::max((double)min_len, max_parent_r * rad_fac)) {
          keep_edge[i] = false;
          pruned = true;
        }
      }
    }
  }

  NumericMatrix out(w, h);
  for (size_t i = 0; i < edges.size(); ++i) {
    if (!keep_edge[i]) continue;
    auto coords = edges[i].coords;
    int n_pts = (int)coords.size();
    if (n_pts < 2) continue;

    if (smooth && n_pts >= 5) {
      std::vector<Point2D> sm = coords;
      int hw = 2;
      for (int k = 1; k < n_pts - 1; ++k) {
        int k_lo = std::max(0, k - hw);
        int k_hi = std::min(n_pts - 1, k + hw);
        double sx = 0.0, sy = 0.0;
        for (int m = k_lo; m <= k_hi; ++m) {
          sx += coords[m].x;
          sy += coords[m].y;
        }
        sm[k].x = sx / (double)(k_hi - k_lo + 1);
        sm[k].y = sy / (double)(k_hi - k_lo + 1);
      }
      coords = sm;
    }

    for (int k = 0; k < n_pts; ++k) {
      int px = (int)std::round(coords[k].x);
      int py = (int)std::round(coords[k].y);
      if (k > 0) {
        int prev_px = (int)std::round(coords[k - 1].x);
        int prev_py = (int)std::round(coords[k - 1].y);
        draw_line_8conn(out, prev_px, prev_py, px, py);
      } else {
        if (px >= 0 && px < w && py >= 0 && py < h) {
          out(px, py) = 1.0;
        }
      }
    }
  }

  // Draw active junction cluster pixels to ensure seamless 8-connectivity at all branching nodes
  std::vector<bool> active_node(total_nodes, false);
  for (size_t i = 0; i < edges.size(); ++i) {
    if (keep_edge[i]) {
      if (edges[i].u >= 0 && edges[i].u < total_nodes) active_node[edges[i].u] = true;
      if (edges[i].v >= 0 && edges[i].v < total_nodes) active_node[edges[i].v] = true;
    }
  }

  for (const auto& p2n : pixel_to_node) {
    int nid = p2n.second;
    if (nid >= 0 && nid < total_nodes && active_node[nid]) {
      int idx = p2n.first;
      int px = idx % w;
      int py = idx / w;
      if (px >= 0 && px < w && py >= 0 && py < h) {
        out(px, py) = 1.0;
      }
    }
  }

  // Thin the connected graph with Guo-Hall to strictly guarantee 1-pixel width without holes
  std::vector<uint8_t> out_bin(npix, 0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (out(x, y) > 0.5) out_bin[x + y * w] = 1;
    }
  }
  std::vector<uint8_t> thinned = guo_hall_skeleton(out_bin.data(), w, h);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      out(x, y) = (double)thinned[x + y * w];
    }
  }

  return out;
}

static void merge_edge_coords(std::vector<Point2D>& c1, std::vector<Point2D> c2) {
  if (c2.empty()) return;
  if (c1.empty()) {
    c1 = c2;
    return;
  }
  double d_back_front  = std::hypot(c1.back().x - c2.front().x, c1.back().y - c2.front().y);
  double d_back_back   = std::hypot(c1.back().x - c2.back().x, c1.back().y - c2.back().y);
  double d_front_front = std::hypot(c1.front().x - c2.front().x, c1.front().y - c2.front().y);
  double d_front_back  = std::hypot(c1.front().x - c2.back().x, c1.front().y - c2.back().y);

  if (d_back_front <= d_back_back && d_back_front <= d_front_front && d_back_front <= d_front_back) {
    c1.insert(c1.end(), c2.begin(), c2.end());
  } else if (d_back_back <= d_back_front && d_back_back <= d_front_front && d_back_back <= d_front_back) {
    std::reverse(c2.begin(), c2.end());
    c1.insert(c1.end(), c2.begin(), c2.end());
  } else if (d_front_front <= d_back_front && d_front_front <= d_back_back && d_front_front <= d_front_back) {
    std::reverse(c1.begin(), c1.end());
    c1.insert(c1.end(), c2.begin(), c2.end());
  } else {
    std::reverse(c1.begin(), c1.end());
    std::reverse(c2.begin(), c2.end());
    c1.insert(c1.end(), c2.begin(), c2.end());
  }
}

// [[Rcpp::export]]
Rcpp::List analyze_root_system_cpp(Rcpp::NumericMatrix binary_mat,
                                  double pixel_size = 1.0,
                                  int prune_length_px = 12,
                                  double cluster_radius_px = 8.0,
                                  double crown_x_in = -1.0,
                                  double crown_y_in = -1.0,
                                  int num_classes = 4,
                                  Rcpp::Nullable<Rcpp::NumericVector> diameter_bins = R_NilValue,
                                  int n_depth_slices = 10,
                                  int max_order = 3,
                                  double gravitropic_alpha = 1.5,
                                  int closing_rad = 1,
                                  int fill_size = 200,
                                  double rad_fac = 1.2,
                                  bool dissolve_cycles = true) {

  int w = binary_mat.nrow();
  int h = binary_mat.ncol();
  int npix = w * h;

  // 1. Binary array (row = x, col = y)
  std::vector<uint8_t> bin(npix, 0);
  int fg_count = 0;
  std::vector<Point2D> fg_pts;
  fg_pts.reserve(npix / 4);

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (binary_mat(x, y) > 0.5) {
        bin[x + y * w] = 1;
        fg_count++;
        fg_pts.push_back({(double)x, (double)y});
      }
    }
  }

  if (fg_count == 0) {
    stop("analyze_root_system_cpp: Input image contains no foreground root pixels.");
  }

  // 1b. Solidify binary mask: fill small cavities and bridge micro bark cracks
  std::vector<uint8_t> filled = (fill_size > 0) ? fill_small_holes(bin, w, h, fill_size) : bin;
  std::vector<uint8_t> closed = (closing_rad > 0) ? morph_closing(filled, w, h, closing_rad) : filled;

  // 2. Exact Euclidean Distance Transform (EDT)
  std::vector<double> dist_map = compute_edt_2d(closed.data(), w, h);

  // 3. Obtain pristine, acyclic Medial Axis Tree skeleton via skeleton_mat_cpp
  NumericMatrix skel_res = skeleton_mat_cpp(binary_mat, closing_rad, fill_size, prune_length_px, rad_fac, dissolve_cycles, false);
  std::vector<uint8_t> skel(npix, 0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (skel_res(x, y) > 0.5) skel[x + y * w] = 1;
    }
  }

  // 4b. Crown detection: find topmost foreground pixel in largest connected component
  double target_crown_x = crown_x_in;
  double target_crown_y = crown_y_in;
  if (target_crown_x < 0 || target_crown_y < 0) {
    std::vector<int> cc_labels(npix, -1);
    std::vector<int> cc_sizes;
    int curr_label = 0;
    for (int i = 0; i < npix; ++i) {
      if (closed[i] && cc_labels[i] == -1) {
        std::queue<int> q;
        q.push(i);
        cc_labels[i] = curr_label;
        int sz = 0;
        while (!q.empty()) {
          int curr = q.front(); q.pop();
          sz++;
          int cx = curr % w;
          int cy = curr / w;
          for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
              if (dx == 0 && dy == 0) continue;
              int nx = cx + dx, ny = cy + dy;
              if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
                int nidx = nx + ny * w;
                if (closed[nidx] && cc_labels[nidx] == -1) {
                  cc_labels[nidx] = curr_label;
                  q.push(nidx);
                }
              }
            }
          }
        }
        cc_sizes.push_back(sz);
        curr_label++;
      }
    }
    int max_cc_label = -1, max_sz = 0;
    for (size_t l = 0; l < cc_sizes.size(); ++l) {
      if (cc_sizes[l] > max_sz) { max_sz = cc_sizes[l]; max_cc_label = (int)l; }
    }
    double min_y = 1e9;
    double best_cx = w / 2.0;
    for (int y = 0; y < h; ++y) {
      bool found = false;
      for (int x = 0; x < w; ++x) {
        int idx = x + y * w;
        if (closed[idx] && cc_labels[idx] == max_cc_label && dist_map[idx] >= 1.5) {
          if ((double)y < min_y) {
            min_y = (double)y;
            best_cx = (double)x;
            found = true;
          }
        }
      }
      if (found) break;
    }
    target_crown_x = best_cx;
    target_crown_y = min_y;
  }

  // 5. Degree calculation on pruned skeleton
  std::vector<int> deg(npix, 0);
  std::vector<int> skel_pixels;
  skel_pixels.reserve(fg_count / 5);

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int idx = x + y * w;
      if (skel[idx] > 0) {
        deg[idx] = get_skel_degree(skel.data(), w, h, x, y);
        skel_pixels.push_back(idx);
      }
    }
  }

  if (skel_pixels.empty()) {
    stop("analyze_root_system_cpp: Skeletonization resulted in empty skeleton. Check threshold or prune_length.");
  }

  // 6. Cluster adjacent junction pixels (deg >= 3) to form discrete junction nodes
  std::vector<int> junction_cluster(npix, -1);
  std::vector<Node> nodes;
  std::map<int, int> pixel_to_node;

  int next_node_id = 0;

  // Find all junction clusters using BFS
  for (int idx : skel_pixels) {
    if (deg[idx] >= 3 && junction_cluster[idx] == -1) {
      std::queue<int> q;
      q.push(idx);
      junction_cluster[idx] = next_node_id;

      double sum_x = 0.0, sum_y = 0.0;
      int count = 0;

      while (!q.empty()) {
        int curr = q.front();
        q.pop();
        int cx = curr % w;
        int cy = curr / w;
        sum_x += cx;
        sum_y += cy;
        count++;
        pixel_to_node[curr] = next_node_id;

        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = cx + dx;
            int ny = cy + dy;
            if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
              int nidx = nx + ny * w;
              if (skel[nidx] > 0 && deg[nidx] >= 3 && junction_cluster[nidx] == -1) {
                junction_cluster[nidx] = next_node_id;
                q.push(nidx);
              }
            }
          }
        }
      }

      Node jnode;
      jnode.id = next_node_id++;
      jnode.x = sum_x / count;
      jnode.y = sum_y / count;
      jnode.type = 2; // Junction
      nodes.push_back(jnode);
    }
  }

  // Find all tip nodes (deg == 1)
  for (int idx : skel_pixels) {
    if (deg[idx] == 1) {
      int tx = idx % w;
      int ty = idx / w;
      Node tnode;
      tnode.id = next_node_id++;
      tnode.x = (double)tx;
      tnode.y = (double)ty;
      tnode.type = 1; // Tip
      pixel_to_node[idx] = tnode.id;
      nodes.push_back(tnode);
    }
  }

  // 7. Early Crown Identification: mark node closest to (target_crown_x, target_crown_y)
  int early_crown_id = -1;
  double min_early_crown_dsq = 1e18;
  for (size_t i = 0; i < nodes.size(); ++i) {
    double dx = nodes[i].x - target_crown_x;
    double dy = nodes[i].y - target_crown_y;
    double dsq = dx * dx + dy * dy;
    if (dsq < min_early_crown_dsq) {
      min_early_crown_dsq = dsq;
      early_crown_id = (int)i;
    }
  }
  if (early_crown_id != -1) {
    nodes[early_crown_id].type = 3; // Crown
  }

  // 8. Trace Edges between Nodes along degree-2 pixels
  std::vector<Edge> edges;
  std::set<std::pair<int, int>> visited_transitions;
  std::vector<bool> visited_edge_pix(npix, false);

  int next_edge_id = 0;

  for (const auto& nd : nodes) {
    // For each node, find skeleton pixels belonging to this node or neighbor degree-2 pixels
    std::vector<int> start_seeds;
    if (nd.type == 1) { // Tip
      int px = (int)std::round(nd.x);
      int py = (int)std::round(nd.y);
      start_seeds.push_back(px + py * w);
    } else { // Junction or Crown
      for (const auto& p2n : pixel_to_node) {
        if (p2n.second == nd.id) {
          start_seeds.push_back(p2n.first);
        }
      }
    }

    for (int sidx : start_seeds) {
      int sx = sidx % w;
      int sy = sidx / w;

      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0) continue;
          int nx = sx + dx;
          int ny = sy + dy;
          if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
          int nidx = nx + ny * w;
          if (skel[nidx] == 0) continue;

          // Don't step inside the same node cluster
          if (pixel_to_node.count(nidx) && pixel_to_node[nidx] == nd.id) continue;

          if (visited_edge_pix[nidx] && deg[nidx] == 2) continue;

          // Trace edge path
          std::vector<Point2D> path;
          path.push_back({(double)sx, (double)sy});

          int curr_x = nx;
          int curr_y = ny;
          int curr_idx = nidx;
          int prev_idx = sidx;

          int target_node_id = -1;

          while (true) {
            path.push_back({(double)curr_x, (double)curr_y});
            visited_edge_pix[curr_idx] = true;

            // Check if we hit another node
            if (pixel_to_node.count(curr_idx)) {
              target_node_id = pixel_to_node[curr_idx];
              break;
            }

            // Find next unvisited step
            int next_x = -1, next_y = -1, next_idx = -1;
            for (int ddy = -1; ddy <= 1; ++ddy) {
              for (int ddx = -1; ddx <= 1; ++ddx) {
                if (ddx == 0 && ddy == 0) continue;
                int tx = curr_x + ddx;
                int ty = curr_y + ddy;
                if (tx < 0 || tx >= w || ty < 0 || ty >= h) continue;
                int tidx = tx + ty * w;
                if (skel[tidx] > 0 && tidx != prev_idx) {
                  next_x = tx;
                  next_y = ty;
                  next_idx = tidx;
                  break;
                }
              }
              if (next_idx != -1) break;
            }

            if (next_idx == -1) {
              // Dead end (isolated segment)
              break;
            }

            prev_idx = curr_idx;
            curr_x = next_x;
            curr_y = next_y;
            curr_idx = next_idx;
          }

          if (target_node_id != -1 && target_node_id != nd.id) {
            int n_a = std::min(nd.id, target_node_id);
            int n_b = std::max(nd.id, target_node_id);
            if (visited_transitions.count({n_a, n_b}) == 0) {
              visited_transitions.insert({n_a, n_b});

              Edge edge;
              edge.id = next_edge_id++;
              edge.node1 = nd.id;
              edge.node2 = target_node_id;
              edge.coords = path;

              // Geodesic length & diameters
              double len_px = 0.0;
              double sum_diam = 0.0;
              double sum_vol = 0.0;
              double sum_surf = 0.0;
              std::vector<double> diam_list;
              double max_d = 0.0;
              double min_d = 1e9;

              for (size_t i = 0; i < path.size(); ++i) {
                int px = (int)std::round(path[i].x);
                int py = (int)std::round(path[i].y);
                px = std::max(0, std::min(w - 1, px));
                py = std::max(0, std::min(h - 1, py));

                double r = dist_map[px + py * w];
                double diam_mm = 2.0 * r * pixel_size;
                diam_list.push_back(diam_mm);
                if (diam_mm > max_d) max_d = diam_mm;
                if (diam_mm < min_d) min_d = diam_mm;
                sum_diam += diam_mm;

                if (i > 0) {
                  double ddx = path[i].x - path[i - 1].x;
                  double ddy = path[i].y - path[i - 1].y;
                  double seg_len_px = std::sqrt(ddx * ddx + ddy * ddy);
                  double seg_len_mm = seg_len_px * pixel_size;
                  len_px += seg_len_px;

                  // Cylinder model integration
                  sum_surf += M_PI * diam_mm * seg_len_mm;
                  sum_vol += 0.25 * M_PI * diam_mm * diam_mm * seg_len_mm;
                }
              }

              if (path.size() == 1) {
                len_px = 1.0;
                double d = 2.0 * dist_map[(int)path[0].x + (int)path[0].y * w] * pixel_size;
                sum_surf = M_PI * d * pixel_size;
                sum_vol = 0.25 * M_PI * d * d * pixel_size;
              }

              edge.length_px = len_px;
              edge.length_mm = len_px * pixel_size;
              edge.avg_diameter_mm = (diam_list.empty()) ? 0.0 : (sum_diam / diam_list.size());
              
              std::sort(diam_list.begin(), diam_list.end());
              edge.median_diameter_mm = (diam_list.empty()) ? 0.0 : diam_list[diam_list.size() / 2];
              edge.max_diameter_mm = (diam_list.empty()) ? 0.0 : max_d;
              edge.min_diameter_mm = (diam_list.empty()) ? 0.0 : min_d;
              edge.surface_area_mm2 = sum_surf;
              edge.volume_mm3 = sum_vol;
              edge.order = 2; // Default lateral
              edge.branching_angle = NA_REAL;
              edge.insertion_depth = std::min(nd.y, nodes[target_node_id].y) * pixel_size;

              nodes[nd.id].edge_ids.push_back(edge.id);
              nodes[target_node_id].edge_ids.push_back(edge.id);
              edges.push_back(edge);
            }
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // 8.5 Graph Refinement: Iterative Spur Pruning, Degree-2 Dissolution & Junction Merging
  // ---------------------------------------------------------------------------
  std::vector<bool> edge_active(edges.size(), true);
  std::vector<bool> node_active(nodes.size(), true);

  auto get_active_edges = [&](int u) {
    std::vector<int> res;
    if (u < 0 || u >= (int)nodes.size()) return res;
    for (int eid : nodes[u].edge_ids) {
      if (eid >= 0 && eid < (int)edges.size() && edge_active[eid]) {
        if (edges[eid].node1 == u || edges[eid].node2 == u) {
          if (std::find(res.begin(), res.end(), eid) == res.end()) {
            res.push_back(eid);
          }
        }
      }
    }
    return res;
  };

  auto get_other_node = [&](int eid, int u) {
    return (edges[eid].node1 == u) ? edges[eid].node2 : edges[eid].node1;
  };

  bool graph_changed = true;
  int iter = 0;
  while (graph_changed && iter++ < 100) {
    graph_changed = false;

    // A. Spur Pruning: edges connecting a tip (type 1) to another node with length < prune_length_px
    for (size_t eid = 0; eid < edges.size(); ++eid) {
      if (!edge_active[eid]) continue;
      int u = edges[eid].node1;
      int v = edges[eid].node2;
      if (!node_active[u] || !node_active[v]) {
        edge_active[eid] = false;
        continue;
      }
      if (u == v) {
        edge_active[eid] = false;
        graph_changed = true;
        continue;
      }

      int type_u = nodes[u].type;
      int type_v = nodes[v].type;

      bool is_spur = false;
      int tip_to_remove = -1;

      if (type_u == 1 && type_v != 1) {
        if (nodes[u].type == 3 || nodes[v].type == 3) continue;
        if (nodes[u].y <= target_crown_y + 15.0) continue;
        is_spur = true;
        tip_to_remove = u;
      } else if (type_v == 1 && type_u != 1) {
        if (nodes[u].type == 3 || nodes[v].type == 3) continue;
        if (nodes[v].y <= target_crown_y + 15.0) continue;
        is_spur = true;
        tip_to_remove = v;
      } else if (type_u == 1 && type_v == 1) {
        // Isolated short fragment where both ends are tips
        if (edges[eid].length_px < (double)prune_length_px && edges[eid].avg_diameter_mm < 2.5 * pixel_size) {
          edge_active[eid] = false;
          node_active[u] = false;
          node_active[v] = false;
          graph_changed = true;
          continue;
        }
      }

      if (is_spur) {
        // Adaptive pruning threshold: a terminal branch shorter than the parent
        // junction's cross-section cannot be a genuine root — it never escaped the trunk.
        //
        // For a junction adjacent to a thick taproot:
        //   parent_radius_px  = sibling_edge_diameter / (2 * pixel_size)
        //   effective_thresh  = max(prune_length_px, 2.0 * parent_radius_px)
        //
        // This removes horizontal "stub" artifacts without touching genuine thin laterals
        // that branch from thick roots (those are long enough to survive the threshold).
        int junction_node = (tip_to_remove == u) ? v : u;
        double max_sibling_diam_mm = 0.0;
        for (int pe : get_active_edges(junction_node)) {
          if (pe != (int)eid && edge_active[pe]) {
            max_sibling_diam_mm = std::max(max_sibling_diam_mm,
                                           edges[pe].avg_diameter_mm);
          }
        }
        double parent_radius_px = max_sibling_diam_mm / (2.0 * std::max(pixel_size, 1e-9));
        double effective_thresh  = std::max((double)prune_length_px,
                                            2.0 * parent_radius_px);

        if (edges[eid].length_px >= effective_thresh) {
          continue; // Long enough relative to parent cross-section → genuine lateral
        }

        int other_node = junction_node;

        auto o_edges = get_active_edges(other_node);
        if (o_edges.size() <= 1 && nodes[other_node].type == 3) {
          // Don't prune sole root connected to crown
          continue;
        }
        edge_active[eid] = false;
        node_active[tip_to_remove] = false;
        graph_changed = true;
      }
    }

    // B. Refresh node active edges and degrees
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i]) continue;
      auto inc_edges = get_active_edges((int)i);
      nodes[i].edge_ids = inc_edges;
      int d = (int)inc_edges.size();
      if (d == 0) {
        node_active[i] = false;
        graph_changed = true;
      } else if (d == 1 && nodes[i].type != 3) {
        if (nodes[i].type != 1) {
          nodes[i].type = 1; // downgraded to Tip
          graph_changed = true;
        }
      }
    }

    // C. Dissolve degree-2 non-crown junction nodes
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i] || nodes[i].type == 3) continue;
      auto inc_edges = get_active_edges((int)i);
      nodes[i].edge_ids = inc_edges;

      if (inc_edges.size() == 2) {
        int e1 = inc_edges[0];
        int e2 = inc_edges[1];
        if (e1 != e2 && edge_active[e1] && edge_active[e2]) {
          int u1 = get_other_node(e1, (int)i);
          int u2 = get_other_node(e2, (int)i);

          if (u1 == u2) {
            // Parallel loop to the same neighbor node
            edge_active[e2] = false;
            node_active[i] = false;
            graph_changed = true;
          } else {
            edges[e1].node1 = u1;
            edges[e1].node2 = u2;
            edges[e1].length_px += edges[e2].length_px;
            edges[e1].length_mm += edges[e2].length_mm;
            edges[e1].surface_area_mm2 += edges[e2].surface_area_mm2;
            edges[e1].volume_mm3 += edges[e2].volume_mm3;
            if (edges[e1].length_mm + edges[e2].length_mm > 1e-6) {
              edges[e1].avg_diameter_mm = (edges[e1].avg_diameter_mm * edges[e1].length_mm +
                                           edges[e2].avg_diameter_mm * edges[e2].length_mm) /
                                          (edges[e1].length_mm + edges[e2].length_mm);
            }

            merge_edge_coords(edges[e1].coords, edges[e2].coords);

            edge_active[e2] = false;
            node_active[i] = false;

            nodes[u1].edge_ids.push_back(e1);
            nodes[u2].edge_ids.push_back(e1);
            graph_changed = true;
          }
        }
      }
    }

    // D. Cluster adjacent junctions connected by short internal edges
    if (cluster_radius_px > 1.0) {
      std::vector<int> parent_node(nodes.size());
      for (size_t i = 0; i < nodes.size(); ++i) parent_node[i] = (int)i;

      std::function<int(int)> find_set = [&](int x) -> int {
        if (parent_node[x] == x) return x;
        return parent_node[x] = find_set(parent_node[x]);
      };

      auto union_set = [&](int x, int y) {
        int rx = find_set(x);
        int ry = find_set(y);
        if (rx != ry) {
          if (nodes[ry].type == 3) parent_node[rx] = ry;
          else parent_node[ry] = rx;
        }
      };

      // Merge junctions connected by short internal edges
      for (size_t eid = 0; eid < edges.size(); ++eid) {
        if (!edge_active[eid]) continue;
        int u = edges[eid].node1;
        int v = edges[eid].node2;
        if (!node_active[u] || !node_active[v]) continue;
        if (nodes[u].type != 1 && nodes[v].type != 1) {
          // Never collapse thick main root segments or crown
          if (nodes[u].type == 3 || nodes[v].type == 3) continue;
          if (edges[eid].avg_diameter_mm > 2.5 * pixel_size) continue;
          if (edges[eid].length_px <= cluster_radius_px) {
            union_set(u, v);
            edge_active[eid] = false;
            graph_changed = true;
          }
        }
      }

      // Group clusters and update centroid
      std::map<int, std::vector<int>> clusters;
      for (size_t i = 0; i < nodes.size(); ++i) {
        if (!node_active[i] || nodes[i].type == 1) continue;
        int r = find_set((int)i);
        clusters[r].push_back((int)i);
      }

      for (const auto& pair : clusters) {
        int rep = pair.first;
        const auto& members = pair.second;
        if (members.size() > 1) {
          double sum_x = 0.0, sum_y = 0.0;
          bool is_crown = false;
          for (int m : members) {
            sum_x += nodes[m].x;
            sum_y += nodes[m].y;
            if (nodes[m].type == 3) is_crown = true;
            if (m != rep) node_active[m] = false;
          }
          nodes[rep].x = sum_x / (double)members.size();
          nodes[rep].y = sum_y / (double)members.size();
          if (is_crown) nodes[rep].type = 3;
        }
      }

      // Rewire edges
      for (size_t eid = 0; eid < edges.size(); ++eid) {
        if (!edge_active[eid]) continue;
        int u = edges[eid].node1;
        int v = edges[eid].node2;
        if (nodes[u].type != 1) edges[eid].node1 = find_set(u);
        if (nodes[v].type != 1) edges[eid].node2 = find_set(v);
        if (edges[eid].node1 == edges[eid].node2) {
          edge_active[eid] = false;
          graph_changed = true;
        }
      }

      // Collapse parallel edges between same pair of nodes
      std::map<std::pair<int, int>, int> edge_pairs;
      for (size_t eid = 0; eid < edges.size(); ++eid) {
        if (!edge_active[eid]) continue;
        int u = std::min(edges[eid].node1, edges[eid].node2);
        int v = std::max(edges[eid].node1, edges[eid].node2);
        std::pair<int, int> p = {u, v};
        if (edge_pairs.count(p)) {
          int old_eid = edge_pairs[p];
          edges[old_eid].length_px = std::max(edges[old_eid].length_px, edges[eid].length_px);
          edges[old_eid].length_mm = std::max(edges[old_eid].length_mm, edges[eid].length_mm);
          edge_active[eid] = false;
          graph_changed = true;
        } else {
          edge_pairs[p] = (int)eid;
        }
      }
    }
  }

  // 8B. Bridge small gaps between collinear endpoints (heal broken disconnected root segments)
  double max_gap_px = std::max(25.0, cluster_radius_px * 2.5);

  // Connected component tracking so we only bridge genuinely disconnected fragments
  std::vector<int> comp_parent(nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) comp_parent[i] = (int)i;
  std::function<int(int)> find_comp = [&](int x) -> int {
    if (comp_parent[x] == x) return x;
    return comp_parent[x] = find_comp(comp_parent[x]);
  };
  auto union_comp = [&](int x, int y) {
    int rx = find_comp(x);
    int ry = find_comp(y);
    if (rx != ry) comp_parent[rx] = ry;
  };

  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (edge_active[eid]) {
      union_comp(edges[eid].node1, edges[eid].node2);
    }
  }

  struct GapPair {
    int u;
    int v;
    double dist;
    double score;
  };
  std::vector<GapPair> candidate_gaps;

  // Find all active tip nodes
  std::vector<int> tip_nodes;
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (node_active[i] && nodes[i].type == 1) {
      tip_nodes.push_back((int)i);
    }
  }

  auto get_tip_tangent = [&](int tip_id, double& tx, double& ty) -> bool {
    auto u_edges = get_active_edges(tip_id);
    if (u_edges.empty()) return false;
    int eid = u_edges[0];
    const auto& c = edges[eid].coords;
    if (c.size() >= 2) {
      if (edges[eid].node1 == tip_id) {
        int idx = std::min((size_t)3, c.size() - 1);
        tx = c[0].x - c[idx].x;
        ty = c[0].y - c[idx].y;
      } else {
        int idx = (c.size() > 3) ? (int)c.size() - 4 : 0;
        tx = c.back().x - c[idx].x;
        ty = c.back().y - c[idx].y;
      }
    } else {
      int u_parent = get_other_node(eid, tip_id);
      tx = nodes[tip_id].x - nodes[u_parent].x;
      ty = nodes[tip_id].y - nodes[u_parent].y;
    }
    double norm_t = std::sqrt(tx * tx + ty * ty);
    if (norm_t > 1e-4) {
      tx /= norm_t;
      ty /= norm_t;
      return true;
    }
    return false;
  };

  for (size_t i = 0; i < tip_nodes.size(); ++i) {
    int u = tip_nodes[i];
    double tu_x = 0.0, tu_y = 0.0;
    if (!get_tip_tangent(u, tu_x, tu_y)) continue;

    for (size_t j = i + 1; j < tip_nodes.size(); ++j) {
      int v = tip_nodes[j];
      // Only bridge if they are in different disconnected components,
      // avoiding spurious cross-connections between adjacent branches of the same root tree
      if (find_comp(u) == find_comp(v)) continue;

      double tv_x = 0.0, tv_y = 0.0;
      if (!get_tip_tangent(v, tv_x, tv_y)) continue;

      double dx = nodes[v].x - nodes[u].x;
      double dy = nodes[v].y - nodes[u].y;
      double dist = std::sqrt(dx * dx + dy * dy);

      if (dist > 1e-4 && dist <= max_gap_px) {
        // Alignment: tu pointing towards v, tv pointing towards u
        double cos_u = (tu_x * dx + tu_y * dy) / dist;
        double cos_v = (tv_x * (-dx) + tv_y * (-dy)) / dist;

        if (cos_u > 0.15 && cos_v > 0.15) {
          double gap_score = (cos_u + cos_v) / (dist + 1.0);
          candidate_gaps.push_back({u, v, dist, gap_score});
        }
      }
    }
  }

  // Sort candidate gaps by best alignment / shortest distance
  std::sort(candidate_gaps.begin(), candidate_gaps.end(), [](const GapPair& a, const GapPair& b) {
    return a.score > b.score;
  });

  std::vector<bool> tip_bridged(nodes.size(), false);
  for (const auto& gap : candidate_gaps) {
    if (tip_bridged[gap.u] || tip_bridged[gap.v]) continue;

    tip_bridged[gap.u] = true;
    tip_bridged[gap.v] = true;

    // Create new bridging edge
    Edge bridge_ed;
    bridge_ed.id = (int)edges.size();
    bridge_ed.node1 = gap.u;
    bridge_ed.node2 = gap.v;
    bridge_ed.length_px = gap.dist;
    bridge_ed.length_mm = gap.dist * pixel_size;
    bridge_ed.order = 2; // Default, hierarchy may promote to 1
    bridge_ed.branching_angle = NA_REAL;
    bridge_ed.insertion_depth = std::min(nodes[gap.u].y, nodes[gap.v].y) * pixel_size;

    // Estimate diameter from adjacent edges
    double diam_u = 0.0, diam_v = 0.0;
    auto eu = get_active_edges(gap.u);
    if (!eu.empty()) diam_u = edges[eu[0]].avg_diameter_mm;
    auto ev = get_active_edges(gap.v);
    if (!ev.empty()) diam_v = edges[ev[0]].avg_diameter_mm;
    bridge_ed.avg_diameter_mm = (diam_u + diam_v) * 0.5;
    bridge_ed.median_diameter_mm = bridge_ed.avg_diameter_mm;
    bridge_ed.max_diameter_mm = std::max(diam_u, diam_v);
    bridge_ed.min_diameter_mm = std::min(diam_u, diam_v);
    bridge_ed.surface_area_mm2 = M_PI * bridge_ed.avg_diameter_mm * bridge_ed.length_mm;
    bridge_ed.volume_mm3 = 0.25 * M_PI * bridge_ed.avg_diameter_mm * bridge_ed.avg_diameter_mm * bridge_ed.length_mm;

    // Generate interpolated coordinates
    int n_steps = std::max(2, (int)std::ceil(gap.dist));
    for (int s = 0; s <= n_steps; ++s) {
      double f = (double)s / (double)n_steps;
      bridge_ed.coords.push_back({nodes[gap.u].x * (1.0 - f) + nodes[gap.v].x * f,
                                  nodes[gap.u].y * (1.0 - f) + nodes[gap.v].y * f});
    }

    edges.push_back(bridge_ed);
    edge_active.push_back(true);

    nodes[gap.u].type = 2; // Downgraded from Tip to Junction
    nodes[gap.v].type = 2;
    nodes[gap.u].edge_ids.push_back(bridge_ed.id);
    nodes[gap.v].edge_ids.push_back(bridge_ed.id);
  }

  // 8C. Cycle & Ladder Dissolution via Maximum Spanning Forest (Kruskal's Algorithm)
  // Roots branch in open trees; cycles in the graph represent either transverse ladder rungs
  // across thick roots or overlapping branches. Keeping the Maximum Spanning Forest (weighting by
  // diameter^2 * length) eliminates 100% of false loops and ladder rungs while preserving
  // true biological root continuity.
  if (dissolve_cycles && !edges.empty()) {
    std::vector<int> edge_order;
    for (size_t eid = 0; eid < edges.size(); ++eid) {
      if (edge_active[eid]) edge_order.push_back((int)eid);
    }

    auto calc_mst_weight = [&](int eid) -> double {
      int u = edges[eid].node1;
      int v = edges[eid].node2;
      double dx = std::abs(nodes[u].x - nodes[v].x);
      double dy = std::abs(nodes[u].y - nodes[v].y);
      double vert_ratio = (dy + 2.0) / (dx + 2.0);
      double D = std::max(0.1, edges[eid].avg_diameter_mm);
      double L = std::max(1.0, edges[eid].length_mm);
      return std::pow(D, 3.5) * std::pow(vert_ratio, 1.5) / std::sqrt(L);
    };

    std::sort(edge_order.begin(), edge_order.end(), [&](int a, int b) {
      return calc_mst_weight(a) > calc_mst_weight(b);
    });

    std::vector<int> dsu_parent(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) dsu_parent[i] = (int)i;
    std::function<int(int)> dsu_find = [&](int i) -> int {
      if (dsu_parent[i] == i) return i;
      return dsu_parent[i] = dsu_find(dsu_parent[i]);
    };

    for (int eid : edge_order) {
      int ru = dsu_find(edges[eid].node1);
      int rv = dsu_find(edges[eid].node2);
      if (ru != rv) {
        dsu_parent[ru] = rv;
      } else {
        // Discard edge that closes a cycle/ladder
        edge_active[eid] = false;
      }
    }

    // Post-MST refinement: prune any short spurs left from broken cycles
    bool post_pruned = true;
    int p_iter = 0;
    while (post_pruned && p_iter < 8) {
      post_pruned = false;
      p_iter++;

      for (size_t eid = 0; eid < edges.size(); ++eid) {
        if (!edge_active[eid]) continue;
        int u = edges[eid].node1;
        int v = edges[eid].node2;
        auto eu = get_active_edges(u);
        auto ev = get_active_edges(v);
        bool u_tip = (eu.size() == 1);
        bool v_tip = (ev.size() == 1);

        if (u_tip && !v_tip && nodes[u].type != 3 && nodes[u].y > target_crown_y + 15.0) {
          if (edges[eid].avg_diameter_mm / (2.0 * std::max(pixel_size, 1e-9)) > 1.8) continue; // Thick root segment, never a spur!
          double max_parent_r = 0.0;
          for (int pe : ev) {
            if (pe != (int)eid && edge_active[pe]) {
              max_parent_r = std::max(max_parent_r, edges[pe].avg_diameter_mm / (2.0 * std::max(pixel_size, 1e-9)));
            }
          }
          if (edges[eid].length_px <= std::max((double)prune_length_px, max_parent_r * rad_fac)) {
            edge_active[eid] = false;
            node_active[u] = false;
            post_pruned = true;
          }
        } else if (v_tip && !u_tip && nodes[v].type != 3 && nodes[v].y > target_crown_y + 15.0) {
          if (edges[eid].avg_diameter_mm / (2.0 * std::max(pixel_size, 1e-9)) > 1.8) continue; // Thick root segment, never a spur!
          double max_parent_r = 0.0;
          for (int pe : eu) {
            if (pe != (int)eid && edge_active[pe]) {
              max_parent_r = std::max(max_parent_r, edges[pe].avg_diameter_mm / (2.0 * std::max(pixel_size, 1e-9)));
            }
          }
          if (edges[eid].length_px <= std::max((double)prune_length_px, max_parent_r * rad_fac)) {
            edge_active[eid] = false;
            node_active[v] = false;
            post_pruned = true;
          }
        }
      }
    }

    // Dissolve degree-2 non-crown junction nodes after cycle breaking
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i] || nodes[i].type == 3) continue;
      auto inc_edges = get_active_edges((int)i);
      if (inc_edges.size() == 2) {
        int e1 = inc_edges[0];
        int e2 = inc_edges[1];
        if (e1 != e2 && edge_active[e1] && edge_active[e2]) {
          int u1 = get_other_node(e1, (int)i);
          int u2 = get_other_node(e2, (int)i);
          if (u1 != u2) {
            edges[e1].node1 = u1;
            edges[e1].node2 = u2;
            edges[e1].length_px += edges[e2].length_px;
            edges[e1].length_mm += edges[e2].length_mm;
            edges[e1].surface_area_mm2 += edges[e2].surface_area_mm2;
            edges[e1].volume_mm3 += edges[e2].volume_mm3;
            if (edges[e1].length_mm + edges[e2].length_mm > 1e-6) {
              edges[e1].avg_diameter_mm = (edges[e1].avg_diameter_mm * edges[e1].length_mm +
                                           edges[e2].avg_diameter_mm * edges[e2].length_mm) /
                                          (edges[e1].length_mm + edges[e2].length_mm);
            }
            merge_edge_coords(edges[e1].coords, edges[e2].coords);
            edge_active[e2] = false;
            node_active[i] = false;
          }
        }
      }
    }
  }

  // Smooth edge coordinates (5-point moving average)
  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (!edge_active[eid]) continue;
    int n_pts = (int)edges[eid].coords.size();
    if (n_pts >= 5) {
      std::vector<Point2D> sm = edges[eid].coords;
      int hw = 2;
      for (int k = 1; k < n_pts - 1; ++k) {
        int k_lo = std::max(0, k - hw);
        int k_hi = std::min(n_pts - 1, k + hw);
        double sx = 0.0, sy = 0.0;
        for (int m = k_lo; m <= k_hi; ++m) {
          sx += edges[eid].coords[m].x;
          sy += edges[eid].coords[m].y;
        }
        sm[k].x = sx / (double)(k_hi - k_lo + 1);
        sm[k].y = sy / (double)(k_hi - k_lo + 1);
      }
      edges[eid].coords = sm;
    }
  }

  // Final edge_ids rebuild and crown re-sync
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (!node_active[i]) continue;
    nodes[i].edge_ids = get_active_edges((int)i);
  }

  // Find connected components among active nodes to ensure crown belongs to the main root system
  std::vector<int> node_comp(nodes.size(), -1);
  std::vector<int> comp_size;
  int current_comp = 0;

  for (size_t i = 0; i < nodes.size(); ++i) {
    if (!node_active[i] || node_comp[i] != -1) continue;
    int c_id = current_comp++;
    int count = 0;
    std::queue<int> cq;
    cq.push((int)i);
    node_comp[i] = c_id;

    while (!cq.empty()) {
      int u = cq.front();
      cq.pop();
      count++;

      for (int eid : nodes[u].edge_ids) {
        if (!edge_active[eid]) continue;
        int v = (edges[eid].node1 == u) ? edges[eid].node2 : edges[eid].node1;
        if (node_active[v] && node_comp[v] == -1) {
          node_comp[v] = c_id;
          cq.push(v);
        }
      }
    }
    comp_size.push_back(count);
  }

  // Identify the largest connected component (main root system)
  int max_comp_id = -1;
  int max_comp_count = 0;
  for (size_t c = 0; c < comp_size.size(); ++c) {
    if (comp_size[c] > max_comp_count) {
      max_comp_count = comp_size[c];
      max_comp_id = (int)c;
    }
  }

  int crown_node_id = -1;
  if (crown_x_in >= 0 && crown_y_in >= 0) {
    // User specified crown: find closest active node in the largest component
    double min_dist_sq = 1e18;
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i]) continue;
      if (max_comp_id != -1 && node_comp[i] != max_comp_id) continue;
      double dsq = (nodes[i].x - crown_x_in) * (nodes[i].x - crown_x_in) + (nodes[i].y - crown_y_in) * (nodes[i].y - crown_y_in);
      if (dsq < min_dist_sq) {
        min_dist_sq = dsq;
        crown_node_id = (int)i;
      }
    }
  } else {
    // Automatic crown: find node in the largest component closest to (target_crown_x, target_crown_y)
    double min_dist_sq = 1e18;
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i]) continue;
      if (max_comp_id != -1 && node_comp[i] != max_comp_id) continue;
      double dsq = (nodes[i].x - target_crown_x) * (nodes[i].x - target_crown_x) +
                   (nodes[i].y - target_crown_y) * (nodes[i].y - target_crown_y);
      if (dsq < min_dist_sq) {
        min_dist_sq = dsq;
        crown_node_id = (int)i;
      }
    }
  }

  // Ensure all node types are reset: exactly one crown, all others tip or junction
  for (size_t i = 0; i < nodes.size(); ++i) {
    if ((int)i == crown_node_id) {
      nodes[i].type = 3;
    } else {
      nodes[i].type = (get_active_edges((int)i).size() <= 1) ? 1 : 2;
    }
  }

  // 9. Root Hierarchy (Primary vs. Lateral vs. Tertiary+) & Branching Angle Classification
  std::vector<int> node_order(nodes.size(), 0);
  int count_sec_axes = 0;
  int count_tert_axes = 0;
  if (crown_node_id != -1 && !edges.empty() && node_active[crown_node_id]) {
    // A. Build directed tree from crown via Dijkstra (shortest geodesic cost path)
    // Minimizes traversal cost, heavily penalizing upward growth, horizontal drift, and thin edges.
    std::priority_queue<std::pair<double, int>, 
                        std::vector<std::pair<double, int>>, 
                        std::greater<std::pair<double, int>>> pq;
    std::vector<double> dist(nodes.size(), 1e18);
    std::vector<int> parent_node(nodes.size(), -1);
    std::vector<int> parent_edge(nodes.size(), -1);
    std::vector<bool> tree_visited(nodes.size(), false);

    dist[crown_node_id] = 0.0;
    pq.push({0.0, crown_node_id});

    while (!pq.empty()) {
      auto top = pq.top();
      pq.pop();
      double d = top.first;
      int u = top.second;

      if (tree_visited[u]) continue;
      tree_visited[u] = true;

      for (int eid : nodes[u].edge_ids) {
        if (!edge_active[eid]) continue;
        int v = (edges[eid].node1 == u) ? edges[eid].node2 : edges[eid].node1;
        if (!node_active[v] || tree_visited[v]) continue;

        double dy = nodes[v].y - nodes[u].y;
        double dx = std::abs(nodes[v].x - nodes[u].x);
        double L = std::max(1.0, edges[eid].length_px);
        double D = std::max(0.1, edges[eid].avg_diameter_mm);

        // Direction penalty: roots grow downwards from crown to tip
        double dir_pen = 1.0;
        if (dy < 0.0) {
          dir_pen = 500.0 * (1.0 + std::abs(dy)); // Heavy penalty for growing upward
        } else {
          double vy = dy / L; // in [0, 1]
          double rxy = dx / std::max(1.0, dy);
          dir_pen = 1.0 + 4.0 * (1.0 - vy) * (1.0 - vy) + 2.0 * rxy;
        }

        // Edge traversal cost: thick, vertical downward edges have lowest cost
        double edge_cost = (L / std::pow(D, 2.5)) * dir_pen;

        if (d + edge_cost < dist[v]) {
          dist[v] = d + edge_cost;
          parent_node[v] = u;
          parent_edge[v] = eid;
          pq.push({dist[v], v});
        }
      }
    }

    // B. PRIMARY ROOT (Taproot) — Globally Optimal Path Selection
    // Evaluates every terminal tip T in the tree.
    // The true primary root maximizes:
    //   Score = (Depth)^2.5 * (MeanDiameter)^3 * (Verticality)^2 * Centrality * Straightness
    int best_tip = -1;
    double best_path_score = -1.0;
    std::vector<int> best_path_edges;

    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i] || !tree_visited[i] || (int)i == crown_node_id) continue;

      auto act_e = get_active_edges((int)i);
      if (act_e.size() > 1) continue; // Junction node, not a terminal tip

      // Trace path from tip i back to crown_node_id
      std::vector<int> cur_path;
      std::vector<int> cur_nodes;
      int curr = (int)i;
      double sum_diam_len = 0.0;
      double total_len_mm = 0.0;
      double total_len_px = 0.0;

      cur_nodes.push_back(curr);
      while (curr != crown_node_id && parent_edge[curr] != -1) {
        int pe = parent_edge[curr];
        cur_path.push_back(pe);
        sum_diam_len += edges[pe].avg_diameter_mm * edges[pe].length_mm;
        total_len_mm += edges[pe].length_mm;
        total_len_px += edges[pe].length_px;
        curr = parent_node[curr];
        cur_nodes.push_back(curr);
      }

      if (curr != crown_node_id || cur_path.empty()) continue;

      // Reverse so it goes from crown to tip
      std::reverse(cur_path.begin(), cur_path.end());
      std::reverse(cur_nodes.begin(), cur_nodes.end());

      double mean_diam = (total_len_mm > 1e-6) ? (sum_diam_len / total_len_mm) : 0.1;
      double dy = nodes[i].y - nodes[crown_node_id].y; // vertical downward depth
      double dx = std::abs(nodes[i].x - nodes[crown_node_id].x); // horizontal offset

      if (dy <= 0.0) continue; // Skip any tips growing above the crown

      double verticality = dy / std::max(1.0, total_len_px); // in [0, 1]
      double centrality = 1.0 / (1.0 + 2.0 * (dx / std::max(1.0, dy)));

      // End-to-end straightness and severe kink penalty
      double min_cos = 1.0;
      for (size_t k = 1; k + 1 < cur_nodes.size(); ++k) {
        int p_prev = cur_nodes[k - 1];
        int p_curr = cur_nodes[k];
        int p_next = cur_nodes[k + 1];

        double v1x = nodes[p_curr].x - nodes[p_prev].x;
        double v1y = nodes[p_curr].y - nodes[p_prev].y;
        double v2x = nodes[p_next].x - nodes[p_curr].x;
        double v2y = nodes[p_next].y - nodes[p_curr].y;

        double n1 = std::sqrt(v1x * v1x + v1y * v1y);
        double n2 = std::sqrt(v2x * v2x + v2y * v2y);
        if (n1 >= 4.0 && n2 >= 4.0) {
          double cos_turn = (v1x * v2x + v1y * v2y) / (n1 * n2);
          if (cos_turn < min_cos) min_cos = cos_turn;
        }
      }
      // Severely penalize paths that make sharp > 75 deg turns or 90-deg kinks (min_cos < 0.25)
      double kink_penalty = (min_cos < 0.25) ? std::max(0.05, 0.4 + 2.4 * min_cos) : 1.0;

      // Global Taproot Dominance Score: heavily favors deeper, thicker, vertically aligned and direct paths
      double score = std::pow(dy, 3.0) * 
                     std::pow(std::max(0.1, mean_diam), 2.5) * 
                     std::pow(std::max(0.01, verticality), 2.0) * 
                     centrality * 
                     kink_penalty;

      if (score > best_path_score) {
        best_path_score = score;
        best_tip = (int)i;
        best_path_edges = cur_path;
      }
    }

    // Reset edge orders to 0 before hierarchy classification
    for (size_t eid = 0; eid < edges.size(); ++eid) {
      edges[eid].order = 0;
    }

    // Apply primary axis assignment (order = 1)
    if (crown_node_id >= 0 && crown_node_id < (int)nodes.size()) {
      node_order[crown_node_id] = 1;
    }
    for (int eid : best_path_edges) {
      edges[eid].order = 1;
      node_order[edges[eid].node1] = 1;
      node_order[edges[eid].node2] = 1;
    }

    // D. Dominant Axis Continuation for Secondary, Tertiary and Higher-Order Roots
    // Secondary roots originate from the primary taproot and follow their dominant
    // (largest diameter and length) continuation downstream to their tips.
    // Tertiary+ roots branch off secondary roots and follow their own dominant continuation.

    // 1. Identify primary root nodes
    std::vector<int> primary_nodes;
    if (crown_node_id >= 0 && crown_node_id < (int)nodes.size()) {
      primary_nodes.push_back(crown_node_id);
    }
    for (int eid : best_path_edges) {
      primary_nodes.push_back(edges[eid].node1);
      primary_nodes.push_back(edges[eid].node2);
    }
    std::sort(primary_nodes.begin(), primary_nodes.end());
    primary_nodes.erase(std::unique(primary_nodes.begin(), primary_nodes.end()), primary_nodes.end());

    std::vector<bool> is_primary_node(nodes.size(), false);
    for (int u : primary_nodes) {
      if (u >= 0 && u < (int)nodes.size()) is_primary_node[u] = true;
    }

    // 2. Queue for iterative dominant axis processing
    struct BranchSeed {
      int eid;
      int from_node;
      int order;
    };
    std::queue<BranchSeed> branch_queue;

    // Seeds from primary root: all active edges emerging from primary_nodes not part of primary root
    for (int u : primary_nodes) {
      for (int eid : nodes[u].edge_ids) {
        if (edge_active[eid] && edges[eid].order == 0) {
          branch_queue.push({eid, u, 2});
        }
      }
    }

    // 3. Process branches iteratively: follow the dominant continuation along each axis
    while (!branch_queue.empty()) {
      BranchSeed seed = branch_queue.front();
      branch_queue.pop();

      int curr_eid = seed.eid;
      int from_node = seed.from_node;
      int curr_ord = std::min(max_order, seed.order);

      if (edges[curr_eid].order > 0) continue; // Already assigned

      if (curr_ord == 2) count_sec_axes++;
      else if (curr_ord >= 3) count_tert_axes++;

      while (curr_eid != -1) {
        if (edges[curr_eid].order > 0) break;
        edges[curr_eid].order = curr_ord;

        int to_node = (edges[curr_eid].node1 == from_node) ? edges[curr_eid].node2 : edges[curr_eid].node1;
        if (to_node < 0 || to_node >= (int)nodes.size() || !node_active[to_node]) break;
        if (is_primary_node[to_node] || nodes[to_node].type == 1) break; // Tip or returned to taproot

        // Find candidate outgoing edges from to_node (excluding curr_eid and edges into primary nodes)
        std::vector<int> candidates;
        for (int ne : nodes[to_node].edge_ids) {
          if (!edge_active[ne] || ne == curr_eid || edges[ne].order > 0) continue;
          int other = (edges[ne].node1 == to_node) ? edges[ne].node2 : edges[ne].node1;
          if (is_primary_node[other]) continue;
          candidates.push_back(ne);
        }

        if (candidates.empty()) break;

        // Dominant continuation: pick candidate with maximum diameter & length score
        // Score = diameter^2.5 * sqrt(length)
        int best_eid = -1;
        double max_score = -1.0;
        for (int ce : candidates) {
          double d = edges[ce].avg_diameter_mm;
          double l = edges[ce].length_mm;
          double score = std::pow(std::max(0.1, d), 2.5) * std::sqrt(std::max(0.1, l));
          if (score > max_score) {
            max_score = score;
            best_eid = ce;
          }
        }

        // All OTHER candidates at this junction branch off as order + 1
        for (int ce : candidates) {
          if (ce != best_eid) {
            branch_queue.push({ce, to_node, curr_ord + 1});
          }
        }

        curr_eid = best_eid;
        from_node = to_node;
      }
    }

    // 4. Fallback: Any active edges not reached (e.g. disconnected components)
    for (size_t eid = 0; eid < edges.size(); ++eid) {
      if (edge_active[eid] && edges[eid].order == 0) {
        edges[eid].order = std::min(max_order, 2);
      }
    }

    // 7. Update node_order: Each node's order is the minimum order among its active connected edges
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (!node_active[i]) continue;
      int min_ord = 999;
      for (int eid : nodes[i].edge_ids) {
        if (edge_active[eid] && edges[eid].order > 0) {
          if (edges[eid].order < min_ord) {
            min_ord = edges[eid].order;
          }
        }
      }
      if (min_ord != 999) {
        node_order[i] = min_ord;
      } else {
        node_order[i] = 1;
      }
    }
    if (crown_node_id >= 0 && crown_node_id < (int)nodes.size()) {
      node_order[crown_node_id] = 1;
    }
  }

  // Calculate Branching Angles for Lateral & Higher-order Roots
  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (!edge_active[eid]) continue;
    auto& ed = edges[eid];
    if (ed.order > 1) {
      int jnode_id = -1;
      // Identify which end is attached to a lower-order parent root
      bool n1_has_parent = false;
      for (int peid : nodes[ed.node1].edge_ids) {
        if (edge_active[peid] && edges[peid].order < ed.order) {
          n1_has_parent = true;
          break;
        }
      }
      bool n2_has_parent = false;
      for (int peid : nodes[ed.node2].edge_ids) {
        if (edge_active[peid] && edges[peid].order < ed.order) {
          n2_has_parent = true;
          break;
        }
      }

      if (n1_has_parent) {
        jnode_id = ed.node1;
      } else if (n2_has_parent) {
        jnode_id = ed.node2;
      }

      if (jnode_id != -1 && node_active[jnode_id]) {
        ed.insertion_depth = nodes[jnode_id].y * pixel_size;

        Point2D p_start = {nodes[jnode_id].x, nodes[jnode_id].y};
        Point2D p_daughter = p_start;

        if (ed.coords.size() > 5) {
          if (ed.node1 == jnode_id) {
            p_daughter = ed.coords[std::min((size_t)10, ed.coords.size() - 1)];
          } else {
            p_daughter = ed.coords[ed.coords.size() - 1 - std::min((size_t)10, ed.coords.size() - 1)];
          }
        } else if (!ed.coords.empty()) {
          p_daughter = (ed.node1 == jnode_id) ? ed.coords.back() : ed.coords.front();
        }

        double vx = p_daughter.x - p_start.x;
        double vy = p_daughter.y - p_start.y;
        double norm_v = std::sqrt(vx * vx + vy * vy);

        double px = 0.0, py = 1.0;

        for (int peid : nodes[jnode_id].edge_ids) {
          if (edge_active[peid] && edges[peid].order < ed.order) {
            int other = (edges[peid].node1 == jnode_id) ? edges[peid].node2 : edges[peid].node1;
            px = nodes[other].x - nodes[jnode_id].x;
            py = nodes[other].y - nodes[jnode_id].y;
            break;
          }
        }

        double norm_p = std::sqrt(px * px + py * py);
        if (norm_v > 1e-4 && norm_p > 1e-4) {
          double dot = (vx * px + vy * py) / (norm_v * norm_p);
          dot = std::max(-1.0, std::min(1.0, dot));
          double angle_deg = std::acos(dot) * (180.0 / M_PI);
          ed.branching_angle = angle_deg;
        }
      }
    }
  }

  // 10. Global Metrics Computation
  double total_length_mm = 0.0;
  double primary_length_mm = 0.0;
  double secondary_length_mm = 0.0;
  double tertiary_length_mm = 0.0;
  double lateral_length_mm = 0.0;
  int n_primary_segments = 0;
  double total_surface_area_mm2 = 0.0;
  double total_volume_mm3 = 0.0;

  std::vector<double> all_angles;
  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (!edge_active[eid]) continue;
    const auto& ed = edges[eid];
    total_length_mm += ed.length_mm;
    total_surface_area_mm2 += ed.surface_area_mm2;
    total_volume_mm3 += ed.volume_mm3;

    if (ed.order == 1) {
      primary_length_mm += ed.length_mm;
      n_primary_segments++;
    } else {
      lateral_length_mm += ed.length_mm;
      if (ed.order == 2) {
        secondary_length_mm += ed.length_mm;
      } else {
        tertiary_length_mm += ed.length_mm;
      }
      if (!R_IsNA(ed.branching_angle)) {
        all_angles.push_back(ed.branching_angle);
      }
    }
  }

  if (count_sec_axes == 0 && secondary_length_mm > 0.0) count_sec_axes = 1;
  if (count_tert_axes == 0 && tertiary_length_mm > 0.0) count_tert_axes = 1;

  int n_primary = (n_primary_segments > 0 ? 1 : 0);
  int n_secondary = count_sec_axes;
  int n_tertiary = count_tert_axes;
  int n_lateral = n_secondary + n_tertiary;

  double mean_branch_angle = NA_REAL;
  double median_branch_angle = NA_REAL;
  if (!all_angles.empty()) {
    double sum_a = 0.0;
    for (double a : all_angles) sum_a += a;
    mean_branch_angle = sum_a / all_angles.size();
    std::sort(all_angles.begin(), all_angles.end());
    median_branch_angle = all_angles[all_angles.size() / 2];
  }

  // Tips & Forks counts (only active nodes)
  int n_tips = 0;
  int n_forks = 0;
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (!node_active[i]) continue;
    if (nodes[i].type == 1) n_tips++;
    else if (nodes[i].type == 2) n_forks++;
  }

  // Projected Area (total foreground pixels * pixel_size^2)
  double projected_area_mm2 = (double)fg_count * pixel_size * pixel_size;

  // Bounding box & Network dimensions
  double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
  double sum_fx = 0.0, sum_fy = 0.0;
  for (const auto& pt : fg_pts) {
    if (pt.x < min_x) min_x = pt.x;
    if (pt.x > max_x) max_x = pt.x;
    if (pt.y < min_y) min_y = pt.y;
    if (pt.y > max_y) max_y = pt.y;
    sum_fx += pt.x;
    sum_fy += pt.y;
  }

  double network_width_mm = (max_x - min_x + 1.0) * pixel_size;
  double network_depth_mm = (max_y - min_y + 1.0) * pixel_size;
  double width_depth_ratio = (network_depth_mm > 0.0) ? (network_width_mm / network_depth_mm) : 0.0;
  double centroid_x = (sum_fx / fg_count) * pixel_size;
  double centroid_y = (sum_fy / fg_count) * pixel_size;

  // Convex Hull & Solidity
  double hull_area_px = compute_convex_hull_area(fg_pts);
  double convex_hull_area_mm2 = hull_area_px * pixel_size * pixel_size;
  double solidity = (convex_hull_area_mm2 > 0.0) ? (projected_area_mm2 / convex_hull_area_mm2) : 1.0;
  double bushiness = (network_depth_mm > 0.0) ? (network_width_mm / network_depth_mm) : 0.0;

  // Average and median diameter across all skeleton pixels
  std::vector<double> skel_diams;
  skel_diams.reserve(skel_pixels.size());
  double sum_skel_diam = 0.0;
  for (int idx : skel_pixels) {
    double d = 2.0 * dist_map[idx] * pixel_size;
    skel_diams.push_back(d);
    sum_skel_diam += d;
  }
  double avg_diameter_mm = (skel_diams.empty()) ? 0.0 : (sum_skel_diam / skel_diams.size());
  std::sort(skel_diams.begin(), skel_diams.end());
  double median_diameter_mm = (skel_diams.empty()) ? 0.0 : skel_diams[skel_diams.size() / 2];

  // Lateral root density (laterals per cm of primary root)
  double primary_len_cm = primary_length_mm / 10.0;
  double lateral_density_per_cm = (primary_len_cm > 0.0) ? ((double)n_lateral / primary_len_cm) : 0.0;

  // 11. Diameter Classes Breakdown (Quantile-based or Custom Bins)
  std::vector<double> cuts;
  if (diameter_bins.isNotNull()) {
    NumericVector db(diameter_bins.get());
    for (int i = 0; i < db.size(); ++i) cuts.push_back(db[i]);
    std::sort(cuts.begin(), cuts.end());
  } else {
    if (num_classes < 2) num_classes = 2;
    std::vector<double> edge_diams;
    for (size_t eid = 0; eid < edges.size(); ++eid) {
      if (edge_active[eid]) edge_diams.push_back(edges[eid].avg_diameter_mm);
    }
    std::sort(edge_diams.begin(), edge_diams.end());

    if (!edge_diams.empty()) {
      for (int k = 1; k < num_classes; ++k) {
        double p = (double)k / (double)num_classes;
        int idx = (int)std::round(p * (edge_diams.size() - 1));
        cuts.push_back(edge_diams[idx]);
      }
      // If duplicates occurred, fall back to equal-width intervals across [min_d, max_d]
      bool has_duplicates = false;
      for (size_t i = 1; i < cuts.size(); ++i) {
        if (cuts[i] <= cuts[i - 1] + 1e-4) {
          has_duplicates = true;
          break;
        }
      }
      if (has_duplicates) {
        cuts.clear();
        double min_d = edge_diams.front();
        double max_d = edge_diams.back();
        if (max_d <= min_d + 1e-4) max_d = min_d + 1.0;
        double step = (max_d - min_d) / (double)num_classes;
        for (int k = 1; k < num_classes; ++k) {
          cuts.push_back(min_d + k * step);
        }
      }
    }
  }

  int n_classes = (int)cuts.size() + 1;
  std::vector<double> class_length(n_classes, 0.0);
  std::vector<double> class_volume(n_classes, 0.0);
  std::vector<double> class_surf(n_classes, 0.0);
  std::vector<std::string> class_names(n_classes);

  for (int c = 0; c < n_classes; ++c) {
    char buf[64];
    if (cuts.empty()) {
      class_names[c] = "All Diameters";
    } else if (c == 0) {
      std::snprintf(buf, sizeof(buf), "< %.2f mm", cuts[0]);
      class_names[c] = std::string(buf);
    } else if (c == n_classes - 1) {
      std::snprintf(buf, sizeof(buf), ">= %.2f mm", cuts[cuts.size() - 1]);
      class_names[c] = std::string(buf);
    } else {
      std::snprintf(buf, sizeof(buf), "%.2f-%.2f mm", cuts[c - 1], cuts[c]);
      class_names[c] = std::string(buf);
    }
  }

  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (!edge_active[eid]) continue;
    double d = edges[eid].avg_diameter_mm;
    int c_idx = 0;
    while (c_idx < (int)cuts.size() && d >= cuts[c_idx]) {
      c_idx++;
    }
    class_length[c_idx] += edges[eid].length_mm;
    class_volume[c_idx] += edges[eid].volume_mm3;
    class_surf[c_idx] += edges[eid].surface_area_mm2;
  }

  // 12. Depth Profile Slicing (Stratified vertical distribution & Root Crossings)
  if (n_depth_slices < 2) n_depth_slices = 10;
  double slice_height_px = (max_y - min_y + 1.0) / (double)n_depth_slices;

  std::vector<double> slice_start_mm(n_depth_slices);
  std::vector<double> slice_end_mm(n_depth_slices);
  std::vector<double> slice_mid_mm(n_depth_slices);
  std::vector<double> slice_length_mm(n_depth_slices, 0.0);
  std::vector<double> slice_volume_mm3(n_depth_slices, 0.0);
  std::vector<double> slice_area_mm2(n_depth_slices, 0.0);
  std::vector<int> slice_crossings(n_depth_slices, 0);

  // Compute depth slices metrics
  for (int s = 0; s < n_depth_slices; ++s) {
    double y_start = min_y + (double)s * slice_height_px;
    double y_end = y_start + slice_height_px;
    double y_mid = 0.5 * (y_start + y_end);

    slice_start_mm[s] = (y_start - min_y) * pixel_size;
    slice_end_mm[s] = (y_end - min_y) * pixel_size;
    slice_mid_mm[s] = (y_mid - min_y) * pixel_size;

    int mid_y_int = (int)std::round(y_mid);
    if (mid_y_int >= 0 && mid_y_int < h) {
      int crossings = 0;
      bool in_root = false;
      for (int x = 0; x < w; ++x) {
        if (skel[x + mid_y_int * w] > 0) {
          if (!in_root) {
            crossings++;
            in_root = true;
          }
        } else {
          in_root = false;
        }
      }
      slice_crossings[s] = crossings;
    }
  }

  for (int idx : skel_pixels) {
    int y = idx / w;
    int s = (int)((y - min_y) / slice_height_px);
    if (s < 0) s = 0;
    if (s >= n_depth_slices) s = n_depth_slices - 1;

    double d = 2.0 * dist_map[idx] * pixel_size;
    slice_length_mm[s] += pixel_size;
    slice_volume_mm3[s] += 0.25 * M_PI * d * d * pixel_size;
    slice_area_mm2[s] += d * pixel_size;
  }

  // 13. Package Output Tables
  DataFrame summary_df = DataFrame::create(
    Named("total_length") = total_length_mm,
    Named("projected_area") = projected_area_mm2,
    Named("surface_area") = total_surface_area_mm2,
    Named("volume") = total_volume_mm3,
    Named("avg_diameter") = avg_diameter_mm,
    Named("median_diameter") = median_diameter_mm,
    Named("max_width") = network_width_mm,
    Named("max_depth") = network_depth_mm,
    Named("width_depth_ratio") = width_depth_ratio,
    Named("network_area") = convex_hull_area_mm2,
    Named("solidity") = solidity,
    Named("bushiness") = bushiness,
    Named("n_tips") = n_tips,
    Named("n_forks") = n_forks,
    Named("n_primary") = (n_primary > 0 ? 1 : 0),
    Named("primary_length") = primary_length_mm,
    Named("n_secondary") = n_secondary,
    Named("secondary_length") = secondary_length_mm,
    Named("mean_secondary_length") = (n_secondary > 0 ? (secondary_length_mm / n_secondary) : 0.0),
    Named("n_tertiary") = n_tertiary,
    Named("tertiary_length") = tertiary_length_mm,
    Named("mean_tertiary_length") = (n_tertiary > 0 ? (tertiary_length_mm / n_tertiary) : 0.0),
    Named("n_lateral") = n_lateral,
    Named("lateral_length") = lateral_length_mm,
    Named("mean_lateral_length") = (n_lateral > 0 ? (lateral_length_mm / n_lateral) : 0.0),
    Named("lateral_density_per_cm") = lateral_density_per_cm,
    Named("mean_branch_angle") = mean_branch_angle,
    Named("median_branch_angle") = median_branch_angle,
    Named("centroid_x") = centroid_x,
    Named("centroid_y") = centroid_y
  );

  // Renumber active nodes
  std::vector<int> old_to_new_node(nodes.size(), -1);
  int new_node_count = 0;
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (node_active[i]) {
      old_to_new_node[i] = ++new_node_count;
    }
  }

  // Pack active edges / roots
  std::vector<int> edge_id_v, node1_v, node2_v, order_v;
  std::vector<double> len_v, diam_v, med_diam_v, surf_v, vol_v, angle_v, depth_v;

  int final_edge_id = 0;
  std::vector<NumericMatrix> path_matrices;
  for (size_t eid = 0; eid < edges.size(); ++eid) {
    if (!edge_active[eid]) continue;
    int n1 = old_to_new_node[edges[eid].node1];
    int n2 = old_to_new_node[edges[eid].node2];
    if (n1 == -1 || n2 == -1 || n1 == n2) continue;

    edge_id_v.push_back(++final_edge_id);
    node1_v.push_back(n1);
    node2_v.push_back(n2);
    order_v.push_back(edges[eid].order);
    len_v.push_back(edges[eid].length_mm);
    diam_v.push_back(edges[eid].avg_diameter_mm);
    med_diam_v.push_back(edges[eid].median_diameter_mm);
    surf_v.push_back(edges[eid].surface_area_mm2);
    vol_v.push_back(edges[eid].volume_mm3);
    angle_v.push_back(edges[eid].branching_angle);
    depth_v.push_back(edges[eid].insertion_depth);

    NumericMatrix pm(edges[eid].coords.size(), 2);
    for (size_t pi = 0; pi < edges[eid].coords.size(); ++pi) {
      pm(pi, 0) = edges[eid].coords[pi].x + 1.0;
      pm(pi, 1) = edges[eid].coords[pi].y + 1.0;
    }
    path_matrices.push_back(pm);
  }

  DataFrame roots_df = DataFrame::create(
    Named("root_id") = wrap(edge_id_v),
    Named("order") = wrap(order_v),
    Named("length") = wrap(len_v),
    Named("avg_diameter") = wrap(diam_v),
    Named("median_diameter") = wrap(med_diam_v),
    Named("surface_area") = wrap(surf_v),
    Named("volume") = wrap(vol_v),
    Named("branching_angle") = wrap(angle_v),
    Named("insertion_depth") = wrap(depth_v),
    Named("node_start") = wrap(node1_v),
    Named("node_end") = wrap(node2_v)
  );

  // Pack active nodes
  std::vector<int> nid_v, ntype_v, ndeg_v, norder_v;
  std::vector<double> nx_v, ny_v;
  std::vector<std::string> nlabel_v;

  for (size_t i = 0; i < nodes.size(); ++i) {
    if (!node_active[i]) continue;
    int nid = old_to_new_node[i];
    nid_v.push_back(nid);
    ntype_v.push_back(nodes[i].type);
    nx_v.push_back(nodes[i].x + 1.0);
    ny_v.push_back(nodes[i].y + 1.0);
    ndeg_v.push_back((int)nodes[i].edge_ids.size());
    norder_v.push_back(std::max(1, node_order[i]));
    if (nodes[i].type == 3) nlabel_v.push_back("crown");
    else if (nodes[i].type == 2) nlabel_v.push_back("junction");
    else nlabel_v.push_back("tip");
  }

  DataFrame nodes_df = DataFrame::create(
    Named("node_id") = wrap(nid_v),
    Named("type") = wrap(nlabel_v),
    Named("x") = wrap(nx_v),
    Named("y") = wrap(ny_v),
    Named("degree") = wrap(ndeg_v),
    Named("order") = wrap(norder_v)
  );

  // Diameter Classes Table
  DataFrame diam_classes_df = DataFrame::create(
    Named("class") = wrap(class_names),
    Named("length") = wrap(class_length),
    Named("surface_area") = wrap(class_surf),
    Named("volume") = wrap(class_volume)
  );

  // Depth Profile Table
  DataFrame depth_profile_df = DataFrame::create(
    Named("slice") = seq(1, n_depth_slices),
    Named("depth_start") = wrap(slice_start_mm),
    Named("depth_end") = wrap(slice_end_mm),
    Named("depth_mid") = wrap(slice_mid_mm),
    Named("length") = wrap(slice_length_mm),
    Named("area") = wrap(slice_area_mm2),
    Named("volume") = wrap(slice_volume_mm3),
    Named("crossings") = wrap(slice_crossings)
  );

  // Build skeleton and distance map matrices to return
  NumericMatrix skel_mat(w, h);
  NumericMatrix dist_mat(w, h);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int idx = x + y * w;
      skel_mat(x, y) = 0.0;
      dist_mat(x, y) = dist_map[idx] * 2.0 * pixel_size; // store diameter
    }
  }

  for (size_t i = 0; i < edges.size(); ++i) {
    if (!edge_active[i]) continue;
    for (const auto& pt : edges[i].coords) {
      int px = (int)std::round(pt.x);
      int py = (int)std::round(pt.y);
      if (px >= 0 && px < w && py >= 0 && py < h) {
        skel_mat(px, py) = 1.0;
      }
    }
  }

  return List::create(
    Named("summary") = summary_df,
    Named("roots") = roots_df,
    Named("nodes") = nodes_df,
    Named("paths") = wrap(path_matrices),
    Named("diameter_classes") = diam_classes_df,
    Named("depth_profile") = depth_profile_df,
    Named("skeleton") = skel_mat,
    Named("diameter_map") = dist_mat,
    Named("crown_node") = (crown_node_id >= 0 && crown_node_id < (int)old_to_new_node.size()) ? old_to_new_node[crown_node_id] : 1,
    Named("pixel_size") = pixel_size
  );
}

