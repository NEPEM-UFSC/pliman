// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <bitset>
#include <cstdint>
#include "stb_image.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static inline double fast_atof(const char* p, const char** next_ptr) {
  while (*p == ' ' || *p == '\t') p++;
  bool neg = false;
  if (*p == '-') { neg = true; p++; }
  else if (*p == '+') { p++; }

  double val = 0.0;
  while (*p >= '0' && *p <= '9') {
    val = val * 10.0 + (*p - '0');
    p++;
  }
  if (*p == '.') {
    p++;
    double frac = 0.1;
    while (*p >= '0' && *p <= '9') {
      val += (*p - '0') * frac;
      frac *= 0.1;
      p++;
    }
  }
  if (*p == 'e' || *p == 'E') {
    p++;
    bool exp_neg = false;
    if (*p == '-') { exp_neg = true; p++; }
    else if (*p == '+') { p++; }
    int exp_val = 0;
    while (*p >= '0' && *p <= '9') {
      exp_val = exp_val * 10 + (*p - '0');
      p++;
    }
    double factor = std::pow(10.0, exp_neg ? -exp_val : exp_val);
    val *= factor;
  }
  if (next_ptr) *next_ptr = p;
  return neg ? -val : val;
}

// ==============================================================================
// HIGH-PERFORMANCE PARALLEL YOLO LABEL PARSING & GEOMETRIC AUDITING (OpenMP)
// ==============================================================================

struct AnomalyItemCpp {
  std::string stem;
  std::string split;
  std::string img;
  std::string lbl;
  std::string category;
  std::string type;
  std::string type_label;
  std::string badge_class;
  std::string detail;
  double xc = 0.0;
  double yc = 0.0;
  double bw = 0.0;
  double bh = 0.0;
};

struct BoxCandidateCpp {
  int cid = 0;
  double x1 = 0.0;
  double y1 = 0.0;
  double x2 = 0.0;
  double y2 = 0.0;
  double xc = 0.0;
  double yc = 0.0;
  double bw = 0.0;
  double bh = 0.0;
};

struct ThreadStatsCpp {
  int total_boxes = 0;
  int total_polys = 0;
  int n_invalid_cids = 0;
  int n_oob = 0;
  int n_zero_area = 0;
  int n_tiny = 0;
  int n_huge = 0;
  int n_extreme_ar = 0;
  int n_duplicate_pairs = 0;

  int n_recorded_tiny = 0;
  int n_recorded_huge = 0;
  int n_recorded_ar = 0;

  std::vector<int> class_split_objs; // [num_classes * num_splits]
  std::vector<int> class_img_counts; // [num_classes]
  std::vector<int> spatial_grid;     // [144] (12x12)
  std::vector<int> dim_grid;         // [100] (10x10)
  std::vector<int> cooc;             // [num_classes * num_classes]

  std::vector<double> all_areas;
  std::vector<double> all_ars;
  std::vector<std::vector<double>> class_box_areas; // [num_classes]
  std::vector<AnomalyItemCpp> anomalies;
};

// [[Rcpp::export]]
Rcpp::List parse_yolo_labels_cpp(
    Rcpp::CharacterVector lbl_files,
    Rcpp::CharacterVector img_files,
    Rcpp::CharacterVector stems,
    Rcpp::CharacterVector splits,
    Rcpp::CharacterVector ptypes,
    Rcpp::CharacterVector classes,
    Rcpp::CharacterVector split_keys,
    int max_anomalies = 500,
    int num_threads = 0
) {
  int n_files = lbl_files.size();
  int num_classes = classes.size();
  int num_splits = split_keys.size();

  // Convert inputs to std::string vectors for thread-safe access
  std::vector<std::string> v_lbl(n_files), v_img(n_files), v_stem(n_files), v_split(n_files), v_ptype(n_files);
  for (int i = 0; i < n_files; ++i) {
    v_lbl[i] = Rcpp::as<std::string>(lbl_files[i]);
    v_img[i] = Rcpp::as<std::string>(img_files[i]);
    v_stem[i] = Rcpp::as<std::string>(stems[i]);
    v_split[i] = Rcpp::as<std::string>(splits[i]);
    v_ptype[i] = Rcpp::as<std::string>(ptypes[i]);
  }

  std::vector<std::string> v_classes(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    v_classes[c] = Rcpp::as<std::string>(classes[c]);
  }

  std::vector<std::string> v_split_keys(num_splits);
  for (int s = 0; s < num_splits; ++s) {
    v_split_keys[s] = Rcpp::as<std::string>(split_keys[s]);
  }

  // Map split string to split index
  std::vector<int> file_split_idx(n_files, 0);
  for (int i = 0; i < n_files; ++i) {
    for (int s = 0; s < num_splits; ++s) {
      if (v_split[i] == v_split_keys[s]) {
        file_split_idx[i] = s;
        break;
      }
    }
  }

  // Determine thread count
  int n_threads_use = 1;
#ifdef _OPENMP
  if (num_threads > 0) {
    n_threads_use = num_threads;
  } else {
    n_threads_use = std::max(1, omp_get_max_threads());
  }
#endif

  std::vector<ThreadStatsCpp> t_stats(n_threads_use);
  for (int t = 0; t < n_threads_use; ++t) {
    t_stats[t].class_split_objs.assign(num_classes * num_splits, 0);
    t_stats[t].class_img_counts.assign(num_classes, 0);
    t_stats[t].spatial_grid.assign(144, 0);
    t_stats[t].dim_grid.assign(100, 0);
    t_stats[t].cooc.assign(num_classes * num_classes, 0);
    t_stats[t].class_box_areas.resize(num_classes);
  }

  std::vector<int> file_n_objs(n_files, 0);
  std::vector<bool> file_is_bg(n_files, true);

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads_use) if(n_files > 8)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    ThreadStatsCpp& ts = t_stats[tid];

    std::vector<double> line_vals;
    line_vals.reserve(32);
    std::vector<BoxCandidateCpp> img_boxes;
    img_boxes.reserve(64);
    std::vector<int> classes_in_img;
    classes_in_img.reserve(16);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 32)
#endif
    for (int i = 0; i < n_files; ++i) {
      const std::string& lpath = v_lbl[i];
      if (lpath.empty()) {
        file_is_bg[i] = true;
        file_n_objs[i] = 0;
        continue;
      }

      FILE* fp = std::fopen(lpath.c_str(), "rb");
      if (!fp) {
        file_is_bg[i] = true;
        file_n_objs[i] = 0;
        continue;
      }

      std::fseek(fp, 0, SEEK_END);
      long fsize = std::ftell(fp);
      std::fseek(fp, 0, SEEK_SET);

      if (fsize <= 0) {
        std::fclose(fp);
        file_is_bg[i] = true;
        file_n_objs[i] = 0;
        continue;
      }

      std::string buf(fsize, '\0');
      size_t read_bytes = std::fread(&buf[0], 1, fsize, fp);
      std::fclose(fp);
      buf.resize(read_bytes);

      if (buf.empty()) {
        file_is_bg[i] = true;
        file_n_objs[i] = 0;
        continue;
      }

      const char* p = buf.c_str();
      const char* end = p + buf.size();

      int sp_idx = file_split_idx[i];
      int n_lines_in_file = 0;
      img_boxes.clear();
      classes_in_img.clear();

      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) p++;
        if (p >= end) break;

        const char* line_end = p;
        while (line_end < end && *line_end != '\r' && *line_end != '\n') line_end++;

        line_vals.clear();
        const char* cur = p;
        while (cur < line_end) {
          while (cur < line_end && (*cur == ' ' || *cur == '\t' || *cur == ',')) cur++;
          if (cur >= line_end) break;
          const char* next_ptr = nullptr;
          double v = fast_atof(cur, &next_ptr);
          if (next_ptr == cur) {
            while (cur < line_end && *cur != ' ' && *cur != '\t' && *cur != ',') cur++;
          } else {
            line_vals.push_back(v);
            cur = next_ptr;
          }
        }
        p = line_end;

        if (line_vals.empty()) continue;
        n_lines_in_file++;

        int raw_cid = static_cast<int>(line_vals[0]);
        bool valid_cid = (raw_cid >= 0 && raw_cid < num_classes);

        if (!valid_cid) {
          ts.n_invalid_cids++;
          if ((int)ts.anomalies.size() < max_anomalies) {
            AnomalyItemCpp an;
            an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
            an.category = v_ptype[i]; an.type = "invalid_class_id";
            an.type_label = "Invalid Class ID"; an.badge_class = "pill-danger";
            char detail_buf[128];
            std::snprintf(detail_buf, sizeof(detail_buf), "Class ID %d exceeds defined classes [0, %d]", raw_cid, num_classes - 1);
            an.detail = detail_buf;
            ts.anomalies.push_back(an);
          }
        } else {
          ts.class_split_objs[raw_cid * num_splits + sp_idx]++;
          classes_in_img.push_back(raw_cid);
        }

        // 1. BOUNDING BOX FORMAT (5 values: cid xc yc bw bh)
        if (line_vals.size() == 5) {
          ts.total_boxes++;
          double xc = line_vals[1];
          double yc = line_vals[2];
          double bw = line_vals[3];
          double bh = line_vals[4];

          if (std::isnan(xc) || std::isnan(yc) || std::isnan(bw) || std::isnan(bh)) continue;

          bool is_oob = (xc < 0.0 || xc > 1.0 || yc < 0.0 || yc > 1.0 || bw < 0.0 || bw > 1.0 || bh < 0.0 || bh > 1.0);
          if (is_oob) {
            ts.n_oob++;
            if ((int)ts.anomalies.size() < max_anomalies) {
              AnomalyItemCpp an;
              an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
              an.category = v_ptype[i]; an.type = "out_of_bounds";
              an.type_label = "Out of Bounds"; an.badge_class = "pill-danger";
              char detail_buf[128];
              std::snprintf(detail_buf, sizeof(detail_buf), "Box [xc=%.2f, yc=%.2f, w=%.2f, h=%.2f] exceeds [0, 1] range", xc, yc, bw, bh);
              an.detail = detail_buf;
              an.xc = xc; an.yc = yc; an.bw = bw; an.bh = bh;
              ts.anomalies.push_back(an);
            }
          }

          if (bw <= 0.0 || bh <= 0.0) {
            ts.n_zero_area++;
            if ((int)ts.anomalies.size() < max_anomalies) {
              AnomalyItemCpp an;
              an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
              an.category = v_ptype[i]; an.type = "zero_area";
              an.type_label = "Zero Area"; an.badge_class = "pill-danger";
              char detail_buf[128];
              std::snprintf(detail_buf, sizeof(detail_buf), "Degenerate box with width=%.3f, height=%.3f", bw, bh);
              an.detail = detail_buf;
              an.xc = xc; an.yc = yc; an.bw = bw; an.bh = bh;
              ts.anomalies.push_back(an);
            }
          } else {
            double area_pct = bw * bh * 100.0;
            if (area_pct < 0.05) {
              ts.n_tiny++;
              if (ts.n_recorded_tiny < 60 && (int)ts.anomalies.size() < max_anomalies) {
                ts.n_recorded_tiny++;
                AnomalyItemCpp an;
                an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
                an.category = v_ptype[i]; an.type = "micro_box";
                an.type_label = "Micro-Box"; an.badge_class = "pill-warn";
                char detail_buf[128];
                std::snprintf(detail_buf, sizeof(detail_buf), "Micro-box area = %.3f%% of canvas (w=%.3f, h=%.3f)", area_pct, bw, bh);
                an.detail = detail_buf;
                an.xc = xc; an.yc = yc; an.bw = bw; an.bh = bh;
                ts.anomalies.push_back(an);
              }
            }
            if (area_pct > 85.0) {
              ts.n_huge++;
              if (ts.n_recorded_huge < 50 && (int)ts.anomalies.size() < max_anomalies) {
                ts.n_recorded_huge++;
                AnomalyItemCpp an;
                an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
                an.category = v_ptype[i]; an.type = "huge_box";
                an.type_label = "Huge Box"; an.badge_class = "pill-warn";
                char detail_buf[128];
                std::snprintf(detail_buf, sizeof(detail_buf), "Huge box occupies %.1f%% of canvas (w=%.3f, h=%.3f)", area_pct, bw, bh);
                an.detail = detail_buf;
                an.xc = xc; an.yc = yc; an.bw = bw; an.bh = bh;
                ts.anomalies.push_back(an);
              }
            }

            double ar = bw / std::max(1e-6, bh);
            if (ar > 8.0 || ar < 0.125) {
              ts.n_extreme_ar++;
              if (ts.n_recorded_ar < 80 && (int)ts.anomalies.size() < max_anomalies) {
                ts.n_recorded_ar++;
                AnomalyItemCpp an;
                an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
                an.category = v_ptype[i]; an.type = "extreme_aspect_ratio";
                an.type_label = "Extreme AR"; an.badge_class = "pill-info";
                char detail_buf[128];
                std::snprintf(detail_buf, sizeof(detail_buf), "Extreme aspect ratio %.1f (%s: w=%.3f, h=%.3f)",
                              ar, (ar > 1.0 ? "Landscape" : "Portrait"), bw, bh);
                an.detail = detail_buf;
                an.xc = xc; an.yc = yc; an.bw = bw; an.bh = bh;
                ts.anomalies.push_back(an);
              }
            }

            ts.all_areas.push_back(area_pct);
            ts.all_ars.push_back(ar);
            if (valid_cid) {
              ts.class_box_areas[raw_cid].push_back(area_pct);
            }

            int gx = std::min(11, std::max(0, static_cast<int>(std::floor(xc * 12.0))));
            int gy = std::min(11, std::max(0, static_cast<int>(std::floor(yc * 12.0))));
            ts.spatial_grid[gy * 12 + gx]++;

            int gw = std::min(9, std::max(0, static_cast<int>(std::floor(bw * 10.0))));
            int gh = std::min(9, std::max(0, static_cast<int>(std::floor(bh * 10.0))));
            ts.dim_grid[gh * 10 + gw]++;

            BoxCandidateCpp bc;
            bc.cid = raw_cid;
            bc.x1 = xc - bw * 0.5; bc.y1 = yc - bh * 0.5;
            bc.x2 = xc + bw * 0.5; bc.y2 = yc + bh * 0.5;
            bc.xc = xc; bc.yc = yc; bc.bw = bw; bc.bh = bh;
            img_boxes.push_back(bc);
          }

        // 2. POLYGON FORMAT (> 5 values: cid x1 y1 x2 y2 ...)
        } else if (line_vals.size() > 5) {
          ts.total_polys++;
          double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
          double sum_x = 0.0, sum_y = 0.0;
          int n_pts = 0;
          bool oob = false;

          for (size_t v = 1; v + 1 < line_vals.size(); v += 2) {
            double vx = line_vals[v];
            double vy = line_vals[v + 1];
            if (vx < 0.0 || vx > 1.0 || vy < 0.0 || vy > 1.0) oob = true;
            if (vx < min_x) min_x = vx;
            if (vx > max_x) max_x = vx;
            if (vy < min_y) min_y = vy;
            if (vy > max_y) max_y = vy;
            sum_x += vx; sum_y += vy;
            n_pts++;
          }

          if (oob) {
            ts.n_oob++;
            if ((int)ts.anomalies.size() < max_anomalies) {
              AnomalyItemCpp an;
              an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
              an.category = v_ptype[i]; an.type = "out_of_bounds";
              an.type_label = "Out of Bounds"; an.badge_class = "pill-danger";
              an.detail = "Polygon has vertices outside [0, 1] range";
              an.xc = (min_x + max_x) * 0.5; an.yc = (min_y + max_y) * 0.5;
              an.bw = max_x - min_x; an.bh = max_y - min_y;
              ts.anomalies.push_back(an);
            }
          }

          if (n_pts > 0) {
            double bw = max_x - min_x;
            double bh = max_y - min_y;
            double xc = sum_x / n_pts;
            double yc = sum_y / n_pts;

            if (bw <= 0.0 || bh <= 0.0) {
              ts.n_zero_area++;
            } else {
              double area_pct = bw * bh * 100.0;
              if (area_pct < 0.05) ts.n_tiny++;
              if (area_pct > 85.0) ts.n_huge++;
              double ar = bw / std::max(1e-6, bh);
              if (ar > 8.0 || ar < 0.125) ts.n_extreme_ar++;

              ts.all_areas.push_back(area_pct);
              ts.all_ars.push_back(ar);
              if (valid_cid) {
                ts.class_box_areas[raw_cid].push_back(area_pct);
              }

              int gx = std::min(11, std::max(0, static_cast<int>(std::floor(xc * 12.0))));
              int gy = std::min(11, std::max(0, static_cast<int>(std::floor(yc * 12.0))));
              ts.spatial_grid[gy * 12 + gx]++;

              int gw = std::min(9, std::max(0, static_cast<int>(std::floor(bw * 10.0))));
              int gh = std::min(9, std::max(0, static_cast<int>(std::floor(bh * 10.0))));
              ts.dim_grid[gh * 10 + gw]++;

              BoxCandidateCpp bc;
              bc.cid = raw_cid;
              bc.x1 = min_x; bc.y1 = min_y;
              bc.x2 = max_x; bc.y2 = max_y;
              bc.xc = xc; bc.yc = yc; bc.bw = bw; bc.bh = bh;
              img_boxes.push_back(bc);
            }
          }
        }
      }

      file_n_objs[i] = n_lines_in_file;
      file_is_bg[i] = (n_lines_in_file == 0);

      // Fast pairwise duplicate box check within same image using sweep-and-prune
      int nb = (int)img_boxes.size();
      if (nb >= 2) {
        std::sort(img_boxes.begin(), img_boxes.end(), [](const BoxCandidateCpp& a, const BoxCandidateCpp& b) {
          return a.x1 < b.x1;
        });
        for (int b1 = 0; b1 < nb - 1; ++b1) {
          const auto& box1 = img_boxes[b1];
          for (int b2 = b1 + 1; b2 < nb; ++b2) {
            const auto& box2 = img_boxes[b2];
            if (box2.x1 > box1.x2) break; // Impossible to overlap beyond right boundary!
            if (box1.cid == box2.cid && box2.y1 < box1.y2 && box2.y2 > box1.y1) {
              double ix1 = std::max(box1.x1, box2.x1);
              double iy1 = std::max(box1.y1, box2.y1);
              double ix2 = std::min(box1.x2, box2.x2);
              double iy2 = std::min(box1.y2, box2.y2);
              if (ix2 > ix1 && iy2 > iy1) {
                double iarea = (ix2 - ix1) * (iy2 - iy1);
                double area1 = (box1.x2 - box1.x1) * (box1.y2 - box1.y1);
                double area2 = (box2.x2 - box2.x1) * (box2.y2 - box2.y1);
                double iou = iarea / (area1 + area2 - iarea);
                if (iou > 0.85) {
                  ts.n_duplicate_pairs++;
                  if ((int)ts.anomalies.size() < max_anomalies) {
                    AnomalyItemCpp an;
                    an.stem = v_stem[i]; an.split = v_split[i]; an.img = v_img[i]; an.lbl = lpath;
                    an.category = v_ptype[i]; an.type = "duplicate_pair";
                    an.type_label = "Duplicate Pair"; an.badge_class = "pill-danger";
                    char detail_buf[128];
                    std::string cname_str = (box1.cid >= 0 && box1.cid < num_classes) ? v_classes[box1.cid] : "class";
                    std::snprintf(detail_buf, sizeof(detail_buf), "Boxes %d and %d share high IoU = %.2f (class: %s)",
                                  b1 + 1, b2 + 1, iou, cname_str.c_str());
                    an.detail = detail_buf;
                    an.xc = box1.xc; an.yc = box1.yc; an.bw = box1.bw; an.bh = box1.bh;
                    ts.anomalies.push_back(an);
                  }
                }
              }
            }
          }
        }
      }

      // Track class presence and co-occurrence
      if (!classes_in_img.empty()) {
        std::sort(classes_in_img.begin(), classes_in_img.end());
        classes_in_img.erase(std::unique(classes_in_img.begin(), classes_in_img.end()), classes_in_img.end());
        for (int cid_val : classes_in_img) {
          if (cid_val >= 0 && cid_val < num_classes) {
            ts.class_img_counts[cid_val]++;
          }
        }
        if (classes_in_img.size() >= 2 && num_classes > 1) {
          for (size_t c1 = 0; c1 < classes_in_img.size(); ++c1) {
            for (size_t c2 = 0; c2 < classes_in_img.size(); ++c2) {
              if (c1 != c2) {
                int id1 = classes_in_img[c1];
                int id2 = classes_in_img[c2];
                if (id1 >= 0 && id1 < num_classes && id2 >= 0 && id2 < num_classes) {
                  ts.cooc[id1 * num_classes + id2]++;
                }
              }
            }
          }
        }
      }
    } // for i in 0..n_files
  } // #pragma omp parallel

  // Master accumulation of thread results
  int total_boxes = 0;
  int total_polys = 0;
  int n_invalid_cids = 0;
  int n_oob = 0;
  int n_zero_area = 0;
  int n_tiny = 0;
  int n_huge = 0;
  int n_extreme_ar = 0;
  int n_duplicate_pairs = 0;

  std::vector<int> master_class_split_objs(num_classes * num_splits, 0);
  std::vector<int> master_class_img_counts(num_classes, 0);
  std::vector<int> master_spatial_grid(144, 0);
  std::vector<int> master_dim_grid(100, 0);
  std::vector<int> master_cooc(num_classes * num_classes, 0);

  size_t total_areas_count = 0;
  for (int t = 0; t < n_threads_use; ++t) {
    total_areas_count += t_stats[t].all_areas.size();
  }

  std::vector<double> master_all_areas;
  master_all_areas.reserve(total_areas_count);
  std::vector<double> master_all_ars;
  master_all_ars.reserve(total_areas_count);

  std::vector<std::vector<double>> master_class_box_areas(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    size_t sz_c = 0;
    for (int t = 0; t < n_threads_use; ++t) {
      sz_c += t_stats[t].class_box_areas[c].size();
    }
    master_class_box_areas[c].reserve(sz_c);
  }

  std::vector<AnomalyItemCpp> master_anomalies;
  master_anomalies.reserve(max_anomalies);

  for (int t = 0; t < n_threads_use; ++t) {
    const ThreadStatsCpp& ts = t_stats[t];
    total_boxes += ts.total_boxes;
    total_polys += ts.total_polys;
    n_invalid_cids += ts.n_invalid_cids;
    n_oob += ts.n_oob;
    n_zero_area += ts.n_zero_area;
    n_tiny += ts.n_tiny;
    n_huge += ts.n_huge;
    n_extreme_ar += ts.n_extreme_ar;
    n_duplicate_pairs += ts.n_duplicate_pairs;

    for (size_t k = 0; k < master_class_split_objs.size(); ++k) {
      master_class_split_objs[k] += ts.class_split_objs[k];
    }
    for (size_t k = 0; k < master_class_img_counts.size(); ++k) {
      master_class_img_counts[k] += ts.class_img_counts[k];
    }
    for (size_t k = 0; k < 144; ++k) {
      master_spatial_grid[k] += ts.spatial_grid[k];
    }
    for (size_t k = 0; k < 100; ++k) {
      master_dim_grid[k] += ts.dim_grid[k];
    }
    for (size_t k = 0; k < master_cooc.size(); ++k) {
      master_cooc[k] += ts.cooc[k];
    }

    master_all_areas.insert(master_all_areas.end(), ts.all_areas.begin(), ts.all_areas.end());
    master_all_ars.insert(master_all_ars.end(), ts.all_ars.begin(), ts.all_ars.end());

    for (int c = 0; c < num_classes; ++c) {
      master_class_box_areas[c].insert(master_class_box_areas[c].end(),
                                      ts.class_box_areas[c].begin(), ts.class_box_areas[c].end());
    }

    for (const auto& an : ts.anomalies) {
      if ((int)master_anomalies.size() < max_anomalies) {
        master_anomalies.push_back(an);
      }
    }
  }

  // Construct R output structures
  Rcpp::IntegerMatrix mat_class_split(num_classes, num_splits);
  for (int c = 0; c < num_classes; ++c) {
    for (int s = 0; s < num_splits; ++s) {
      mat_class_split(c, s) = master_class_split_objs[c * num_splits + s];
    }
  }

  Rcpp::IntegerVector vec_class_imgs(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    vec_class_imgs[c] = master_class_img_counts[c];
  }

  Rcpp::List list_class_areas(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    list_class_areas[c] = Rcpp::wrap(master_class_box_areas[c]);
  }

  Rcpp::IntegerMatrix mat_spatial(12, 12);
  for (int r = 0; r < 12; ++r) {
    for (int col = 0; col < 12; ++col) {
      mat_spatial(r, col) = master_spatial_grid[r * 12 + col];
    }
  }

  Rcpp::IntegerMatrix mat_dim(10, 10);
  for (int r = 0; r < 10; ++r) {
    for (int col = 0; col < 10; ++col) {
      mat_dim(r, col) = master_dim_grid[r * 10 + col];
    }
  }

  Rcpp::IntegerMatrix mat_cooc(num_classes, num_classes);
  for (int r = 0; r < num_classes; ++r) {
    for (int col = 0; col < num_classes; ++col) {
      mat_cooc(r, col) = master_cooc[r * num_classes + col];
    }
  }

  Rcpp::List anom_list(master_anomalies.size());
  for (size_t a = 0; a < master_anomalies.size(); ++a) {
    const auto& an = master_anomalies[a];
    Rcpp::NumericVector coords = Rcpp::NumericVector::create(an.xc, an.yc, an.bw, an.bh);
    anom_list[a] = Rcpp::List::create(
      Rcpp::Named("stem") = an.stem,
      Rcpp::Named("split") = an.split,
      Rcpp::Named("img") = an.img,
      Rcpp::Named("lbl") = an.lbl,
      Rcpp::Named("category") = an.category,
      Rcpp::Named("type") = an.type,
      Rcpp::Named("type_label") = an.type_label,
      Rcpp::Named("badge_class") = an.badge_class,
      Rcpp::Named("detail") = an.detail,
      Rcpp::Named("coords") = coords
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("class_split_objs") = mat_class_split,
    Rcpp::Named("class_img_counts") = vec_class_imgs,
    Rcpp::Named("class_box_areas") = list_class_areas,
    Rcpp::Named("all_box_areas") = Rcpp::wrap(master_all_areas),
    Rcpp::Named("all_box_ars") = Rcpp::wrap(master_all_ars),
    Rcpp::Named("file_n_objs") = Rcpp::wrap(file_n_objs),
    Rcpp::Named("file_is_bg") = Rcpp::wrap(file_is_bg),
    Rcpp::Named("n_boxes_total") = total_boxes,
    Rcpp::Named("n_polys_total") = total_polys,
    Rcpp::Named("spatial_grid") = mat_spatial,
    Rcpp::Named("dim_grid") = mat_dim,
    Rcpp::Named("class_cooccurrence") = mat_cooc,
    Rcpp::Named("n_invalid_cids") = n_invalid_cids,
    Rcpp::Named("n_oob") = n_oob,
    Rcpp::Named("n_zero_area") = n_zero_area,
    Rcpp::Named("n_tiny") = n_tiny,
    Rcpp::Named("n_huge") = n_huge,
    Rcpp::Named("n_extreme_ar") = n_extreme_ar,
    Rcpp::Named("n_duplicate_pairs") = n_duplicate_pairs,
    Rcpp::Named("anomalies_records") = anom_list
  );
}

// ==============================================================================
// HIGH-PERFORMANCE PARALLEL DATASET MERGE & LABEL REMAP (OpenMP)
// ==============================================================================

// [[Rcpp::export]]
Rcpp::LogicalVector parallel_copy_files_cpp(
    Rcpp::CharacterVector src_files,
    Rcpp::CharacterVector dest_files,
    bool copy = true,
    int num_threads = 0
) {
  int n = src_files.size();
  if (dest_files.size() != n) {
    Rcpp::stop("src_files and dest_files must have the same length.");
  }

  Rcpp::LogicalVector res(n, true);

  std::vector<std::string> v_src(n), v_dst(n);
  for (int i = 0; i < n; ++i) {
    v_src[i] = Rcpp::as<std::string>(src_files[i]);
    v_dst[i] = Rcpp::as<std::string>(dest_files[i]);
  }

  int threads = (num_threads > 0) ? num_threads : 0;
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_max_threads();
#else
  threads = 1;
#endif

#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
  for (int i = 0; i < n; ++i) {
    const std::string& src = v_src[i];
    const std::string& dst = v_dst[i];
    std::error_code ec;

    if (!copy) {
      // Move / rename
      std::filesystem::rename(src, dst, ec);
      if (ec) {
        // Fallback to copy + remove (e.g. across drives or filesystem permissions)
        ec.clear();
        std::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing, ec);
        if (!ec) {
          std::filesystem::remove(src, ec);
        }
      }
      res[i] = !ec;
    } else {
      // Copy
      std::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing, ec);
      res[i] = !ec;
    }
  }

  return res;
}

// [[Rcpp::export]]
int parallel_remap_labels_cpp(
    Rcpp::CharacterVector src_files,
    Rcpp::CharacterVector dest_files,
    Rcpp::IntegerVector class_map,
    int num_threads = 0
) {
  int n = src_files.size();
  if (dest_files.size() != n) {
    Rcpp::stop("src_files and dest_files must have the same length.");
  }

  std::vector<std::string> v_src(n), v_dst(n);
  for (int i = 0; i < n; ++i) {
    v_src[i] = Rcpp::as<std::string>(src_files[i]);
    v_dst[i] = Rcpp::as<std::string>(dest_files[i]);
  }

  std::vector<int> v_map(class_map.size());
  bool is_identity = true;
  for (int i = 0; i < class_map.size(); ++i) {
    v_map[i] = class_map[i];
    if (v_map[i] != i) is_identity = false;
  }

  int threads = (num_threads > 0) ? num_threads : 0;
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_max_threads();
#else
  threads = 1;
#endif

  int total_objects = 0;

#pragma omp parallel for schedule(dynamic, 64) num_threads(threads) reduction(+:total_objects)
  for (int i = 0; i < n; ++i) {
    const std::string& src = v_src[i];
    const std::string& dst = v_dst[i];

    FILE* fin = fopen(src.c_str(), "rb");
    if (!fin) {
      // Create empty destination file
      FILE* fout = fopen(dst.c_str(), "wb");
      if (fout) fclose(fout);
      continue;
    }

    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    if (fsize <= 0) {
      fclose(fin);
      FILE* fout = fopen(dst.c_str(), "wb");
      if (fout) fclose(fout);
      continue;
    }
    fseek(fin, 0, SEEK_SET);

    if (is_identity) {
      // Fast path: Identity class mapping - direct copy
      fclose(fin);
      std::error_code ec;
      std::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing, ec);
      continue;
    }

    std::string content;
    content.resize(fsize);
    size_t bytes_read = fread(&content[0], 1, fsize, fin);
    fclose(fin);
    content.resize(bytes_read);

    std::string out_content;
    out_content.reserve(content.size() + 64);

    size_t pos = 0;
    size_t len = content.size();

    while (pos < len) {
      while (pos < len && (content[pos] == ' ' || content[pos] == '\t' || content[pos] == '\r' || content[pos] == '\n')) {
        pos++;
      }
      if (pos >= len) break;

      size_t line_end = pos;
      while (line_end < len && content[line_end] != '\r' && content[line_end] != '\n') {
        line_end++;
      }

      size_t token_end = pos;
      while (token_end < line_end && content[token_end] != ' ' && content[token_end] != '\t') {
        token_end++;
      }

      if (token_end > pos) {
        int old_cid = 0;
        bool valid_cid = true;
        try {
          old_cid = std::stoi(content.substr(pos, token_end - pos));
        } catch (...) {
          valid_cid = false;
        }

        if (valid_cid) {
          int new_cid = (old_cid >= 0 && old_cid < (int)v_map.size()) ? v_map[old_cid] : old_cid;
          out_content += std::to_string(new_cid);
          out_content.append(content, token_end, line_end - token_end);
          out_content += '\n';
          total_objects++;
        }
      }

      pos = line_end;
    }

    FILE* fout = fopen(dst.c_str(), "wb");
    if (fout) {
      if (!out_content.empty()) {
        fwrite(out_content.data(), 1, out_content.size(), fout);
      }
      fclose(fout);
    }
  }

  return total_objects;
}

// [[Rcpp::export]]
double fast_dir_size_cpp(std::string dir_path) {
  std::error_code ec;
  std::filesystem::path p(dir_path);
  if (!std::filesystem::exists(p, ec) || !std::filesystem::is_directory(p, ec)) {
    return 0.0;
  }
  double total = 0.0;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(p, std::filesystem::directory_options::skip_permission_denied, ec)) {
    if (ec) break;
    if (entry.is_regular_file(ec)) {
      total += static_cast<double>(entry.file_size(ec));
    }
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix fast_image_dimensions_cpp(Rcpp::CharacterVector img_files) {
  int n = img_files.size();
  std::vector<std::string> paths(n);
  for (int i = 0; i < n; ++i) {
    paths[i] = Rcpp::as<std::string>(img_files[i]);
  }
  std::vector<int> widths(n, 0);
  std::vector<int> heights(n, 0);

#if defined(_OPENMP)
  #pragma omp parallel for if(n > 32) schedule(dynamic, 16)
#endif
  for (int i = 0; i < n; ++i) {
    int w = 0, h = 0, comp = 0;
    int ok = stbi_info(paths[i].c_str(), &w, &h, &comp);
    if (ok) {
      widths[i] = w;
      heights[i] = h;
    }
  }

  Rcpp::IntegerMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    out(i, 0) = widths[i];
    out(i, 1) = heights[i];
  }
  return out;
}

// ==============================================================================
// HIGH-PERFORMANCE IMAGE SCRAMBLING / PUZZLE PERMUTATION (OpenMP)
// ==============================================================================

// [[Rcpp::export]]
SEXP image_scramble_cpp(
    SEXP img_sexp,
    int grid_r = 8,
    int grid_c = 8,
    bool rotate_tiles = true,
    bool flip_tiles = true,
    std::string mode = "grid",
    double swap_prob = 0.5,
    bool shuffle_pixels = false,
    int seed = 42
) {
  RawVector raw_vec;
  bool is_raw = (TYPEOF(img_sexp) == RAWSXP);
  NumericVector num_vec;

  if (is_raw) {
    raw_vec = img_sexp;
  } else {
    num_vec = img_sexp;
  }

  SEXP dim_attr = Rf_getAttrib(img_sexp, R_DimSymbol);
  if (Rf_isNull(dim_attr)) {
    Rcpp::stop("Image must have a 'dim' attribute.");
  }
  IntegerVector dims(dim_attr);
  if (dims.size() < 2) {
    Rcpp::stop("Image dimensions must be at least 2D (w x h).");
  }

  int w = dims[0];
  int h = dims[1];
  int channels = (dims.size() >= 3) ? dims[2] : 1;

  std::size_t total_elements = static_cast<std::size_t>(w) * h * channels;
  std::size_t plane_size = static_cast<std::size_t>(w) * h;

  uint64_t rng_state = (seed > 0) ? static_cast<uint64_t>(seed) : 123456789ULL;
  auto fast_rand = [&rng_state]() -> uint32_t {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<uint32_t>(rng_state >> 32);
  };
  auto fast_unif = [&fast_rand]() -> double {
    return static_cast<double>(fast_rand()) / 4294967296.0;
  };

  // ============================================================================
  // MODE: "pixel" - Complete global spatial shuffle of all pixels
  // ============================================================================
  if (mode == "pixel") {
    std::vector<std::size_t> p_indices(plane_size);
    for (std::size_t i = 0; i < plane_size; ++i) {
      p_indices[i] = i;
    }
    for (std::size_t i = plane_size - 1; i > 0; --i) {
      std::size_t j = static_cast<std::size_t>(fast_rand() % (i + 1));
      std::swap(p_indices[i], p_indices[j]);
    }

    if (is_raw) {
      RawVector out(total_elements);
      const unsigned char* p_in = reinterpret_cast<const unsigned char*>(raw_vec.begin());
      unsigned char* p_out = reinterpret_cast<unsigned char*>(out.begin());

      for (int ch = 0; ch < channels; ++ch) {
        std::size_t ch_offset = static_cast<std::size_t>(ch) * plane_size;
        for (std::size_t i = 0; i < plane_size; ++i) {
          p_out[ch_offset + i] = p_in[ch_offset + p_indices[i]];
        }
      }

      out.attr("dim") = raw_vec.attr("dim");
      if (raw_vec.hasAttribute("class")) out.attr("class") = raw_vec.attr("class");
      if (raw_vec.hasAttribute("colormode")) out.attr("colormode") = raw_vec.attr("colormode");
      return out;
    } else {
      NumericVector out(total_elements);
      const double* p_in = num_vec.begin();
      double* p_out = out.begin();

      for (int ch = 0; ch < channels; ++ch) {
        std::size_t ch_offset = static_cast<std::size_t>(ch) * plane_size;
        for (std::size_t i = 0; i < plane_size; ++i) {
          p_out[ch_offset + i] = p_in[ch_offset + p_indices[i]];
        }
      }

      out.attr("dim") = num_vec.attr("dim");
      if (num_vec.hasAttribute("class")) out.attr("class") = num_vec.attr("class");
      if (num_vec.hasAttribute("colormode")) out.attr("colormode") = num_vec.attr("colormode");
      return out;
    }
  }

  // ============================================================================
  // MODES: "grid" and "block_swap"
  // ============================================================================
  if (grid_r < 1) grid_r = 1;
  if (grid_c < 1) grid_c = 1;
  if (grid_r > h) grid_r = h;
  if (grid_c > w) grid_c = w;

  int num_tiles = grid_r * grid_c;

  // Compute tile bounds
  std::vector<int> col_starts(grid_c);
  std::vector<int> col_ends(grid_c);
  for (int c = 0; c < grid_c; ++c) {
    col_starts[c] = static_cast<int>(c * static_cast<double>(w) / grid_c);
    col_ends[c] = static_cast<int>((c + 1) * static_cast<double>(w) / grid_c);
  }

  std::vector<int> row_starts(grid_r);
  std::vector<int> row_ends(grid_r);
  for (int r = 0; r < grid_r; ++r) {
    row_starts[r] = static_cast<int>(r * static_cast<double>(h) / grid_r);
    row_ends[r] = static_cast<int>((r + 1) * static_cast<double>(h) / grid_r);
  }

  struct TileCoord {
    int r, c;
  };
  std::vector<TileCoord> tiles(num_tiles);
  for (int r = 0; r < grid_r; ++r) {
    for (int c = 0; c < grid_c; ++c) {
      tiles[r * grid_c + c] = {r, c};
    }
  }

  std::vector<TileCoord> permuted = tiles;

  if (mode == "block_swap") {
    double p_swap = (swap_prob > 0.0 && swap_prob <= 1.0) ? swap_prob : 0.5;
    int num_swaps = std::max(1, static_cast<int>(num_tiles * p_swap));

    for (int s = 0; s < num_swaps; ++s) {
      int t1 = fast_rand() % num_tiles;
      int r1 = tiles[t1].r;
      int c1 = tiles[t1].c;

      int dr = (fast_rand() % 3) - 1; // -1, 0, 1
      int dc = (fast_rand() % 3) - 1; // -1, 0, 1
      if (dr == 0 && dc == 0) {
        dr = (fast_rand() % 2 == 0) ? 1 : -1;
      }

      int r2 = std::max(0, std::min(grid_r - 1, r1 + dr));
      int c2 = std::max(0, std::min(grid_c - 1, c1 + dc));
      int t2 = r2 * grid_c + c2;

      std::swap(permuted[t1], permuted[t2]);
    }
  } else {
    // Mode "grid" (default)
    for (int i = num_tiles - 1; i > 0; --i) {
      int j = fast_rand() % (i + 1);
      std::swap(permuted[i], permuted[j]);
    }
  }

  std::vector<int> tile_rot(num_tiles, 0);
  std::vector<bool> tile_flip_x(num_tiles, false);
  std::vector<bool> tile_flip_y(num_tiles, false);

  if (rotate_tiles) {
    for (int i = 0; i < num_tiles; ++i) {
      if (mode == "block_swap" && fast_unif() > 0.5) {
        tile_rot[i] = 0;
      } else {
        tile_rot[i] = fast_rand() % 4;
      }
    }
  }
  if (flip_tiles) {
    for (int i = 0; i < num_tiles; ++i) {
      if (mode == "block_swap" && fast_unif() > 0.5) {
        tile_flip_x[i] = false;
        tile_flip_y[i] = false;
      } else {
        tile_flip_x[i] = (fast_rand() % 2) == 1;
        tile_flip_y[i] = (fast_rand() % 2) == 1;
      }
    }
  }

  if (is_raw) {
    RawVector out(total_elements);
    const unsigned char* p_in = reinterpret_cast<const unsigned char*>(raw_vec.begin());
    unsigned char* p_out = reinterpret_cast<unsigned char*>(out.begin());

    for (int t = 0; t < num_tiles; ++t) {
      int dst_r = tiles[t].r;
      int dst_c = tiles[t].c;
      int dst_x0 = col_starts[dst_c];
      int dst_x1 = col_ends[dst_c];
      int dst_y0 = row_starts[dst_r];
      int dst_y1 = row_ends[dst_r];
      int dst_tw = dst_x1 - dst_x0;
      int dst_th = dst_y1 - dst_y0;
      int tile_pixels = dst_tw * dst_th;

      int src_r = permuted[t].r;
      int src_c = permuted[t].c;
      int src_x0 = col_starts[src_c];
      int src_x1 = col_ends[src_c];
      int src_y0 = row_starts[src_r];
      int src_y1 = row_ends[src_r];
      int src_tw = src_x1 - src_x0;
      int src_th = src_y1 - src_y0;

      int rot = tile_rot[t];
      bool fx = tile_flip_x[t];
      bool fy = tile_flip_y[t];

      // Intra-tile pixel permutation (shuffles all pixels inside this block!)
      std::vector<int> pix_perm;
      if (shuffle_pixels && tile_pixels > 1) {
        pix_perm.resize(tile_pixels);
        for (int p = 0; p < tile_pixels; ++p) pix_perm[p] = p;
        for (int p = tile_pixels - 1; p > 0; --p) {
          int q = fast_rand() % (p + 1);
          std::swap(pix_perm[p], pix_perm[q]);
        }
      }

      for (int ch = 0; ch < channels; ++ch) {
        std::size_t ch_offset = static_cast<std::size_t>(ch) * plane_size;

        for (int dy = 0; dy < dst_th; ++dy) {
          int out_y = dst_y0 + dy;
          std::size_t out_row_idx = ch_offset + static_cast<std::size_t>(out_y) * w;

          for (int dx = 0; dx < dst_tw; ++dx) {
            int out_x = dst_x0 + dx;

            int eff_dx = dx;
            int eff_dy = dy;
            if (shuffle_pixels && tile_pixels > 1) {
              int p_idx = pix_perm[dy * dst_tw + dx];
              eff_dx = p_idx % dst_tw;
              eff_dy = p_idx / dst_tw;
            }

            double u = (dst_tw > 1) ? static_cast<double>(eff_dx) / (dst_tw - 1) : 0.0;
            double v = (dst_th > 1) ? static_cast<double>(eff_dy) / (dst_th - 1) : 0.0;

            if (fx) u = 1.0 - u;
            if (fy) v = 1.0 - v;

            double su = u;
            double sv = v;
            if (rot == 1) { // 90 deg
              su = 1.0 - v;
              sv = u;
            } else if (rot == 2) { // 180 deg
              su = 1.0 - u;
              sv = 1.0 - v;
            } else if (rot == 3) { // 270 deg
              su = v;
              sv = 1.0 - u;
            }

            int sx = src_x0 + std::min(src_tw - 1, static_cast<int>(su * (src_tw - 1) + 0.5));
            int sy = src_y0 + std::min(src_th - 1, static_cast<int>(sv * (src_th - 1) + 0.5));

            p_out[out_row_idx + out_x] = p_in[ch_offset + static_cast<std::size_t>(sy) * w + sx];
          }
        }
      }
    }

    out.attr("dim") = raw_vec.attr("dim");
    if (raw_vec.hasAttribute("class")) out.attr("class") = raw_vec.attr("class");
    if (raw_vec.hasAttribute("colormode")) out.attr("colormode") = raw_vec.attr("colormode");
    return out;

  } else {
    NumericVector out(total_elements);
    const double* p_in = num_vec.begin();
    double* p_out = out.begin();

    for (int t = 0; t < num_tiles; ++t) {
      int dst_r = tiles[t].r;
      int dst_c = tiles[t].c;
      int dst_x0 = col_starts[dst_c];
      int dst_x1 = col_ends[dst_c];
      int dst_y0 = row_starts[dst_r];
      int dst_y1 = row_ends[dst_r];
      int dst_tw = dst_x1 - dst_x0;
      int dst_th = dst_y1 - dst_y0;
      int tile_pixels = dst_tw * dst_th;

      int src_r = permuted[t].r;
      int src_c = permuted[t].c;
      int src_x0 = col_starts[src_c];
      int src_x1 = col_ends[src_c];
      int src_y0 = row_starts[src_r];
      int src_y1 = row_ends[src_r];
      int src_tw = src_x1 - src_x0;
      int src_th = src_y1 - src_y0;

      int rot = tile_rot[t];
      bool fx = tile_flip_x[t];
      bool fy = tile_flip_y[t];

      std::vector<int> pix_perm;
      if (shuffle_pixels && tile_pixels > 1) {
        pix_perm.resize(tile_pixels);
        for (int p = 0; p < tile_pixels; ++p) pix_perm[p] = p;
        for (int p = tile_pixels - 1; p > 0; --p) {
          int q = fast_rand() % (p + 1);
          std::swap(pix_perm[p], pix_perm[q]);
        }
      }

      for (int ch = 0; ch < channels; ++ch) {
        std::size_t ch_offset = static_cast<std::size_t>(ch) * plane_size;

        for (int dy = 0; dy < dst_th; ++dy) {
          int out_y = dst_y0 + dy;
          std::size_t out_row_idx = ch_offset + static_cast<std::size_t>(out_y) * w;

          for (int dx = 0; dx < dst_tw; ++dx) {
            int out_x = dst_x0 + dx;

            int eff_dx = dx;
            int eff_dy = dy;
            if (shuffle_pixels && tile_pixels > 1) {
              int p_idx = pix_perm[dy * dst_tw + dx];
              eff_dx = p_idx % dst_tw;
              eff_dy = p_idx / dst_tw;
            }

            double u = (dst_tw > 1) ? static_cast<double>(eff_dx) / (dst_tw - 1) : 0.0;
            double v = (dst_th > 1) ? static_cast<double>(eff_dy) / (dst_th - 1) : 0.0;

            if (fx) u = 1.0 - u;
            if (fy) v = 1.0 - v;

            double su = u;
            double sv = v;
            if (rot == 1) {
              su = 1.0 - v;
              sv = u;
            } else if (rot == 2) {
              su = 1.0 - u;
              sv = 1.0 - v;
            } else if (rot == 3) {
              su = v;
              sv = 1.0 - u;
            }

            int sx = src_x0 + std::min(src_tw - 1, static_cast<int>(su * (src_tw - 1) + 0.5));
            int sy = src_y0 + std::min(src_th - 1, static_cast<int>(sv * (src_th - 1) + 0.5));

            p_out[out_row_idx + out_x] = p_in[ch_offset + static_cast<std::size_t>(sy) * w + sx];
          }
        }
      }
    }

    out.attr("dim") = num_vec.attr("dim");
    if (num_vec.hasAttribute("class")) out.attr("class") = num_vec.attr("class");
    if (num_vec.hasAttribute("colormode")) out.attr("colormode") = num_vec.attr("colormode");
    return out;
  }
}

// ==============================================================================
// HIGH-PERFORMANCE PARALLEL YOLO LABEL AUDIT & OBJECT COUNTING (OpenMP)
// ==============================================================================

struct AuditIssueCpp {
  std::string file;
  std::string type;
  std::string message;
};

struct AuditThreadStateCpp {
  std::vector<AuditIssueCpp> issues;
  int repaired_count = 0;
  int n_objects = 0;
};

static inline bool audit_parse_double_tokens(const char* p, const char* end, std::vector<double>& out_vals) {
  out_vals.clear();
  while (p < end) {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) p++;
    if (p >= end) break;

    const char* token_start = p;
    bool has_sign = (*p == '+' || *p == '-');
    if (has_sign) p++;
    bool has_digits = false;
    while (p < end && *p >= '0' && *p <= '9') { has_digits = true; p++; }
    if (p < end && *p == '.') {
      p++;
      while (p < end && *p >= '0' && *p <= '9') { has_digits = true; p++; }
    }
    if (has_digits && p < end && (*p == 'e' || *p == 'E')) {
      p++;
      if (p < end && (*p == '+' || *p == '-')) p++;
      while (p < end && *p >= '0' && *p <= '9') p++;
    }
    if (!has_digits || (p < end && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n')) {
      return false; // Malformed token
    }
    const char* next_p = nullptr;
    double val = fast_atof(token_start, &next_p);
    out_vals.push_back(val);
  }
  return true;
}

// [[Rcpp::export]]
int fast_count_yolo_objects_cpp(Rcpp::CharacterVector lbl_files, int num_threads = 0) {
  int n_files = lbl_files.size();
  if (n_files == 0) return 0;

  std::vector<std::string> v_lbl(n_files);
  for (int i = 0; i < n_files; ++i) {
    v_lbl[i] = Rcpp::as<std::string>(lbl_files[i]);
  }

  int n_threads_use = 1;
#ifdef _OPENMP
  if (num_threads > 0) n_threads_use = num_threads;
  else n_threads_use = std::max(1, omp_get_max_threads());
#endif

  int total_objs = 0;
#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads_use) reduction(+:total_objs) if(n_files > 16)
#endif
  {
#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 64)
#endif
    for (int i = 0; i < n_files; ++i) {
      const std::string& lpath = v_lbl[i];
      if (lpath.empty()) continue;
      FILE* fp = std::fopen(lpath.c_str(), "rb");
      if (!fp) continue;
      std::fseek(fp, 0, SEEK_END);
      long fsize = std::ftell(fp);
      std::fseek(fp, 0, SEEK_SET);
      if (fsize <= 0) {
        std::fclose(fp);
        continue;
      }
      std::string buf(fsize, '\0');
      size_t read_bytes = std::fread(&buf[0], 1, fsize, fp);
      std::fclose(fp);

      const char* p = buf.c_str();
      const char* end = p + read_bytes;
      int objs_in_file = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) p++;
        if (p < end && *p != '\r' && *p != '\n') {
          objs_in_file++;
          while (p < end && *p != '\r' && *p != '\n') p++;
        }
        while (p < end && (*p == '\r' || *p == '\n')) p++;
      }
      total_objs += objs_in_file;
    }
  }
  return total_objs;
}

// [[Rcpp::export]]
Rcpp::List fast_audit_yolo_labels_cpp(
    Rcpp::CharacterVector lbl_files,
    int n_classes = 0,
    bool repair = false,
    int num_threads = 0
) {
  int n_files = lbl_files.size();
  std::vector<std::string> v_lbl(n_files);
  for (int i = 0; i < n_files; ++i) {
    v_lbl[i] = Rcpp::as<std::string>(lbl_files[i]);
  }

  int n_threads_use = 1;
#ifdef _OPENMP
  if (num_threads > 0) n_threads_use = num_threads;
  else n_threads_use = std::max(1, omp_get_max_threads());
#endif

  std::vector<AuditThreadStateCpp> t_states(n_threads_use);

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads_use) if(n_files > 16)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    AuditThreadStateCpp& ts = t_states[tid];
    std::vector<double> line_vals;
    line_vals.reserve(32);
    std::vector<std::string> fixed_lines;
    fixed_lines.reserve(16);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 32)
#endif
    for (int i = 0; i < n_files; ++i) {
      const std::string& lpath = v_lbl[i];
      if (lpath.empty()) continue;

      FILE* fp = std::fopen(lpath.c_str(), "rb");
      if (!fp) continue;

      std::fseek(fp, 0, SEEK_END);
      long fsize = std::ftell(fp);
      std::fseek(fp, 0, SEEK_SET);

      if (fsize <= 0) {
        std::fclose(fp);
        continue;
      }

      std::string buf(fsize, '\0');
      size_t read_bytes = std::fread(&buf[0], 1, fsize, fp);
      std::fclose(fp);
      buf.resize(read_bytes);

      if (buf.empty()) continue;

      fixed_lines.clear();
      bool file_modified = false;

      const char* p = buf.c_str();
      const char* end = p + buf.size();

      while (p < end) {
        const char* l_start = p;
        while (p < end && *p != '\r' && *p != '\n') p++;
        const char* l_end = p;
        while (p < end && (*p == '\r' || *p == '\n')) p++;

        // Trim leading and trailing whitespace
        while (l_start < l_end && (*l_start == ' ' || *l_start == '\t')) l_start++;
        while (l_end > l_start && (*(l_end - 1) == ' ' || *(l_end - 1) == '\t')) l_end--;

        if (l_start >= l_end) continue; // blank line

        std::string line_str(l_start, l_end);

        bool parse_ok = audit_parse_double_tokens(l_start, l_end, line_vals);
        if (!parse_ok || line_vals.size() < 5) {
          ts.issues.push_back({lpath, "malformed_annotation", "Corrupt annotation line: '" + line_str + "'"});
          file_modified = true;
          continue;
        }

        int cid = static_cast<int>(line_vals[0]);
        if (cid < 0 || (n_classes > 0 && cid >= n_classes)) {
          ts.issues.push_back({
            lpath, "invalid_class_id",
            "Class index " + std::to_string(cid) + " is out of range (0 to " + std::to_string(n_classes - 1) + ")"
          });
        }

        if (line_vals.size() == 5) {
          // Bounding Box: class xc yc w h
          double xc = line_vals[1];
          double yc = line_vals[2];
          double w = line_vals[3];
          double h = line_vals[4];
          bool has_bad_geom = (w <= 0.0 || h <= 0.0 || xc < 0.0 || xc > 1.0 || yc < 0.0 || yc > 1.0);

          if (has_bad_geom) {
            ts.issues.push_back({lpath, "geometry_out_of_bounds", "Invalid bounding box: " + line_str});
            if (repair) {
              if (w > 0.0 && h > 0.0) {
                xc = std::max(0.001, std::min(0.999, xc));
                yc = std::max(0.001, std::min(0.999, yc));
                w = std::min(w, 2.0 * std::min(xc, 1.0 - xc));
                h = std::min(h, 2.0 * std::min(yc, 1.0 - yc));
                char fbuf[128];
                std::snprintf(fbuf, sizeof(fbuf), "%d %.6f %.6f %.6f %.6f", cid, xc, yc, w, h);
                fixed_lines.push_back(std::string(fbuf));
                file_modified = true;
                ts.repaired_count++;
                ts.n_objects++;
              } else {
                // Discard zero-area box
                file_modified = true;
                ts.repaired_count++;
              }
              continue;
            }
          }
          fixed_lines.push_back(line_str);
          ts.n_objects++;
        } else {
          // Segmentation Polygon: class x1 y1 x2 y2 ...
          size_t n_coords = line_vals.size() - 1;
          if (n_coords % 2 != 0 || n_coords < 6) {
            ts.issues.push_back({lpath, "malformed_polygon", "Polygon has odd number or < 3 vertices: " + line_str});
            file_modified = true;
            continue;
          }

          bool oob = false;
          for (size_t k = 1; k < line_vals.size(); ++k) {
            if (line_vals[k] < 0.0 || line_vals[k] > 1.0) {
              oob = true;
              break;
            }
          }

          if (oob) {
            ts.issues.push_back({lpath, "polygon_out_of_bounds", "Polygon vertices leak past [0, 1] range: " + line_str});
            if (repair) {
              std::string rep_str = std::to_string(cid);
              char cbuf[32];
              for (size_t k = 1; k < line_vals.size(); ++k) {
                double c = std::max(0.0, std::min(1.0, line_vals[k]));
                std::snprintf(cbuf, sizeof(cbuf), " %.6f", c);
                rep_str += cbuf;
              }
              fixed_lines.push_back(rep_str);
              file_modified = true;
              ts.repaired_count++;
              ts.n_objects++;
              continue;
            }
          }
          fixed_lines.push_back(line_str);
          ts.n_objects++;
        }
      }

      if (repair && file_modified) {
        FILE* out_fp = std::fopen(lpath.c_str(), "wb");
        if (out_fp) {
          for (size_t l = 0; l < fixed_lines.size(); ++l) {
            std::fputs(fixed_lines[l].c_str(), out_fp);
            std::fputc('\n', out_fp);
          }
          std::fclose(out_fp);
        }
      }
    }
  }

  int total_issues = 0;
  int total_repaired = 0;
  int total_objects = 0;
  for (int t = 0; t < n_threads_use; ++t) {
    total_issues += static_cast<int>(t_states[t].issues.size());
    total_repaired += t_states[t].repaired_count;
    total_objects += t_states[t].n_objects;
  }

  Rcpp::CharacterVector res_files(total_issues);
  Rcpp::CharacterVector res_types(total_issues);
  Rcpp::CharacterVector res_messages(total_issues);

  int cur = 0;
  for (int t = 0; t < n_threads_use; ++t) {
    for (const auto& iss : t_states[t].issues) {
      res_files[cur] = iss.file;
      res_types[cur] = iss.type;
      res_messages[cur] = iss.message;
      cur++;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("files") = res_files,
    Rcpp::Named("types") = res_types,
    Rcpp::Named("messages") = res_messages,
    Rcpp::Named("repaired_count") = total_repaired,
    Rcpp::Named("n_objects") = total_objects
  );
}

// ==============================================================================
// HIGH-SPEED VECTORIZED COCO EVALUATION (mAP50, mAP75, mAP50-95, PR-CURVE, CM)
// ==============================================================================

inline double calc_fast_box_iou(double ax1, double ay1, double ax2, double ay2,
                                double bx1, double by1, double bx2, double by2) {
  double x1 = (ax1 > bx1) ? ax1 : bx1;
  double y1 = (ay1 > by1) ? ay1 : by1;
  double x2 = (ax2 < bx2) ? ax2 : bx2;
  double y2 = (ay2 < by2) ? ay2 : by2;
  double w = (x2 > x1) ? (x2 - x1) : 0.0;
  double h = (y2 > y1) ? (y2 - y1) : 0.0;
  double inter = w * h;
  if (inter <= 0.0) return 0.0;
  double areaA = (ax2 - ax1) * (ay2 - ay1);
  double areaB = (bx2 - bx1) * (by2 - by1);
  double uni = areaA + areaB - inter;
  return (uni > 0.0) ? (inter / uni) : 0.0;
}

struct GtBoxInfo {
  int cls;
  double x1, y1, x2, y2;
};

// [[Rcpp::export]]
Rcpp::List fast_coco_eval_cpp(
    Rcpp::IntegerVector pred_img_id,
    Rcpp::IntegerVector pred_cls_id,
    Rcpp::NumericVector pred_score,
    Rcpp::NumericVector pred_x1,
    Rcpp::NumericVector pred_y1,
    Rcpp::NumericVector pred_x2,
    Rcpp::NumericVector pred_y2,
    Rcpp::IntegerVector gt_img_id,
    Rcpp::IntegerVector gt_cls_id,
    Rcpp::NumericVector gt_x1,
    Rcpp::NumericVector gt_y1,
    Rcpp::NumericVector gt_x2,
    Rcpp::NumericVector gt_y2,
    int num_images,
    int num_classes,
    Rcpp::NumericVector iou_thresholds,
    double eval_iou_th = 0.50
) {
  int n_preds = pred_score.size();
  int n_gts = gt_img_id.size();
  int n_th = iou_thresholds.size();

  // Group GT by image
  std::vector<std::vector<GtBoxInfo>> img_gts(num_images);
  for (int g = 0; g < n_gts; ++g) {
    int img = gt_img_id[g];
    if (img >= 0 && img < num_images) {
      img_gts[img].push_back({gt_cls_id[g], gt_x1[g], gt_y1[g], gt_x2[g], gt_y2[g]});
    }
  }

  // Sort prediction indices descending by score
  std::vector<int> sort_order(n_preds);
  std::iota(sort_order.begin(), sort_order.end(), 0);
  std::sort(sort_order.begin(), sort_order.end(), [&](int a, int b) {
    return pred_score[a] > pred_score[b];
  });

  Rcpp::NumericVector ap_values(n_th);
  std::vector<int> eval_cum_tp(n_preds, 0);
  std::vector<int> eval_cum_fp(n_preds, 0);

  // Multi-threshold mAP computation
  for (int t = 0; t < n_th; ++t) {
    double th = iou_thresholds[t];
    std::vector<std::vector<char>> gt_matched(num_images);
    for (int i = 0; i < num_images; ++i) {
      gt_matched[i].resize(img_gts[i].size(), 0);
    }

    std::vector<char> tp(n_preds, 0);
    std::vector<char> fp(n_preds, 0);

    for (int p = 0; p < n_preds; ++p) {
      int idx = sort_order[p];
      int img = pred_img_id[idx];
      int pcls = pred_cls_id[idx];
      double px1 = pred_x1[idx], py1 = pred_y1[idx], px2 = pred_x2[idx], py2 = pred_y2[idx];

      if (img < 0 || img >= num_images) {
        fp[p] = 1;
        continue;
      }

      const auto& gts = img_gts[img];
      double best_iou = 0.0;
      int best_g = -1;

      for (size_t g = 0; g < gts.size(); ++g) {
        if (gts[g].cls == pcls && !gt_matched[img][g]) {
          double iou = calc_fast_box_iou(px1, py1, px2, py2, gts[g].x1, gts[g].y1, gts[g].x2, gts[g].y2);
          if (iou > best_iou) {
            best_iou = iou;
            best_g = static_cast<int>(g);
          }
        }
      }

      if (best_iou >= th && best_g >= 0) {
        tp[p] = 1;
        gt_matched[img][best_g] = 1;
      } else {
        fp[p] = 1;
      }
    }

    // Cumulative sums
    std::vector<int> cum_tp(n_preds, 0);
    std::vector<int> cum_fp(n_preds, 0);
    int c_tp = 0, c_fp = 0;
    for (int p = 0; p < n_preds; ++p) {
      c_tp += tp[p];
      c_fp += fp[p];
      cum_tp[p] = c_tp;
      cum_fp[p] = c_fp;
    }

    if (std::abs(th - eval_iou_th) < 1e-4) {
      eval_cum_tp = cum_tp;
      eval_cum_fp = cum_fp;
    }

    // 101-point COCO interpolated AP
    double ap = 0.0;
    if (n_gts > 0 && n_preds > 0) {
      int p_ptr = n_preds - 1;
      double max_p = 0.0;
      for (int r_step = 100; r_step >= 0; --r_step) {
        double r_th = r_step / 100.0;
        while (p_ptr >= 0 && (cum_tp[p_ptr] / static_cast<double>(n_gts)) >= r_th) {
          double p_val = cum_tp[p_ptr] / static_cast<double>(cum_tp[p_ptr] + cum_fp[p_ptr]);
          if (p_val > max_p) max_p = p_val;
          --p_ptr;
        }
        ap += max_p;
      }
      ap /= 101.0;
    }
    ap_values[t] = ap;
  }

  // Bipartite matching at eval_iou_th for Confusion Matrix and per-image counts
  int bg_idx = num_classes;
  int cm_dim = num_classes + 1;
  std::vector<int> cm(cm_dim * cm_dim, 0);

  std::vector<int> tp_per_img(num_images, 0);
  std::vector<int> fp_per_img(num_images, 0);
  std::vector<int> fn_per_img(num_images, 0);
  std::vector<double> sum_iou_per_img(num_images, 0.0);
  std::vector<int> count_tp_per_img(num_images, 0);
  std::vector<double> all_tp_ious;
  all_tp_ious.reserve(std::min(n_preds, n_gts));

  std::vector<std::vector<char>> eval_gt_matched(num_images);
  for (int i = 0; i < num_images; ++i) {
    eval_gt_matched[i].resize(img_gts[i].size(), 0);
  }

  for (int p = 0; p < n_preds; ++p) {
    int idx = sort_order[p];
    int img = pred_img_id[idx];
    int pcls = pred_cls_id[idx];
    double px1 = pred_x1[idx], py1 = pred_y1[idx], px2 = pred_x2[idx], py2 = pred_y2[idx];

    if (img < 0 || img >= num_images) continue;

    const auto& gts = img_gts[img];
    double best_iou = 0.0;
    int best_g = -1;

    for (size_t g = 0; g < gts.size(); ++g) {
      if (gts[g].cls == pcls && !eval_gt_matched[img][g]) {
        double iou = calc_fast_box_iou(px1, py1, px2, py2, gts[g].x1, gts[g].y1, gts[g].x2, gts[g].y2);
        if (iou > best_iou) {
          best_iou = iou;
          best_g = static_cast<int>(g);
        }
      }
    }

    if (best_iou >= eval_iou_th && best_g >= 0) {
      eval_gt_matched[img][best_g] = 1;
      tp_per_img[img]++;
      count_tp_per_img[img]++;
      sum_iou_per_img[img] += best_iou;
      all_tp_ious.push_back(best_iou);
      if (pcls >= 0 && pcls < num_classes) {
        cm[pcls * cm_dim + pcls]++;
      }
    } else {
      // Check if overlapped another class (misclassification)
      double best_other_iou = 0.0;
      int best_other_g = -1;
      for (size_t g = 0; g < gts.size(); ++g) {
        if (!eval_gt_matched[img][g]) {
          double iou = calc_fast_box_iou(px1, py1, px2, py2, gts[g].x1, gts[g].y1, gts[g].x2, gts[g].y2);
          if (iou > best_other_iou) {
            best_other_iou = iou;
            best_other_g = static_cast<int>(g);
          }
        }
      }
      if (best_other_iou >= eval_iou_th && best_other_g >= 0) {
        eval_gt_matched[img][best_other_g] = 1;
        fp_per_img[img]++;
        int act_cls = gts[best_other_g].cls;
        if (act_cls >= 0 && act_cls < num_classes && pcls >= 0 && pcls < num_classes) {
          cm[act_cls * cm_dim + pcls]++;
        }
      } else {
        // False positive on background
        fp_per_img[img]++;
        if (pcls >= 0 && pcls < num_classes) {
          cm[bg_idx * cm_dim + pcls]++;
        }
      }
    }
  }

  // Unmatched GT are False Negatives
  for (int i = 0; i < num_images; ++i) {
    const auto& gts = img_gts[i];
    int n_unmatched = 0;
    for (size_t g = 0; g < gts.size(); ++g) {
      if (!eval_gt_matched[i][g]) {
        n_unmatched++;
        int act_cls = gts[g].cls;
        if (act_cls >= 0 && act_cls < num_classes) {
          cm[act_cls * cm_dim + bg_idx]++;
        }
      }
    }
    fn_per_img[i] = n_unmatched;
    if (gts.empty() && count_tp_per_img[i] == 0 && fp_per_img[i] == 0) {
      cm[bg_idx * cm_dim + bg_idx]++;
    }
  }

  std::vector<double> mean_iou_per_img(num_images, 0.0);
  for (int i = 0; i < num_images; ++i) {
    if (count_tp_per_img[i] > 0) {
      mean_iou_per_img[i] = sum_iou_per_img[i] / count_tp_per_img[i];
    }
  }

  Rcpp::IntegerMatrix cm_out(cm_dim, cm_dim);
  for (int r = 0; r < cm_dim; ++r) {
    for (int c = 0; c < cm_dim; ++c) {
      cm_out(r, c) = cm[r * cm_dim + c];
    }
  }

  // Subsample PR curve points (up to 500 points)
  int pr_step = std::max(1, static_cast<int>(n_preds / 500));
  std::vector<double> pr_r, pr_p, pr_c;
  pr_r.push_back(0.0); pr_p.push_back(1.0); pr_c.push_back(1.0);
  for (int p = 0; p < n_preds; p += pr_step) {
    pr_r.push_back(n_gts > 0 ? eval_cum_tp[p] / static_cast<double>(n_gts) : 0.0);
    pr_p.push_back(eval_cum_tp[p] / static_cast<double>(std::max(1, eval_cum_tp[p] + eval_cum_fp[p])));
    pr_c.push_back(pred_score[sort_order[p]]);
  }
  pr_r.push_back(1.0); pr_p.push_back(0.0); pr_c.push_back(0.0);

  // F1 curve points across 19 confidence thresholds (0.05 to 0.95)
  std::vector<double> f1_confs(19);
  std::vector<double> f1_scores(19, 0.0);
  for (int ci = 0; ci < 19; ++ci) {
    double cf = 0.05 + ci * 0.05;
    f1_confs[ci] = cf;
    int last_p = -1;
    for (int p = n_preds - 1; p >= 0; --p) {
      if (pred_score[sort_order[p]] >= cf) {
        last_p = p;
        break;
      }
    }
    if (last_p >= 0 && (eval_cum_tp[last_p] + eval_cum_fp[last_p]) > 0 && n_gts > 0) {
      double r = eval_cum_tp[last_p] / static_cast<double>(n_gts);
      double prec = eval_cum_tp[last_p] / static_cast<double>(eval_cum_tp[last_p] + eval_cum_fp[last_p]);
      if (r + prec > 0.0) {
        f1_scores[ci] = 2.0 * (r * prec) / (r + prec);
      }
    }
  }

  double map50 = ap_values.size() > 0 ? ap_values[0] : 0.0;
  double map75 = ap_values.size() > 5 ? ap_values[5] : map50;
  double mean_map = (ap_values.size() > 0) ? Rcpp::mean(ap_values) : 0.0;

  return Rcpp::List::create(
    Rcpp::Named("map50") = map50,
    Rcpp::Named("map75") = map75,
    Rcpp::Named("map50_95") = mean_map,
    Rcpp::Named("map_all") = ap_values,
    Rcpp::Named("tp_per_img") = Rcpp::wrap(tp_per_img),
    Rcpp::Named("fp_per_img") = Rcpp::wrap(fp_per_img),
    Rcpp::Named("fn_per_img") = Rcpp::wrap(fn_per_img),
    Rcpp::Named("mean_iou_per_img") = Rcpp::wrap(mean_iou_per_img),
    Rcpp::Named("confusion_matrix") = cm_out,
    Rcpp::Named("matched_ious") = Rcpp::wrap(all_tp_ious),
    Rcpp::Named("pr_recall") = Rcpp::wrap(pr_r),
    Rcpp::Named("pr_precision") = Rcpp::wrap(pr_p),
    Rcpp::Named("pr_confidence") = Rcpp::wrap(pr_c),
    Rcpp::Named("f1_confidence") = Rcpp::wrap(f1_confs),
    Rcpp::Named("f1_values") = Rcpp::wrap(f1_scores)
  );
}

// -----------------------------------------------------------------------------
// 10. FAST PERCEPTUAL HASH (dHash 64-bit) & CROSS-SPLIT DATA LEAKAGE DETECTION
// -----------------------------------------------------------------------------

static inline uint64_t compute_single_dhash(const std::string& path) {
  int w = 0, h = 0, channels = 0;
  unsigned char* img = stbi_load(path.c_str(), &w, &h, &channels, 1);
  if (!img || w < 9 || h < 8) {
    if (img) stbi_image_free(img);
    return 0ULL;
  }

  // Downsample to 9 x 8 via box averaging
  unsigned char small_grid[8][9];
  for (int y = 0; y < 8; y++) {
    int y_start = (y * h) / 8;
    int y_end   = ((y + 1) * h) / 8;
    if (y_end <= y_start) y_end = y_start + 1;

    for (int x = 0; x < 9; x++) {
      int x_start = (x * w) / 9;
      int x_end   = ((x + 1) * w) / 9;
      if (x_end <= x_start) x_end = x_start + 1;

      uint64_t sum = 0;
      int count = 0;
      for (int sy = y_start; sy < y_end; sy++) {
        const unsigned char* row_ptr = img + (size_t)sy * w;
        for (int sx = x_start; sx < x_end; sx++) {
          sum += row_ptr[sx];
          count++;
        }
      }
      small_grid[y][x] = (unsigned char)(count > 0 ? (sum / count) : 0);
    }
  }

  stbi_image_free(img);

  // Compute 64-bit horizontal difference hash
  uint64_t hash = 0ULL;
  int bit_idx = 0;
  for (int y = 0; y < 8; y++) {
    for (int x = 0; x < 8; x++) {
      if (small_grid[y][x] > small_grid[y][x + 1]) {
        hash |= (1ULL << bit_idx);
      }
      bit_idx++;
    }
  }
  return hash;
}

// [[Rcpp::export]]
List fast_dhash_leakage_cpp(CharacterVector img_files,
                            CharacterVector splits,
                            int max_distance = 4,
                            int num_threads = 0) {
  int n = img_files.size();
  std::vector<std::string> paths(n);
  std::vector<std::string> s_splits(n);
  for (int i = 0; i < n; i++) {
    paths[i] = Rcpp::as<std::string>(img_files[i]);
    s_splits[i] = (splits.size() == n) ? Rcpp::as<std::string>(splits[i]) : "train";
  }

  std::vector<uint64_t> hashes(n, 0ULL);

  int nthreads = 1;
#ifdef _OPENMP
  if (num_threads <= 0) {
    nthreads = omp_get_max_threads();
  } else {
    nthreads = num_threads;
  }
#endif

#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (int i = 0; i < n; i++) {
    hashes[i] = compute_single_dhash(paths[i]);
  }

  CharacterVector hex_hashes(n);
  char hex_buf[24];
  for (int i = 0; i < n; i++) {
    std::snprintf(hex_buf, sizeof(hex_buf), "%016llx", (unsigned long long)hashes[i]);
    hex_hashes[i] = hex_buf;
  }

  std::vector<int> train_idx;
  std::vector<int> eval_idx;
  for (int i = 0; i < n; i++) {
    if (hashes[i] == 0ULL) continue;
    if (s_splits[i] == "train") {
      train_idx.push_back(i);
    } else if (s_splits[i] == "val" || s_splits[i] == "test") {
      eval_idx.push_back(i);
    }
  }

  std::vector<std::string> leak_tr;
  std::vector<std::string> leak_ev;
  std::vector<std::string> leak_sp;
  std::vector<int> leak_dist;
  std::vector<double> leak_sim;

  for (size_t i = 0; i < train_idx.size(); i++) {
    int t_i = train_idx[i];
    uint64_t h_tr = hashes[t_i];
    for (size_t j = 0; j < eval_idx.size(); j++) {
      int e_j = eval_idx[j];
      uint64_t h_ev = hashes[e_j];
      int dist = (int)std::bitset<64>(h_tr ^ h_ev).count();
      if (dist <= max_distance) {
        leak_tr.push_back(paths[t_i]);
        leak_ev.push_back(paths[e_j]);
        leak_sp.push_back(s_splits[e_j]);
        leak_dist.push_back(dist);
        leak_sim.push_back(std::round((1.0 - (double)dist / 64.0) * 1000.0) / 10.0);
      }
    }
  }

  std::vector<std::string> dup_f1;
  std::vector<std::string> dup_f2;
  std::vector<std::string> dup_sp1;
  std::vector<std::string> dup_sp2;
  std::vector<int> dup_dist;
  std::vector<double> dup_sim;

  int max_check = std::min(n, 15000);
  for (int i = 0; i < max_check; i++) {
    if (hashes[i] == 0ULL) continue;
    for (int j = i + 1; j < max_check; j++) {
      if (hashes[j] == 0ULL) continue;
      bool is_leak = (s_splits[i] == "train" && (s_splits[j] == "val" || s_splits[j] == "test")) ||
                     (s_splits[j] == "train" && (s_splits[i] == "val" || s_splits[i] == "test"));
      if (is_leak) continue;

      int dist = (int)std::bitset<64>(hashes[i] ^ hashes[j]).count();
      if (dist <= max_distance) {
        dup_f1.push_back(paths[i]);
        dup_f2.push_back(paths[j]);
        dup_sp1.push_back(s_splits[i]);
        dup_sp2.push_back(s_splits[j]);
        dup_dist.push_back(dist);
        dup_sim.push_back(std::round((1.0 - (double)dist / 64.0) * 1000.0) / 10.0);
      }
    }
  }

  DataFrame leaks_df = DataFrame::create(
    Named("img_train") = wrap(leak_tr),
    Named("img_eval") = wrap(leak_ev),
    Named("split_eval") = wrap(leak_sp),
    Named("hamming_dist") = wrap(leak_dist),
    Named("similarity_pct") = wrap(leak_sim),
    Named("stringsAsFactors") = false
  );

  DataFrame dups_df = DataFrame::create(
    Named("img_1") = wrap(dup_f1),
    Named("img_2") = wrap(dup_f2),
    Named("split_1") = wrap(dup_sp1),
    Named("split_2") = wrap(dup_sp2),
    Named("hamming_dist") = wrap(dup_dist),
    Named("similarity_pct") = wrap(dup_sim),
    Named("stringsAsFactors") = false
  );

  return List::create(
    Named("hashes") = hex_hashes,
    Named("leaks") = leaks_df,
    Named("duplicates") = dups_df,
    Named("leak_count") = (int)leak_tr.size(),
    Named("duplicate_count") = (int)dup_f1.size()
  );
}







