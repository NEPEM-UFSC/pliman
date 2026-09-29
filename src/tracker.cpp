#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

using namespace Rcpp;

// 2D Cross-product of vector AB and AC
static inline double cross_product_2d(double ax, double ay, double bx, double by, double cx, double cy) {
  return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

// Check if segment P1->P2 intersects segment L1->L2
// Returns 1 for positive crossing, -1 for negative crossing, 0 for no crossing
static int check_segment_intersection(
    double px1, double py1, double px2, double py2,
    double lx1, double ly1, double lx2, double ly2
) {
  // 1. Fast exact path for vertical counting line
  if (std::abs(lx2 - lx1) < 1.0) {
    double min_y = std::min(ly1, ly2) - 16.0;
    double max_y = std::max(ly1, ly2) + 16.0;
    double dx = px2 - px1;
    if (std::abs(dx) > 1e-5) {
      bool crossed = (px1 <= lx1 && px2 > lx1) || (px1 >= lx1 && px2 < lx1);
      if (crossed) {
        double t = (lx1 - px1) / dx;
        double y_cross = py1 + t * (py2 - py1);
        if (y_cross >= min_y && y_cross <= max_y) {
          return (px2 > px1) ? 1 : -1;
        }
      }
    }
  }

  // 2. Fast exact path for horizontal counting line
  if (std::abs(ly2 - ly1) < 1.0) {
    double min_x = std::min(lx1, lx2) - 16.0;
    double max_x = std::max(lx1, lx2) + 16.0;
    double dy = py2 - py1;
    if (std::abs(dy) > 1e-5) {
      bool crossed = (py1 <= ly1 && py2 > ly1) || (py1 >= ly1 && py2 < ly1);
      if (crossed) {
        double t = (ly1 - py1) / dy;
        double x_cross = px1 + t * (px2 - px1);
        if (x_cross >= min_x && x_cross <= max_x) {
          return (py2 > py1) ? 1 : -1;
        }
      }
    }
  }

  // 3. General 2D line segment intersection
  double dx = px2 - px1;
  double dy = py2 - py1;
  double lx = lx2 - lx1;
  double ly = ly2 - ly1;
  double det = lx * dy - ly * dx;
  if (std::abs(det) > 1e-7) {
    double t = (lx * (py1 - ly1) - ly * (px1 - lx1)) / (-det);
    double u = (dx * (ly1 - py1) - dy * (lx1 - px1)) / det;
    if (t > 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
      return (det > 0.0) ? 1 : -1;
    }
  }

  return 0;
}

// Compute IoU between two bounding boxes
static inline double compute_iou(double ax1, double ay1, double ax2, double ay2,
                                 double bx1, double by1, double bx2, double by2) {
  double ix1 = std::max(ax1, bx1);
  double iy1 = std::max(ay1, by1);
  double ix2 = std::min(ax2, bx2);
  double iy2 = std::min(ay2, by2);

  double iw = std::max(0.0, ix2 - ix1);
  double ih = std::max(0.0, iy2 - iy1);
  double iarea = iw * ih;

  double area_a = std::max(0.0, ax2 - ax1) * std::max(0.0, ay2 - ay1);
  double area_b = std::max(0.0, bx2 - bx1) * std::max(0.0, by2 - by1);
  double uarea = area_a + area_b - iarea;

  return (uarea > 1e-6) ? (iarea / uarea) : 0.0;
}

// [[Rcpp::export]]
Rcpp::List update_tracker_cpp(
    Rcpp::NumericVector xmin,
    Rcpp::NumericVector ymin,
    Rcpp::NumericVector xmax,
    Rcpp::NumericVector ymax,
    Rcpp::CharacterVector labels,
    Rcpp::NumericVector scores,
    Rcpp::List state,
    Rcpp::Nullable<Rcpp::NumericVector> count_line = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> roi = R_NilValue,
    double max_dist = 120.0,
    double min_iou = 0.25,
    int max_lost = 15,
    int max_history = 30
) {
  // Extract state variables
  int next_id = 1;
  if (state.containsElementNamed("next_id")) {
    next_id = Rcpp::as<int>(state["next_id"]);
  }

  int total_count = 0;
  if (state.containsElementNamed("total_count")) {
    total_count = Rcpp::as<int>(state["total_count"]);
  }

  Rcpp::List active_tracks;
  if (state.containsElementNamed("tracks")) {
    active_tracks = Rcpp::as<Rcpp::List>(state["tracks"]);
  }

  Rcpp::List counts_by_class;
  if (state.containsElementNamed("counts_by_class")) {
    counts_by_class = Rcpp::as<Rcpp::List>(state["counts_by_class"]);
  }

  Rcpp::List crossing_events;
  if (state.containsElementNamed("crossing_events")) {
    crossing_events = Rcpp::as<Rcpp::List>(state["crossing_events"]);
  }

  int num_dets = xmin.size();
  int num_tracks = active_tracks.size();

  // Internal structure for C++ tracking
  struct TrackObj {
    int id;
    std::string label;
    double xmin, ymin, xmax, ymax;
    double cx, cy;
    double prev_cx, prev_cy;
    double vx, vy;
    std::vector<double> history_x;
    std::vector<double> history_y;
    int lost;
    int age;
    bool counted;
    int count_flash; // visual flash counter
    std::string last_dir;
  };

  std::vector<TrackObj> tracks(num_tracks);
  for (int i = 0; i < num_tracks; ++i) {
    Rcpp::List t = active_tracks[i];
    tracks[i].id = Rcpp::as<int>(t["id"]);
    tracks[i].label = Rcpp::as<std::string>(t["label"]);
    tracks[i].xmin = Rcpp::as<double>(t["xmin"]);
    tracks[i].ymin = Rcpp::as<double>(t["ymin"]);
    tracks[i].xmax = Rcpp::as<double>(t["xmax"]);
    tracks[i].ymax = Rcpp::as<double>(t["ymax"]);
    tracks[i].cx = Rcpp::as<double>(t["cx"]);
    tracks[i].cy = Rcpp::as<double>(t["cy"]);
    tracks[i].prev_cx = Rcpp::as<double>(t["prev_cx"]);
    tracks[i].prev_cy = Rcpp::as<double>(t["prev_cy"]);
    tracks[i].vx = t.containsElementNamed("vx") ? Rcpp::as<double>(t["vx"]) : 0.0;
    tracks[i].vy = t.containsElementNamed("vy") ? Rcpp::as<double>(t["vy"]) : 0.0;

    Rcpp::NumericVector hx = t["history_x"];
    Rcpp::NumericVector hy = t["history_y"];
    tracks[i].history_x.assign(hx.begin(), hx.end());
    tracks[i].history_y.assign(hy.begin(), hy.end());

    tracks[i].lost = Rcpp::as<int>(t["lost"]);
    tracks[i].age = Rcpp::as<int>(t["age"]);
    tracks[i].counted = Rcpp::as<bool>(t["counted"]);
    tracks[i].count_flash = t.containsElementNamed("count_flash") ? Rcpp::as<int>(t["count_flash"]) : 0;
    tracks[i].last_dir = t.containsElementNamed("last_dir") ? Rcpp::as<std::string>(t["last_dir"]) : "";
  }

  // Precompute detection centroids
  std::vector<double> ncx(num_dets), ncy(num_dets);
  for (int j = 0; j < num_dets; ++j) {
    ncx[j] = (xmin[j] + xmax[j]) / 2.0;
    ncy[j] = (ymin[j] + ymax[j]) / 2.0;
  }

  // Precompute track predicted positions using linear velocity
  std::vector<double> pred_cx(num_tracks), pred_cy(num_tracks);
  std::vector<double> pred_xmin(num_tracks), pred_ymin(num_tracks);
  std::vector<double> pred_xmax(num_tracks), pred_ymax(num_tracks);
  for (int i = 0; i < num_tracks; ++i) {
    double step_factor = (double)(tracks[i].lost + 1);
    pred_cx[i] = tracks[i].cx + tracks[i].vx * step_factor;
    pred_cy[i] = tracks[i].cy + tracks[i].vy * step_factor;
    double tw = tracks[i].xmax - tracks[i].xmin;
    double th = tracks[i].ymax - tracks[i].ymin;
    pred_xmin[i] = pred_cx[i] - tw / 2.0;
    pred_xmax[i] = pred_cx[i] + tw / 2.0;
    pred_ymin[i] = pred_cy[i] - th / 2.0;
    pred_ymax[i] = pred_cy[i] + th / 2.0;
  }

  std::vector<bool> det_matched(num_dets, false);
  std::vector<bool> track_matched(num_tracks, false);
  std::vector<int> det_to_track(num_dets, -1);

  struct MatchPair {
    int track_idx;
    int det_idx;
    double score;
  };

  // Two-step Global Optimal Matching:
  if (num_tracks > 0 && num_dets > 0) {
    // Step 1: IoU matching with predicted bounding boxes (highest IoU first)
    std::vector<MatchPair> iou_pairs;
    for (int i = 0; i < num_tracks; ++i) {
      std::string t_lbl = tracks[i].label;
      for (int j = 0; j < num_dets; ++j) {
        if (t_lbl != Rcpp::as<std::string>(labels[j])) continue;
        double iou = compute_iou(
          pred_xmin[i], pred_ymin[i], pred_xmax[i], pred_ymax[i],
          xmin[j], ymin[j], xmax[j], ymax[j]
        );
        if (iou >= min_iou) {
          iou_pairs.push_back({i, j, iou});
        }
      }
    }

    std::sort(iou_pairs.begin(), iou_pairs.end(), [](const MatchPair& a, const MatchPair& b) {
      return a.score > b.score;
    });

    for (const auto& pair : iou_pairs) {
      if (!track_matched[pair.track_idx] && !det_matched[pair.det_idx]) {
        track_matched[pair.track_idx] = true;
        det_matched[pair.det_idx] = true;
        det_to_track[pair.det_idx] = pair.track_idx;
      }
    }

    // Step 2: Distance matching for remaining unmatched pairs (closest distance first)
    std::vector<MatchPair> dist_pairs;
    double effective_max_dist = std::max(max_dist, 150.0);
    for (int i = 0; i < num_tracks; ++i) {
      if (track_matched[i]) continue;
      std::string t_lbl = tracks[i].label;
      for (int j = 0; j < num_dets; ++j) {
        if (det_matched[j]) continue;
        if (t_lbl != Rcpp::as<std::string>(labels[j])) continue;

        double dx = pred_cx[i] - ncx[j];
        double dy = pred_cy[i] - ncy[j];
        double d = std::sqrt(dx * dx + dy * dy);
        if (d <= effective_max_dist) {
          dist_pairs.push_back({i, j, d});
        }
      }
    }

    std::sort(dist_pairs.begin(), dist_pairs.end(), [](const MatchPair& a, const MatchPair& b) {
      return a.score < b.score;
    });

    for (const auto& pair : dist_pairs) {
      if (!track_matched[pair.track_idx] && !det_matched[pair.det_idx]) {
        track_matched[pair.track_idx] = true;
        det_matched[pair.det_idx] = true;
        det_to_track[pair.det_idx] = pair.track_idx;
      }
    }
  }

  // Update matched tracks state and velocity
  for (int j = 0; j < num_dets; ++j) {
    if (det_matched[j]) {
      int i = det_to_track[j];
      double cur_vx = ncx[j] - tracks[i].cx;
      double cur_vy = ncy[j] - tracks[i].cy;
      if (tracks[i].age <= 1) {
        tracks[i].vx = cur_vx;
        tracks[i].vy = cur_vy;
      } else {
        tracks[i].vx = 0.7 * cur_vx + 0.3 * tracks[i].vx;
        tracks[i].vy = 0.7 * cur_vy + 0.3 * tracks[i].vy;
      }

      tracks[i].prev_cx = tracks[i].cx;
      tracks[i].prev_cy = tracks[i].cy;
      tracks[i].cx = ncx[j];
      tracks[i].cy = ncy[j];
      tracks[i].xmin = xmin[j];
      tracks[i].ymin = ymin[j];
      tracks[i].xmax = xmax[j];
      tracks[i].ymax = ymax[j];
      tracks[i].lost = 0;
      tracks[i].age += 1;
      if (tracks[i].count_flash > 0) tracks[i].count_flash -= 1;

      tracks[i].history_x.push_back(ncx[j]);
      tracks[i].history_y.push_back(ncy[j]);
      if ((int)tracks[i].history_x.size() > max_history) {
        tracks[i].history_x.erase(tracks[i].history_x.begin());
        tracks[i].history_y.erase(tracks[i].history_y.begin());
      }
    }
  }

  // Increment lost count for unmatched tracks
  for (int i = 0; i < num_tracks; ++i) {
    if (!track_matched[i]) {
      tracks[i].lost += 1;
      if (tracks[i].count_flash > 0) tracks[i].count_flash -= 1;
    }
  }

  // Parse count line if supplied: c(lx1, ly1, lx2, ly2)
  bool has_line = false;
  double lx1 = 0, ly1 = 0, lx2 = 0, ly2 = 0;
  if (count_line.isNotNull()) {
    Rcpp::NumericVector cl(count_line.get());
    if (cl.size() >= 4) {
      has_line = true;
      lx1 = cl[0]; ly1 = cl[1]; lx2 = cl[2]; ly2 = cl[3];
    }
  }

  // Parse ROI if supplied: c(rx1, ry1, rx2, ry2)
  bool has_roi = false;
  double rx1 = 0, ry1 = 0, rx2 = 0, ry2 = 0;
  if (roi.isNotNull()) {
    Rcpp::NumericVector rbox(roi.get());
    if (rbox.size() >= 4) {
      has_roi = true;
      rx1 = rbox[0]; ry1 = rbox[1]; rx2 = rbox[2]; ry2 = rbox[3];
      if (rx1 > rx2) std::swap(rx1, rx2);
      if (ry1 > ry2) std::swap(ry1, ry2);
    }
  }

  // Check line crossing and ROI entry for matched tracks
  int new_crossings = 0;
  for (size_t i = 0; i < tracks.size(); ++i) {
    if (tracks[i].counted) continue;
    if (tracks[i].lost > 1) continue;

    bool triggered = false;
    std::string trigger_dir = "";

    // A) Check virtual line intersection if count_line is present
    if (has_line) {
      if (tracks[i].lost == 0 && tracks[i].history_x.size() >= 2) {
        int crossing_dir = check_segment_intersection(
          tracks[i].prev_cx, tracks[i].prev_cy,
          tracks[i].cx, tracks[i].cy,
          lx1, ly1, lx2, ly2
        );
        if (crossing_dir != 0) {
          triggered = true;
          trigger_dir = (crossing_dir > 0) ? "A_to_B" : "B_to_A";
        }
      } else if (tracks[i].lost == 1 && tracks[i].age >= 2) {
        // If track just got lost (e.g. exited screen or cropped area right at the count line),
        // check if its predicted trajectory crossed the line
        double pred_x = tracks[i].cx + tracks[i].vx;
        double pred_y = tracks[i].cy + tracks[i].vy;
        int crossing_dir = check_segment_intersection(
          tracks[i].cx, tracks[i].cy,
          pred_x, pred_y,
          lx1, ly1, lx2, ly2
        );
        if (crossing_dir != 0) {
          triggered = true;
          trigger_dir = (crossing_dir > 0) ? "A_to_B" : "B_to_A";
        }
      }
    }

    // B) Or check entry into ROI zone if object crossed ROI boundary (ONLY when no count_line is defined)
    if (!has_line && has_roi && tracks[i].lost == 0 && tracks[i].history_x.size() >= 2) {
      bool prev_in = (tracks[i].prev_cx >= rx1 && tracks[i].prev_cx <= rx2 &&
                      tracks[i].prev_cy >= ry1 && tracks[i].prev_cy <= ry2);
      bool curr_in = (tracks[i].cx >= rx1 && tracks[i].cx <= rx2 &&
                      tracks[i].cy >= ry1 && tracks[i].cy <= ry2);
      if (!prev_in && curr_in) {
        triggered = true;
        trigger_dir = "enter_roi";
      }
    }

    if (triggered && !tracks[i].counted) {
      tracks[i].counted = true;
      tracks[i].count_flash = 12; // flash green for 12 frames
      tracks[i].last_dir = trigger_dir;
      total_count += 1;
      new_crossings += 1;

      std::string lbl = tracks[i].label;
      int cur_cls_count = counts_by_class.containsElementNamed(lbl.c_str()) ?
        Rcpp::as<int>(counts_by_class[lbl]) : 0;
      counts_by_class[lbl] = cur_cls_count + 1;

      // Record crossing event
      Rcpp::List ev = Rcpp::List::create(
        Rcpp::Named("id") = tracks[i].id,
        Rcpp::Named("label") = lbl,
        Rcpp::Named("direction") = tracks[i].last_dir,
        Rcpp::Named("x") = tracks[i].cx,
        Rcpp::Named("y") = tracks[i].cy
      );
      crossing_events.push_back(ev);
    }
  }

  // Create new tracks for unmatched detections
  for (int j = 0; j < num_dets; ++j) {
    if (!det_matched[j]) {
      TrackObj nt;
      nt.id = next_id++;
      nt.label = Rcpp::as<std::string>(labels[j]);
      nt.xmin = xmin[j];
      nt.ymin = ymin[j];
      nt.xmax = xmax[j];
      nt.ymax = ymax[j];
      nt.cx = ncx[j];
      nt.cy = ncy[j];
      nt.prev_cx = nt.cx;
      nt.prev_cy = nt.cy;
      nt.vx = 0.0;
      nt.vy = 0.0;
      nt.history_x.push_back(nt.cx);
      nt.history_y.push_back(nt.cy);
      nt.lost = 0;
      nt.age = 1;
      nt.counted = false;
      nt.count_flash = 0;
      nt.last_dir = "";

      tracks.push_back(nt);
      det_to_track[j] = (int)tracks.size() - 1;
    }
  }

  // Filter out expired tracks for persistence
  std::vector<TrackObj> surviving_tracks;
  for (size_t i = 0; i < tracks.size(); ++i) {
    if (tracks[i].lost <= max_lost) {
      surviving_tracks.push_back(tracks[i]);
    }
  }

  // Output tracked detections STRICTLY corresponding to input detections 0..num_dets-1
  Rcpp::IntegerVector out_track_ids(num_dets);
  Rcpp::NumericVector out_xmin(num_dets);
  Rcpp::NumericVector out_ymin(num_dets);
  Rcpp::NumericVector out_xmax(num_dets);
  Rcpp::NumericVector out_ymax(num_dets);
  Rcpp::CharacterVector out_labels(num_dets);
  Rcpp::LogicalVector out_counted(num_dets);
  Rcpp::IntegerVector out_flash(num_dets);
  Rcpp::LogicalVector out_in_zone(num_dets);

  Rcpp::List out_history_x(num_dets);
  Rcpp::List out_history_y(num_dets);

  for (int j = 0; j < num_dets; ++j) {
    int t_idx = det_to_track[j];
    out_track_ids[j] = tracks[t_idx].id;
    out_xmin[j] = xmin[j];
    out_ymin[j] = ymin[j];
    out_xmax[j] = xmax[j];
    out_ymax[j] = ymax[j];
    out_labels[j] = tracks[t_idx].label;
    out_counted[j] = tracks[t_idx].counted;
    out_flash[j] = tracks[t_idx].count_flash;

    bool is_in = true;
    if (has_roi) {
      is_in = (ncx[j] >= rx1 && ncx[j] <= rx2 && ncy[j] >= ry1 && ncy[j] <= ry2);
    }
    out_in_zone[j] = is_in;

    out_history_x[j] = Rcpp::wrap(tracks[t_idx].history_x);
    out_history_y[j] = Rcpp::wrap(tracks[t_idx].history_y);
  }

  // Serialize surviving tracks for state persistence
  Rcpp::List new_tracks_list(surviving_tracks.size());
  for (size_t i = 0; i < surviving_tracks.size(); ++i) {
    new_tracks_list[i] = Rcpp::List::create(
      Rcpp::Named("id") = surviving_tracks[i].id,
      Rcpp::Named("label") = surviving_tracks[i].label,
      Rcpp::Named("xmin") = surviving_tracks[i].xmin,
      Rcpp::Named("ymin") = surviving_tracks[i].ymin,
      Rcpp::Named("xmax") = surviving_tracks[i].xmax,
      Rcpp::Named("ymax") = surviving_tracks[i].ymax,
      Rcpp::Named("cx") = surviving_tracks[i].cx,
      Rcpp::Named("cy") = surviving_tracks[i].cy,
      Rcpp::Named("prev_cx") = surviving_tracks[i].prev_cx,
      Rcpp::Named("prev_cy") = surviving_tracks[i].prev_cy,
      Rcpp::Named("vx") = surviving_tracks[i].vx,
      Rcpp::Named("vy") = surviving_tracks[i].vy,
      Rcpp::Named("history_x") = Rcpp::wrap(surviving_tracks[i].history_x),
      Rcpp::Named("history_y") = Rcpp::wrap(surviving_tracks[i].history_y),
      Rcpp::Named("lost") = surviving_tracks[i].lost,
      Rcpp::Named("age") = surviving_tracks[i].age,
      Rcpp::Named("counted") = surviving_tracks[i].counted,
      Rcpp::Named("count_flash") = surviving_tracks[i].count_flash,
      Rcpp::Named("last_dir") = surviving_tracks[i].last_dir
    );
  }

  Rcpp::List new_state = Rcpp::List::create(
    Rcpp::Named("next_id") = next_id,
    Rcpp::Named("total_count") = total_count,
    Rcpp::Named("tracks") = new_tracks_list,
    Rcpp::Named("counts_by_class") = counts_by_class,
    Rcpp::Named("crossing_events") = crossing_events
  );

  int n_in_zone = 0;
  for (size_t i = 0; i < surviving_tracks.size(); ++i) {
    if (surviving_tracks[i].lost == 0) {
      if (has_roi) {
        if (surviving_tracks[i].cx >= rx1 && surviving_tracks[i].cx <= rx2 &&
            surviving_tracks[i].cy >= ry1 && surviving_tracks[i].cy <= ry2) {
          n_in_zone++;
        }
      } else {
        n_in_zone++;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("state") = new_state,
    Rcpp::Named("track_ids") = out_track_ids,
    Rcpp::Named("xmin") = out_xmin,
    Rcpp::Named("ymin") = out_ymin,
    Rcpp::Named("xmax") = out_xmax,
    Rcpp::Named("ymax") = out_ymax,
    Rcpp::Named("labels") = out_labels,
    Rcpp::Named("counted") = out_counted,
    Rcpp::Named("flash") = out_flash,
    Rcpp::Named("history_x") = out_history_x,
    Rcpp::Named("history_y") = out_history_y,
    Rcpp::Named("total_count") = total_count,
    Rcpp::Named("new_crossings") = new_crossings,
    Rcpp::Named("in_zone") = n_in_zone,
    Rcpp::Named("in_zone_vec") = out_in_zone
  );
}
