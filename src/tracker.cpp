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
    int max_history = 30,
    double current_time = -1.0,
    int frame_w = 0,
    int frame_h = 0,
    Rcpp::Nullable<Rcpp::IntegerVector> class_ids = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> mask_ids = R_NilValue,
    bool interpolate_box = true,
    bool interpolate_mask = true,
    int min_hits = 1
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

  // Welford accumulators for online mean and variance of counted object area (W * H)
  double welford_n = state.containsElementNamed("welford_n") ? Rcpp::as<double>(state["welford_n"]) : 0.0;
  double welford_mean = state.containsElementNamed("welford_mean") ? Rcpp::as<double>(state["welford_mean"]) : 0.0;
  double welford_m2 = state.containsElementNamed("welford_m2") ? Rcpp::as<double>(state["welford_m2"]) : 0.0;

  // 30-frame ring buffer for rolling flow rate
  Rcpp::NumericVector ring_time;
  Rcpp::IntegerVector ring_count;
  int ring_head = 0;
  int ring_total = 0;
  if (state.containsElementNamed("ring_time")) {
    ring_time = Rcpp::as<Rcpp::NumericVector>(state["ring_time"]);
    ring_count = Rcpp::as<Rcpp::IntegerVector>(state["ring_count"]);
    ring_head = Rcpp::as<int>(state["ring_head"]);
    ring_total = Rcpp::as<int>(state["ring_total"]);
  } else {
    ring_time = Rcpp::NumericVector(30, -1.0);
    ring_count = Rcpp::IntegerVector(30, 0);
    ring_head = 0;
    ring_total = 0;
  }

  // 20-step flow rate history for sparkline
  Rcpp::NumericVector flow_history;
  if (state.containsElementNamed("flow_history")) {
    flow_history = Rcpp::as<Rcpp::NumericVector>(state["flow_history"]);
  } else {
    flow_history = Rcpp::NumericVector(20, 0.0);
  }

  int num_dets = xmin.size();
  int num_tracks = active_tracks.size();

  // Internal structure for C++ tracking
  struct TrackObj {
    int id;
    std::string label;
    int class_id;
    double score;
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
    int mask_id;
    double mask_shift_x;
    double mask_shift_y;
    int hits;
    int consecutive_hits;
    bool confirmed;
  };

  std::vector<TrackObj> tracks(num_tracks);
  for (int i = 0; i < num_tracks; ++i) {
    Rcpp::List t = active_tracks[i];
    tracks[i].id = Rcpp::as<int>(t["id"]);
    tracks[i].label = Rcpp::as<std::string>(t["label"]);
    tracks[i].class_id = t.containsElementNamed("class_id") ? Rcpp::as<int>(t["class_id"]) : 0;
    tracks[i].score = t.containsElementNamed("score") ? Rcpp::as<double>(t["score"]) : 0.8;
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
    tracks[i].mask_id = t.containsElementNamed("mask_id") ? Rcpp::as<int>(t["mask_id"]) : -1;
    tracks[i].mask_shift_x = t.containsElementNamed("mask_shift_x") ? Rcpp::as<double>(t["mask_shift_x"]) : 0.0;
    tracks[i].mask_shift_y = t.containsElementNamed("mask_shift_y") ? Rcpp::as<double>(t["mask_shift_y"]) : 0.0;
    tracks[i].hits = t.containsElementNamed("hits") ? Rcpp::as<int>(t["hits"]) : tracks[i].age;
    tracks[i].consecutive_hits = t.containsElementNamed("consecutive_hits") ? Rcpp::as<int>(t["consecutive_hits"]) : (tracks[i].lost == 0 ? 1 : 0);
    tracks[i].confirmed = t.containsElementNamed("confirmed") ? Rcpp::as<bool>(t["confirmed"]) : (min_hits <= 1);
  }

  // Precompute detection centroids
  bool has_cids = class_ids.isNotNull();
  Rcpp::IntegerVector cids;
  if (has_cids) cids = Rcpp::as<Rcpp::IntegerVector>(class_ids);

  bool has_mids = mask_ids.isNotNull();
  Rcpp::IntegerVector mids;
  if (has_mids) mids = Rcpp::as<Rcpp::IntegerVector>(mask_ids);

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
      tracks[i].hits += 1;
      tracks[i].consecutive_hits += 1;
      if (tracks[i].hits >= min_hits || tracks[i].consecutive_hits >= min_hits) {
        tracks[i].confirmed = true;
      }
      tracks[i].score = (!Rcpp::NumericVector::is_na(scores[j])) ? scores[j] : 0.8;
      tracks[i].class_id = (has_cids && j < cids.size()) ? cids[j] : 0;
      tracks[i].mask_id = (has_mids && j < mids.size()) ? mids[j] : (j + 1);
      tracks[i].mask_shift_x = 0.0;
      tracks[i].mask_shift_y = 0.0;
      if (tracks[i].count_flash > 0) tracks[i].count_flash -= 1;

      tracks[i].history_x.push_back(ncx[j]);
      tracks[i].history_y.push_back(ncy[j]);
      if ((int)tracks[i].history_x.size() > max_history) {
        tracks[i].history_x.erase(tracks[i].history_x.begin());
        tracks[i].history_y.erase(tracks[i].history_y.begin());
      }
    }
  }

  // Increment lost count for unmatched tracks and linearly propagate position
  for (int i = 0; i < num_tracks; ++i) {
    if (!track_matched[i]) {
      tracks[i].lost += 1;
      tracks[i].consecutive_hits = 0;
      if (tracks[i].count_flash > 0) tracks[i].count_flash -= 1;
      tracks[i].prev_cx = tracks[i].cx;
      tracks[i].prev_cy = tracks[i].cy;
      tracks[i].cx += tracks[i].vx;
      tracks[i].cy += tracks[i].vy;
      tracks[i].xmin += tracks[i].vx;
      tracks[i].xmax += tracks[i].vx;
      tracks[i].ymin += tracks[i].vy;
      tracks[i].ymax += tracks[i].vy;
      tracks[i].history_x.push_back(tracks[i].cx);
      tracks[i].history_y.push_back(tracks[i].cy);
      if ((int)tracks[i].history_x.size() > max_history) {
        tracks[i].history_x.erase(tracks[i].history_x.begin());
        tracks[i].history_y.erase(tracks[i].history_y.begin());
      }
      tracks[i].mask_shift_x += tracks[i].vx;
      tracks[i].mask_shift_y += tracks[i].vy;
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

    if (triggered && tracks[i].confirmed && !tracks[i].counted) {
      tracks[i].counted = true;
      tracks[i].count_flash = 12; // flash green for 12 frames
      tracks[i].last_dir = trigger_dir;
      total_count += 1;
      new_crossings += 1;

      // Update Welford accumulator for counted object area (W * H)
      double bw = std::max(0.0, tracks[i].xmax - tracks[i].xmin);
      double bh = std::max(0.0, tracks[i].ymax - tracks[i].ymin);
      double obj_area = bw * bh;
      welford_n += 1.0;
      double delta = obj_area - welford_mean;
      welford_mean += delta / welford_n;
      double delta2 = obj_area - welford_mean;
      welford_m2 += delta * delta2;

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
      nt.hits = 1;
      nt.consecutive_hits = 1;
      nt.confirmed = (min_hits <= 1);
      nt.counted = false;
      nt.count_flash = 0;
      nt.last_dir = "";

      nt.class_id = (has_cids && j < cids.size()) ? cids[j] : 0;
      nt.score = (!Rcpp::NumericVector::is_na(scores[j])) ? scores[j] : 0.8;
      nt.mask_id = (has_mids && j < mids.size()) ? mids[j] : (j + 1);
      nt.mask_shift_x = 0.0;
      nt.mask_shift_y = 0.0;

      tracks.push_back(nt);
      det_to_track[j] = (int)tracks.size() - 1;
    }
  }

  // Filter out expired tracks for persistence
  std::vector<TrackObj> surviving_tracks;
  std::vector<int> interp_track_indices;
  for (size_t i = 0; i < tracks.size(); ++i) {
    if (tracks[i].lost <= max_lost) {
      surviving_tracks.push_back(tracks[i]);
      if (interpolate_box && tracks[i].lost > 0 && tracks[i].confirmed) {
        interp_track_indices.push_back((int)i);
      }
    }
  }

  // Which matched detections belong to confirmed tracks?
  std::vector<int> confirmed_det_indices;
  for (int j = 0; j < num_dets; ++j) {
    int t_idx = det_to_track[j];
    if (t_idx >= 0 && tracks[t_idx].confirmed) {
      confirmed_det_indices.push_back(j);
    }
  }

  int total_out = (int)confirmed_det_indices.size() + (int)interp_track_indices.size();

  Rcpp::IntegerVector out_track_ids(total_out);
  Rcpp::NumericVector out_xmin(total_out);
  Rcpp::NumericVector out_ymin(total_out);
  Rcpp::NumericVector out_xmax(total_out);
  Rcpp::NumericVector out_ymax(total_out);
  Rcpp::CharacterVector out_labels(total_out);
  Rcpp::NumericVector out_scores(total_out);
  Rcpp::IntegerVector out_class_ids(total_out);
  Rcpp::LogicalVector out_counted(total_out);
  Rcpp::IntegerVector out_flash(total_out);
  Rcpp::LogicalVector out_in_zone(total_out);
  Rcpp::IntegerVector out_lost(total_out);
  Rcpp::IntegerVector out_mask_ids(total_out);
  Rcpp::NumericVector out_mask_shift_x(total_out);
  Rcpp::NumericVector out_mask_shift_y(total_out);

  Rcpp::List out_history_x(total_out);
  Rcpp::List out_history_y(total_out);

  // 1. Output matched detections for confirmed tracks
  for (size_t k = 0; k < confirmed_det_indices.size(); ++k) {
    int j = confirmed_det_indices[k];
    int t_idx = det_to_track[j];
    out_track_ids[k] = tracks[t_idx].id;
    out_xmin[k] = xmin[j];
    out_ymin[k] = ymin[j];
    out_xmax[k] = xmax[j];
    out_ymax[k] = ymax[j];
    out_labels[k] = tracks[t_idx].label;
    out_scores[k] = tracks[t_idx].score;
    out_class_ids[k] = tracks[t_idx].class_id;
    out_counted[k] = tracks[t_idx].counted;
    out_flash[k] = tracks[t_idx].count_flash;
    out_lost[k] = 0;
    out_mask_ids[k] = tracks[t_idx].mask_id;
    out_mask_shift_x[k] = 0.0;
    out_mask_shift_y[k] = 0.0;

    bool is_in = true;
    if (has_roi) {
      is_in = (ncx[j] >= rx1 && ncx[j] <= rx2 && ncy[j] >= ry1 && ncy[j] <= ry2);
    }
    out_in_zone[k] = is_in;

    out_history_x[k] = Rcpp::wrap(tracks[t_idx].history_x);
    out_history_y[k] = Rcpp::wrap(tracks[t_idx].history_y);
  }

  // 2. Output interpolated tracks (lost <= max_lost and confirmed)
  for (size_t k = 0; k < interp_track_indices.size(); ++k) {
    int t_idx = interp_track_indices[k];
    size_t out_idx = confirmed_det_indices.size() + k;

    out_track_ids[out_idx] = tracks[t_idx].id;
    out_xmin[out_idx] = tracks[t_idx].xmin;
    out_ymin[out_idx] = tracks[t_idx].ymin;
    out_xmax[out_idx] = tracks[t_idx].xmax;
    out_ymax[out_idx] = tracks[t_idx].ymax;
    out_labels[out_idx] = tracks[t_idx].label;
    out_scores[out_idx] = tracks[t_idx].score;
    out_class_ids[out_idx] = tracks[t_idx].class_id;
    out_counted[out_idx] = tracks[t_idx].counted;
    out_flash[out_idx] = tracks[t_idx].count_flash;
    out_lost[out_idx] = tracks[t_idx].lost;
    out_mask_ids[out_idx] = interpolate_mask ? tracks[t_idx].mask_id : -1;
    out_mask_shift_x[out_idx] = tracks[t_idx].mask_shift_x;
    out_mask_shift_y[out_idx] = tracks[t_idx].mask_shift_y;

    bool is_in = true;
    if (has_roi) {
      is_in = (tracks[t_idx].cx >= rx1 && tracks[t_idx].cx <= rx2 &&
               tracks[t_idx].cy >= ry1 && tracks[t_idx].cy <= ry2);
    }
    out_in_zone[out_idx] = is_in;

    out_history_x[out_idx] = Rcpp::wrap(tracks[t_idx].history_x);
    out_history_y[out_idx] = Rcpp::wrap(tracks[t_idx].history_y);
  }

  // Serialize surviving tracks for state persistence
  Rcpp::List new_tracks_list(surviving_tracks.size());
  for (size_t i = 0; i < surviving_tracks.size(); ++i) {
    new_tracks_list[i] = Rcpp::List::create(
      Rcpp::Named("id") = surviving_tracks[i].id,
      Rcpp::Named("label") = surviving_tracks[i].label,
      Rcpp::Named("class_id") = surviving_tracks[i].class_id,
      Rcpp::Named("score") = surviving_tracks[i].score,
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
      Rcpp::Named("last_dir") = surviving_tracks[i].last_dir,
      Rcpp::Named("mask_id") = surviving_tracks[i].mask_id,
      Rcpp::Named("mask_shift_x") = surviving_tracks[i].mask_shift_x,
      Rcpp::Named("mask_shift_y") = surviving_tracks[i].mask_shift_y,
      Rcpp::Named("hits") = surviving_tracks[i].hits,
      Rcpp::Named("consecutive_hits") = surviving_tracks[i].consecutive_hits,
      Rcpp::Named("confirmed") = surviving_tracks[i].confirmed
    );
  }

  // Compute Welford summary metrics
  double area_var = (welford_n > 1.0) ? (welford_m2 / (welford_n - 1.0)) : 0.0;
  double area_sd = std::sqrt(area_var);
  double area_cv = (welford_mean > 1e-4) ? (area_sd / welford_mean * 100.0) : 0.0;

  // Rolling Flow Rate via 30-frame ring buffer
  double flow_rate = 0.0;
  if (current_time >= 0.0 && ring_time.size() >= 30 && ring_count.size() >= 30) {
    ring_total -= ring_count[ring_head];
    ring_time[ring_head] = current_time;
    ring_count[ring_head] = new_crossings;
    ring_total += new_crossings;
    ring_head = (ring_head + 1) % 30;

    double oldest_t = ring_time[ring_head];
    double dt = (oldest_t >= 0.0) ? (current_time - oldest_t) : 0.0;
    flow_rate = (dt > 0.05) ? ((double)ring_total / dt) : 0.0;
  }

  // Update 20-step flow history for sparkline
  if (flow_history.size() >= 20) {
    for (int k = 0; k < 19; ++k) {
      flow_history[k] = flow_history[k + 1];
    }
    flow_history[19] = flow_rate;
  }

  // Single-pass aggregations over active tracks in ROI (Kinematics, Occupancy, Active count)
  int n_in_zone = 0;
  int n_active_roi = 0;
  double vel_sum = 0.0;
  double vel_sq_sum = 0.0;
  int vel_n = 0;
  double roi_box_area_sum = 0.0;

  for (size_t i = 0; i < surviving_tracks.size(); ++i) {
    bool in_roi = true;
    if (has_roi) {
      in_roi = (surviving_tracks[i].cx >= rx1 && surviving_tracks[i].cx <= rx2 &&
                surviving_tracks[i].cy >= ry1 && surviving_tracks[i].cy <= ry2);
    }
    bool is_vis = surviving_tracks[i].confirmed && ((surviving_tracks[i].lost == 0) || (interpolate_box && surviving_tracks[i].lost <= max_lost));
    if (is_vis && in_roi) {
      n_in_zone++;
      n_active_roi++;

      double spd = std::sqrt(surviving_tracks[i].vx * surviving_tracks[i].vx +
                             surviving_tracks[i].vy * surviving_tracks[i].vy);
      vel_sum += spd;
      vel_sq_sum += spd * spd;
      vel_n++;

      double bw = std::max(0.0, surviving_tracks[i].xmax - surviving_tracks[i].xmin);
      double bh = std::max(0.0, surviving_tracks[i].ymax - surviving_tracks[i].ymin);
      roi_box_area_sum += bw * bh;
    }
  }

  double vel_mean = (vel_n > 0) ? (vel_sum / (double)vel_n) : 0.0;
  double vel_var = (vel_n > 1) ? std::max(0.0, (vel_sq_sum - (vel_sum * vel_sum) / (double)vel_n) / (double)(vel_n - 1)) : 0.0;
  double vel_jit = std::sqrt(vel_var);

  // Confidence metrics over all visible detections/tracks in ROI
  double conf_sum = 0.0;
  double conf_sq_sum = 0.0;
  double conf_min = (total_out > 0) ? 1.0 : 0.0;
  int conf_n = 0;
  for (int j = 0; j < total_out; ++j) {
    if (out_in_zone[j] && !Rcpp::NumericVector::is_na(out_scores[j])) {
      double sc = out_scores[j];
      conf_sum += sc;
      conf_sq_sum += sc * sc;
      if (sc < conf_min) conf_min = sc;
      conf_n++;
    }
  }
  if (conf_n == 0) conf_min = 0.0;
  double conf_mean = (conf_n > 0) ? (conf_sum / (double)conf_n) : 0.0;
  double conf_var = (conf_n > 1) ? std::max(0.0, (conf_sq_sum - (conf_sum * conf_sum) / (double)conf_n) / (double)(conf_n - 1)) : 0.0;
  double conf_sd = std::sqrt(conf_var);

  // ROI Occupancy percentage (fast bounding box sum approximation)
  double a_roi = 1.0;
  if (has_roi) {
    a_roi = std::max(1.0, (rx2 - rx1) * (ry2 - ry1));
  } else if (frame_w > 0 && frame_h > 0) {
    a_roi = (double)(frame_w * frame_h);
  }
  double roi_occ = std::min(100.0, (roi_box_area_sum / a_roi) * 100.0);

  Rcpp::List new_state = Rcpp::List::create(
    Rcpp::Named("next_id") = next_id,
    Rcpp::Named("total_count") = total_count,
    Rcpp::Named("tracks") = new_tracks_list,
    Rcpp::Named("counts_by_class") = counts_by_class,
    Rcpp::Named("crossing_events") = crossing_events,
    Rcpp::Named("welford_n") = welford_n,
    Rcpp::Named("welford_mean") = welford_mean,
    Rcpp::Named("welford_m2") = welford_m2,
    Rcpp::Named("ring_time") = ring_time,
    Rcpp::Named("ring_count") = ring_count,
    Rcpp::Named("ring_head") = ring_head,
    Rcpp::Named("ring_total") = ring_total,
    Rcpp::Named("flow_history") = flow_history
  );

  return Rcpp::List::create(
    Rcpp::Named("state") = new_state,
    Rcpp::Named("track_ids") = out_track_ids,
    Rcpp::Named("xmin") = out_xmin,
    Rcpp::Named("ymin") = out_ymin,
    Rcpp::Named("xmax") = out_xmax,
    Rcpp::Named("ymax") = out_ymax,
    Rcpp::Named("labels") = out_labels,
    Rcpp::Named("scores") = out_scores,
    Rcpp::Named("class_ids") = out_class_ids,
    Rcpp::Named("counted") = out_counted,
    Rcpp::Named("flash") = out_flash,
    Rcpp::Named("history_x") = out_history_x,
    Rcpp::Named("history_y") = out_history_y,
    Rcpp::Named("lost") = out_lost,
    Rcpp::Named("mask_ids") = out_mask_ids,
    Rcpp::Named("mask_shift_x") = out_mask_shift_x,
    Rcpp::Named("mask_shift_y") = out_mask_shift_y,
    Rcpp::Named("total_count") = total_count,
    Rcpp::Named("new_crossings") = new_crossings,
    Rcpp::Named("in_zone") = n_in_zone,
    Rcpp::Named("in_zone_vec") = out_in_zone,
    Rcpp::Named("active") = n_active_roi,
    Rcpp::Named("conf_mean") = conf_mean,
    Rcpp::Named("conf_sd") = conf_sd,
    Rcpp::Named("conf_min") = conf_min,
    Rcpp::Named("flow") = flow_rate,
    Rcpp::Named("flow_history") = flow_history,
    Rcpp::Named("vel_mean") = vel_mean,
    Rcpp::Named("vel_jit") = vel_jit,
    Rcpp::Named("welford_mean") = welford_mean,
    Rcpp::Named("welford_cv") = area_cv,
    Rcpp::Named("roi_occ") = roi_occ
  );
}
