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

// ---------------------------------------------------------------------------
// Tracker state kept in C++.
// The historical API (update_tracker_cpp) round-trips the whole state through nested
// R lists on every frame, which costs O(tracks x fields + crossing events) per frame
// (crossing events were appended with List::push_back, which copies the list).
// tracker_create_cpp()/tracker_update_cpp() keep the same state behind an external
// pointer and only return the per-frame outputs; update_tracker_cpp() is preserved and
// now runs the same step function through a decode/encode wrapper.
// ---------------------------------------------------------------------------
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

struct CrossingEvent {
  int id;
  std::string label;
  std::string direction;
  double x, y;
};

struct TrackerState {
  int next_id = 1;
  int total_count = 0;
  std::vector<TrackObj> tracks;
  std::vector<std::pair<std::string, int> > counts_by_class;  // insertion order preserved
  std::vector<CrossingEvent> crossing_events;
  double welford_n = 0.0, welford_mean = 0.0, welford_m2 = 0.0;
  std::vector<double> ring_time = std::vector<double>(30, -1.0);
  std::vector<int> ring_count = std::vector<int>(30, 0);
  int ring_head = 0;
  int ring_total = 0;
  std::vector<double> flow_history = std::vector<double>(20, 0.0);
};

static Rcpp::List tracker_counts_to_list(const TrackerState& S) {
  Rcpp::List out(S.counts_by_class.size());
  Rcpp::CharacterVector nm(S.counts_by_class.size());
  for (size_t i = 0; i < S.counts_by_class.size(); ++i) {
    out[i] = S.counts_by_class[i].second;
    nm[i] = S.counts_by_class[i].first;
  }
  if (!S.counts_by_class.empty()) out.attr("names") = nm;
  return out;
}

static Rcpp::List tracker_events_to_list(const TrackerState& S) {
  Rcpp::List out(S.crossing_events.size());
  for (size_t i = 0; i < S.crossing_events.size(); ++i) {
    const CrossingEvent& e = S.crossing_events[i];
    out[i] = Rcpp::List::create(
      Rcpp::Named("id") = e.id,
      Rcpp::Named("label") = e.label,
      Rcpp::Named("direction") = e.direction,
      Rcpp::Named("x") = e.x,
      Rcpp::Named("y") = e.y
    );
  }
  return out;
}

static Rcpp::List tracker_state_to_list(const TrackerState& S) {
  Rcpp::List new_tracks_list(S.tracks.size());
  for (size_t i = 0; i < S.tracks.size(); ++i) {
    const TrackObj& t = S.tracks[i];
    new_tracks_list[i] = Rcpp::List::create(
      Rcpp::Named("id") = t.id,
      Rcpp::Named("label") = t.label,
      Rcpp::Named("class_id") = t.class_id,
      Rcpp::Named("score") = t.score,
      Rcpp::Named("xmin") = t.xmin,
      Rcpp::Named("ymin") = t.ymin,
      Rcpp::Named("xmax") = t.xmax,
      Rcpp::Named("ymax") = t.ymax,
      Rcpp::Named("cx") = t.cx,
      Rcpp::Named("cy") = t.cy,
      Rcpp::Named("prev_cx") = t.prev_cx,
      Rcpp::Named("prev_cy") = t.prev_cy,
      Rcpp::Named("vx") = t.vx,
      Rcpp::Named("vy") = t.vy,
      Rcpp::Named("history_x") = Rcpp::wrap(t.history_x),
      Rcpp::Named("history_y") = Rcpp::wrap(t.history_y),
      Rcpp::Named("lost") = t.lost,
      Rcpp::Named("age") = t.age,
      Rcpp::Named("counted") = t.counted,
      Rcpp::Named("count_flash") = t.count_flash,
      Rcpp::Named("last_dir") = t.last_dir,
      Rcpp::Named("mask_id") = t.mask_id,
      Rcpp::Named("mask_shift_x") = t.mask_shift_x,
      Rcpp::Named("mask_shift_y") = t.mask_shift_y,
      Rcpp::Named("hits") = t.hits,
      Rcpp::Named("consecutive_hits") = t.consecutive_hits,
      Rcpp::Named("confirmed") = t.confirmed
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("next_id") = S.next_id,
    Rcpp::Named("total_count") = S.total_count,
    Rcpp::Named("tracks") = new_tracks_list,
    Rcpp::Named("counts_by_class") = tracker_counts_to_list(S),
    Rcpp::Named("crossing_events") = tracker_events_to_list(S),
    Rcpp::Named("welford_n") = S.welford_n,
    Rcpp::Named("welford_mean") = S.welford_mean,
    Rcpp::Named("welford_m2") = S.welford_m2,
    Rcpp::Named("ring_time") = Rcpp::wrap(S.ring_time),
    Rcpp::Named("ring_count") = Rcpp::wrap(S.ring_count),
    Rcpp::Named("ring_head") = S.ring_head,
    Rcpp::Named("ring_total") = S.ring_total,
    Rcpp::Named("flow_history") = Rcpp::wrap(S.flow_history)
  );
}

static void tracker_state_from_list(Rcpp::List state, TrackerState& S, int min_hits) {
  S = TrackerState();
  if (state.containsElementNamed("next_id")) S.next_id = Rcpp::as<int>(state["next_id"]);
  if (state.containsElementNamed("total_count")) S.total_count = Rcpp::as<int>(state["total_count"]);

  if (state.containsElementNamed("counts_by_class")) {
    Rcpp::List cbc = Rcpp::as<Rcpp::List>(state["counts_by_class"]);
    if (cbc.size() > 0) {
      Rcpp::CharacterVector nm = cbc.names();
      for (R_xlen_t i = 0; i < cbc.size(); ++i) {
        S.counts_by_class.push_back(std::make_pair(Rcpp::as<std::string>(nm[i]), Rcpp::as<int>(cbc[i])));
      }
    }
  }
  if (state.containsElementNamed("crossing_events")) {
    Rcpp::List evs = Rcpp::as<Rcpp::List>(state["crossing_events"]);
    for (R_xlen_t i = 0; i < evs.size(); ++i) {
      Rcpp::List e = evs[i];
      CrossingEvent ce;
      ce.id = Rcpp::as<int>(e["id"]);
      ce.label = Rcpp::as<std::string>(e["label"]);
      ce.direction = Rcpp::as<std::string>(e["direction"]);
      ce.x = Rcpp::as<double>(e["x"]);
      ce.y = Rcpp::as<double>(e["y"]);
      S.crossing_events.push_back(ce);
    }
  }

  S.welford_n = state.containsElementNamed("welford_n") ? Rcpp::as<double>(state["welford_n"]) : 0.0;
  S.welford_mean = state.containsElementNamed("welford_mean") ? Rcpp::as<double>(state["welford_mean"]) : 0.0;
  S.welford_m2 = state.containsElementNamed("welford_m2") ? Rcpp::as<double>(state["welford_m2"]) : 0.0;

  if (state.containsElementNamed("ring_time")) {
    S.ring_time = Rcpp::as<std::vector<double> >(state["ring_time"]);
    S.ring_count = Rcpp::as<std::vector<int> >(state["ring_count"]);
    S.ring_head = Rcpp::as<int>(state["ring_head"]);
    S.ring_total = Rcpp::as<int>(state["ring_total"]);
  }
  if (state.containsElementNamed("flow_history")) {
    S.flow_history = Rcpp::as<std::vector<double> >(state["flow_history"]);
  }

  if (state.containsElementNamed("tracks")) {
    Rcpp::List active_tracks = Rcpp::as<Rcpp::List>(state["tracks"]);
    int num_tracks = active_tracks.size();
    S.tracks.resize(num_tracks);
    std::vector<TrackObj>& tracks = S.tracks;
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
  }
}

static Rcpp::List tracker_step(
    TrackerState& S,
    Rcpp::NumericVector xmin,
    Rcpp::NumericVector ymin,
    Rcpp::NumericVector xmax,
    Rcpp::NumericVector ymax,
    Rcpp::CharacterVector labels,
    Rcpp::NumericVector scores,
    Rcpp::Nullable<Rcpp::NumericVector> count_line,
    Rcpp::Nullable<Rcpp::NumericVector> roi,
    double max_dist,
    double min_iou,
    int max_lost,
    int max_history,
    double current_time,
    int frame_w,
    int frame_h,
    Rcpp::Nullable<Rcpp::IntegerVector> class_ids,
    Rcpp::Nullable<Rcpp::IntegerVector> mask_ids,
    bool interpolate_box,
    bool interpolate_mask,
    int min_hits,
    bool include_state
) {
  int& next_id = S.next_id;
  int& total_count = S.total_count;
  double& welford_n = S.welford_n;
  double& welford_mean = S.welford_mean;
  double& welford_m2 = S.welford_m2;
  std::vector<double>& ring_time = S.ring_time;
  std::vector<int>& ring_count = S.ring_count;
  int& ring_head = S.ring_head;
  int& ring_total = S.ring_total;
  std::vector<double>& flow_history = S.flow_history;
  std::vector<TrackObj>& tracks = S.tracks;

  int num_dets = xmin.size();
  int num_tracks = (int)tracks.size();

  // Labels converted once per frame (not once per track x detection pair)
  std::vector<std::string> det_labels(num_dets);
  for (int j = 0; j < num_dets; ++j) det_labels[j] = Rcpp::as<std::string>(labels[j]);

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
        if (t_lbl != det_labels[j]) continue;
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
        if (t_lbl != det_labels[j]) continue;

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
      bool found_cls = false;
      for (size_t q = 0; q < S.counts_by_class.size(); ++q) {
        if (S.counts_by_class[q].first == lbl) {
          S.counts_by_class[q].second += 1;
          found_cls = true;
          break;
        }
      }
      if (!found_cls) S.counts_by_class.push_back(std::make_pair(lbl, 1));

      // Record crossing event (C++ vector: O(1) amortised, exported to R only on demand)
      CrossingEvent ev;
      ev.id = tracks[i].id;
      ev.label = lbl;
      ev.direction = tracks[i].last_dir;
      ev.x = tracks[i].cx;
      ev.y = tracks[i].cy;
      S.crossing_events.push_back(ev);
    }
  }

  // Create new tracks for unmatched detections
  for (int j = 0; j < num_dets; ++j) {
    if (!det_matched[j]) {
      TrackObj nt;
      nt.id = next_id++;
      nt.label = det_labels[j];
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


  // Persist only surviving tracks (lost <= max_lost) for the next frame.
  S.tracks.swap(surviving_tracks);

  Rcpp::List res = Rcpp::List::create(
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
    Rcpp::Named("flow_history") = Rcpp::wrap(flow_history),
    Rcpp::Named("vel_mean") = vel_mean,
    Rcpp::Named("vel_jit") = vel_jit,
    Rcpp::Named("welford_mean") = welford_mean,
    Rcpp::Named("welford_cv") = area_cv,
    Rcpp::Named("roi_occ") = roi_occ
  );
  if (include_state) {
    // Same element order as the historical API: "state" first.
    Rcpp::List with_state(res.size() + 1);
    Rcpp::CharacterVector nm_old = res.names();
    Rcpp::CharacterVector nm(res.size() + 1);
    with_state[0] = tracker_state_to_list(S);
    nm[0] = "state";
    for (R_xlen_t k = 0; k < res.size(); ++k) {
      with_state[k + 1] = res[k];
      nm[k + 1] = nm_old[k];
    }
    with_state.attr("names") = nm;
    return with_state;
  }
  return res;
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
  TrackerState S;
  tracker_state_from_list(state, S, min_hits);
  return tracker_step(S, xmin, ymin, xmax, ymax, labels, scores, count_line, roi, max_dist, min_iou,
                      max_lost, max_history, current_time, frame_w, frame_h, class_ids, mask_ids,
                      interpolate_box, interpolate_mask, min_hits, true);
}

static TrackerState* tracker_from_xptr(SEXP tracker) {
  if (TYPEOF(tracker) != EXTPTRSXP) Rcpp::stop("tracker must be created with tracker_create_cpp().");
  TrackerState* S = (TrackerState*)R_ExternalPtrAddr(tracker);
  if (S == NULL) Rcpp::stop("tracker pointer is no longer valid (was it saved and reloaded?).");
  return S;
}

// [[Rcpp::export]]
SEXP tracker_create_cpp() {
  Rcpp::XPtr<TrackerState> ptr(new TrackerState(), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List tracker_update_cpp(
    SEXP tracker,
    Rcpp::NumericVector xmin,
    Rcpp::NumericVector ymin,
    Rcpp::NumericVector xmax,
    Rcpp::NumericVector ymax,
    Rcpp::CharacterVector labels,
    Rcpp::NumericVector scores,
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
  TrackerState* S = tracker_from_xptr(tracker);
  return tracker_step(*S, xmin, ymin, xmax, ymax, labels, scores, count_line, roi, max_dist, min_iou,
                      max_lost, max_history, current_time, frame_w, frame_h, class_ids, mask_ids,
                      interpolate_box, interpolate_mask, min_hits, false);
}

// [[Rcpp::export]]
Rcpp::List tracker_velocities_cpp(SEXP tracker) {
  TrackerState* S = tracker_from_xptr(tracker);
  size_t n = S->tracks.size();
  Rcpp::IntegerVector id(n);
  Rcpp::NumericVector vx(n), vy(n);
  for (size_t i = 0; i < n; ++i) {
    id[i] = S->tracks[i].id;
    vx[i] = S->tracks[i].vx;
    vy[i] = S->tracks[i].vy;
  }
  return Rcpp::List::create(Rcpp::Named("id") = id, Rcpp::Named("vx") = vx, Rcpp::Named("vy") = vy);
}

// [[Rcpp::export]]
int tracker_n_tracks_cpp(SEXP tracker) {
  return (int)tracker_from_xptr(tracker)->tracks.size();
}

// [[Rcpp::export]]
Rcpp::List tracker_export_state_cpp(SEXP tracker) {
  return tracker_state_to_list(*tracker_from_xptr(tracker));
}

// [[Rcpp::export]]
void tracker_import_state_cpp(SEXP tracker, Rcpp::List state, int min_hits = 1) {
  tracker_state_from_list(state, *tracker_from_xptr(tracker), min_hits);
}

// [[Rcpp::export]]
void tracker_reset_cpp(SEXP tracker) {
  *tracker_from_xptr(tracker) = TrackerState();
}

// ---------------------------------------------------------------------------
// Nested-box filter used by video inference (replaces the per-frame O(n^2) R loop
// `.filter_nested_mat`). Box i is dropped when a smaller, confident box lies (mostly)
// inside it: child c qualifies when score[c] >= conf_floor, area[c] > 0,
// area[c] < area[i] * max_area_ratio and the boxes overlap strictly; i is removed if
// intersection / area[c] >= containment_thresh for any child. `class_ids` is accepted
// for API symmetry and ignored, exactly like the R implementation. When `group` is
// supplied only boxes of the same group (e.g. the same SAHI tile) are compared.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::LogicalVector filter_nested_boxes_cpp(
    Rcpp::NumericMatrix boxes,
    Rcpp::NumericVector scores,
    Rcpp::Nullable<Rcpp::IntegerVector> class_ids = R_NilValue,
    double containment_thresh = 0.85,
    double max_area_ratio = 0.70,
    double conf_floor = 0.50,
    Rcpp::Nullable<Rcpp::IntegerVector> group = R_NilValue
) {
  Rcpp::IntegerVector grp;
  if (group.isNotNull()) grp = Rcpp::IntegerVector(group.get());
  const int n = boxes.nrow();
  Rcpp::LogicalVector keep(n, true);
  if (n <= 1 || boxes.ncol() < 4) return keep;
  std::vector<double> x1(n), y1(n), x2(n), y2(n), area(n);
  for (int i = 0; i < n; ++i) {
    x1[i] = boxes(i, 0); y1[i] = boxes(i, 1); x2[i] = boxes(i, 2); y2[i] = boxes(i, 3);
    area[i] = std::max(0.0, x2[i] - x1[i]) * std::max(0.0, y2[i] - y1[i]);
  }
  std::vector<int> valid;
  valid.reserve(n);
  for (int j = 0; j < n; ++j) {
    double s = (j < scores.size()) ? scores[j] : NA_REAL;
    if (!ISNAN(s) && s >= conf_floor && area[j] > 0.0) valid.push_back(j);
  }
  if (valid.empty()) return keep;
  const bool use_group = (grp.size() == n);
  std::vector<char> k(n, 1);
  #pragma omp parallel for schedule(dynamic, 16) if(n > 256)
  for (int i = 0; i < n; ++i) {
    double ai = area[i];
    if (ai <= 0.0) continue;
    for (size_t q = 0; q < valid.size(); ++q) {
      int c = valid[q];
      if (c == i) continue;
      if (use_group && grp[c] != grp[i]) continue;
      if (!(area[c] < ai * max_area_ratio)) continue;
      if (!(x2[c] > x1[i] && x1[c] < x2[i] && y2[c] > y1[i] && y1[c] < y2[i])) continue;
      double iw = std::min(x2[i], x2[c]) - std::max(x1[i], x1[c]);
      double ih = std::min(y2[i], y2[c]) - std::max(y1[i], y1[c]);
      double inter = std::max(0.0, iw) * std::max(0.0, ih);
      if (inter / area[c] >= containment_thresh) {
        k[i] = 0;
        break;
      }
    }
  }
  for (int i = 0; i < n; ++i) keep[i] = (k[i] != 0);
  (void)class_ids;
  return keep;
}
