// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// ==============================================================================
// 1. FAST TILE EXTRACTION & CROP KERNEL (OpenMP)
// ==============================================================================
template <int RTYPE>
SEXP do_crop_tile(SEXP img_sexp, int src_w, int src_h, int nch,
                  int x1, int y1, int tile_w, int tile_h,
                  NumericVector bg_val) {
  typedef typename Rcpp::Vector<RTYPE> VecType;
  typedef typename VecType::stored_type T;

  VecType src_vec(img_sexp);
  const T* src = src_vec.begin();

  VecType out_vec(tile_w * tile_h * nch);
  T* dst = out_vec.begin();

  std::vector<T> bg(nch);
  for (int c = 0; c < nch; ++c) {
    if (RTYPE == RAWSXP) {
      bg[c] = static_cast<T>(std::max(0.0, std::min(255.0, bg_val[c % bg_val.size()])));
    } else if (RTYPE == REALSXP) {
      bg[c] = static_cast<T>(bg_val[c % bg_val.size()]);
    } else {
      bg[c] = static_cast<T>(bg_val[c % bg_val.size()]);
    }
  }

  #pragma omp parallel for collapse(2) schedule(static)
  for (int c = 0; c < nch; ++c) {
    for (int ty = 0; ty < tile_h; ++ty) {
      int sy = (y1 - 1) + ty;
      T bg_c = bg[c];
      size_t dst_row_offset = (static_cast<size_t>(c) * tile_h + ty) * tile_w;

      if (sy < 0 || sy >= src_h) {
        for (int tx = 0; tx < tile_w; ++tx) {
          dst[dst_row_offset + tx] = bg_c;
        }
      } else {
        size_t src_row_offset = (static_cast<size_t>(c) * src_h + sy) * src_w;
        for (int tx = 0; tx < tile_w; ++tx) {
          int sx = (x1 - 1) + tx;
          if (sx >= 0 && sx < src_w) {
            dst[dst_row_offset + tx] = src[src_row_offset + sx];
          } else {
            dst[dst_row_offset + tx] = bg_c;
          }
        }
      }
    }
  }

  IntegerVector out_dims(3);
  out_dims[0] = tile_w;
  out_dims[1] = tile_h;
  out_dims[2] = nch;
  out_vec.attr("dim") = out_dims;
  out_vec.attr("class") = CharacterVector::create("image", "array");

  return out_vec;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP crop_tile_cpp(SEXP img, int x1, int y1, int tile_w, int tile_h, NumericVector bg_val) {
  IntegerVector dims = Rf_getAttrib(img, R_DimSymbol);
  int W = dims[0];
  int H = dims[1];
  int nch = (dims.size() >= 3) ? dims[2] : 1;

  SEXP res = R_NilValue;
  switch (TYPEOF(img)) {
    case RAWSXP:
      res = do_crop_tile<RAWSXP>(img, W, H, nch, x1, y1, tile_w, tile_h, bg_val);
      break;
    case REALSXP:
      res = do_crop_tile<REALSXP>(img, W, H, nch, x1, y1, tile_w, tile_h, bg_val);
      break;
    case INTSXP:
      res = do_crop_tile<INTSXP>(img, W, H, nch, x1, y1, tile_w, tile_h, bg_val);
      break;
    default:
      stop("Unsupported image type for crop_tile_cpp. Must be raw, real, or integer.");
  }

  SEXP cm = Rf_getAttrib(img, Rf_install("colormode"));
  if (!Rf_isNull(cm)) {
    Rf_setAttrib(res, Rf_install("colormode"), cm);
  } else {
    Rf_setAttrib(res, Rf_install("colormode"), wrap(nch >= 3 ? "Color" : "Grayscale"));
  }

  return res;
}


// ==============================================================================
// 2. HIGH-PERFORMANCE BOUNDING BOX CLIPPING TO TILE (YOLO FORMAT)
// ==============================================================================
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix clip_yolo_boxes_cpp(NumericMatrix boxes,
                                  double tile_x1, double tile_y1,
                                  double tile_x2, double tile_y2,
                                  double tile_w, double tile_h,
                                  double min_area_ratio = 0.50,
                                  int mode = 0,
                                  double min_box_size = 2.0) {
  int n_boxes = boxes.nrow();
  if (n_boxes == 0) {
    NumericMatrix empty(0, 5);
    return empty;
  }

  std::vector<double> out_cls;
  std::vector<double> out_xc;
  std::vector<double> out_yc;
  std::vector<double> out_bw;
  std::vector<double> out_bh;

  for (int i = 0; i < n_boxes; ++i) {
    double cid  = boxes(i, 0);
    double b_x1 = boxes(i, 1);
    double b_y1 = boxes(i, 2);
    double b_x2 = boxes(i, 3);
    double b_y2 = boxes(i, 4);

    double orig_w = b_x2 - b_x1;
    double orig_h = b_y2 - b_y1;
    if (orig_w <= 0.0 || orig_h <= 0.0) continue;
    double orig_area = orig_w * orig_h;

    // Intersection with tile boundary
    double ix1 = std::max(b_x1, tile_x1);
    double ix2 = std::min(b_x2, tile_x2);
    double iy1 = std::max(b_y1, tile_y1);
    double iy2 = std::min(b_y2, tile_y2);

    double iw = ix2 - ix1;
    double ih = iy2 - iy1;
    if (iw <= 0.0 || ih <= 0.0) continue;

    double inter_area = iw * ih;
    double ratio = inter_area / orig_area;

    double cx = (b_x1 + b_x2) * 0.5;
    double cy = (b_y1 + b_y2) * 0.5;
    bool center_inside = (cx >= tile_x1 && cx <= tile_x2 && cy >= tile_y1 && cy <= tile_y2);
    bool area_preserved = (ratio >= min_area_ratio);

    bool keep = false;
    if (mode == 0) {
      // center_or_overlap
      keep = center_inside || area_preserved;
    } else if (mode == 1) {
      // center only
      keep = center_inside;
    } else if (mode == 2) {
      // overlap only
      keep = area_preserved;
    }

    if (!keep) continue;

    // Discard tiny clipped slivers
    if (iw < min_box_size || ih < min_box_size) continue;

    // Local normalized coordinates relative to tile
    double lx1 = ix1 - tile_x1;
    double lx2 = ix2 - tile_x1;
    double ly1 = iy1 - tile_y1;
    double ly2 = iy2 - tile_y1;

    double norm_xc = (lx1 + lx2) * 0.5 / tile_w;
    double norm_yc = (ly1 + ly2) * 0.5 / tile_h;
    double norm_w  = (lx2 - lx1) / tile_w;
    double norm_h  = (ly2 - ly1) / tile_h;

    norm_xc = std::max(0.0, std::min(1.0, norm_xc));
    norm_yc = std::max(0.0, std::min(1.0, norm_yc));
    norm_w  = std::max(0.0, std::min(1.0, norm_w));
    norm_h  = std::max(0.0, std::min(1.0, norm_h));

    out_cls.push_back(cid);
    out_xc.push_back(norm_xc);
    out_yc.push_back(norm_yc);
    out_bw.push_back(norm_w);
    out_bh.push_back(norm_h);
  }

  int n_out = out_cls.size();
  NumericMatrix res(n_out, 5);
  for (int i = 0; i < n_out; ++i) {
    res(i, 0) = out_cls[i];
    res(i, 1) = out_xc[i];
    res(i, 2) = out_yc[i];
    res(i, 3) = out_bw[i];
    res(i, 4) = out_bh[i];
  }

  return res;
}


// ==============================================================================
// 3. SUTHERLAND-HODGMAN POLYGON CLIPPING (SEGMENTATION FORMAT)
// ==============================================================================
struct Point2D {
  double x;
  double y;
};

inline double polygon_area_shoelace(const std::vector<Point2D>& poly) {
  int n = poly.size();
  if (n < 3) return 0.0;
  double area = 0.0;
  for (int i = 0; i < n; ++i) {
    int next = (i + 1) % n;
    area += poly[i].x * poly[next].y - poly[next].x * poly[i].y;
  }
  return std::abs(area) * 0.5;
}

inline Point2D polygon_centroid(const std::vector<Point2D>& poly, double area) {
  int n = poly.size();
  if (n < 3 || area <= 1e-9) {
    Point2D p = {0.0, 0.0};
    for (int i = 0; i < n; ++i) { p.x += poly[i].x; p.y += poly[i].y; }
    p.x /= n; p.y /= n;
    return p;
  }
  double cx = 0.0, cy = 0.0;
  for (int i = 0; i < n; ++i) {
    int next = (i + 1) % n;
    double cross = poly[i].x * poly[next].y - poly[next].x * poly[i].y;
    cx += (poly[i].x + poly[next].x) * cross;
    cy += (poly[i].y + poly[next].y) * cross;
  }
  cx /= (6.0 * area);
  cy /= (6.0 * area);
  return {std::abs(cx), std::abs(cy)};
}

inline std::vector<Point2D> clip_polygon_sutherland_hodgman(const std::vector<Point2D>& poly,
                                                            double x1, double y1, double x2, double y2) {
  std::vector<Point2D> output = poly;

  for (int edge = 0; edge < 4; ++edge) {
    if (output.empty()) break;
    std::vector<Point2D> input = output;
    output.clear();

    Point2D s = input.back();
    for (size_t i = 0; i < input.size(); ++i) {
      Point2D p = input[i];

      bool p_inside = false;
      bool s_inside = false;
      Point2D inter = {0.0, 0.0};

      if (edge == 0) { // Left: x >= x1
        p_inside = (p.x >= x1);
        s_inside = (s.x >= x1);
        if (s.x != p.x) {
          inter.x = x1;
          inter.y = s.y + (p.y - s.y) * (x1 - s.x) / (p.x - s.x);
        }
      } else if (edge == 1) { // Right: x <= x2
        p_inside = (p.x <= x2);
        s_inside = (s.x <= x2);
        if (s.x != p.x) {
          inter.x = x2;
          inter.y = s.y + (p.y - s.y) * (x2 - s.x) / (p.x - s.x);
        }
      } else if (edge == 2) { // Top: y >= y1
        p_inside = (p.y >= y1);
        s_inside = (s.y >= y1);
        if (s.y != p.y) {
          inter.y = y1;
          inter.x = s.x + (p.x - s.x) * (y1 - s.y) / (p.y - s.y);
        }
      } else if (edge == 3) { // Bottom: y <= y2
        p_inside = (p.y <= y2);
        s_inside = (s.y <= y2);
        if (s.y != p.y) {
          inter.y = y2;
          inter.x = s.x + (p.x - s.x) * (y2 - s.y) / (p.y - s.y);
        }
      }

      if (p_inside) {
        if (!s_inside) {
          output.push_back(inter);
        }
        output.push_back(p);
      } else if (s_inside) {
        output.push_back(inter);
      }
      s = p;
    }
  }

  return output;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List clip_yolo_polygons_cpp(List polygons,
                            NumericVector class_ids,
                            double tile_x1, double tile_y1,
                            double tile_x2, double tile_y2,
                            double tile_w, double tile_h,
                            double min_area_ratio = 0.50,
                            int mode = 0,
                            int max_vertices = 32) {
  int n_polys = polygons.size();
  std::vector<int> out_cls;
  std::vector<NumericVector> out_coords;

  for (int i = 0; i < n_polys; ++i) {
    NumericVector coords = polygons[i];
    int n_pts = coords.size() / 2;
    if (n_pts < 3) continue;

    std::vector<Point2D> poly(n_pts);
    for (int p = 0; p < n_pts; ++p) {
      poly[p].x = coords[2 * p];
      poly[p].y = coords[2 * p + 1];
    }

    double orig_area = polygon_area_shoelace(poly);
    if (orig_area <= 1.0) continue;

    Point2D centroid = polygon_centroid(poly, orig_area);

    std::vector<Point2D> clipped = clip_polygon_sutherland_hodgman(poly, tile_x1, tile_y1, tile_x2, tile_y2);
    if (clipped.size() < 3) continue;

    double clipped_area = polygon_area_shoelace(clipped);
    double ratio = clipped_area / orig_area;

    bool center_inside = (centroid.x >= tile_x1 && centroid.x <= tile_x2 &&
                          centroid.y >= tile_y1 && centroid.y <= tile_y2);
    bool area_preserved = (ratio >= min_area_ratio);

    bool keep = false;
    if (mode == 0) keep = center_inside || area_preserved;
    else if (mode == 1) keep = center_inside;
    else if (mode == 2) keep = area_preserved;

    if (!keep) continue;

    int cur_n = clipped.size();
    std::vector<Point2D> final_poly;
    if (max_vertices > 0 && cur_n > max_vertices) {
      for (int v = 0; v < max_vertices; ++v) {
        int idx = std::round(v * static_cast<double>(cur_n) / max_vertices);
        if (idx >= cur_n) idx = cur_n - 1;
        final_poly.push_back(clipped[idx]);
      }
    } else {
      final_poly = clipped;
    }

    int out_n = final_poly.size();
    NumericVector norm_vec(out_n * 2);
    for (int v = 0; v < out_n; ++v) {
      double nx = (final_poly[v].x - tile_x1) / tile_w;
      double ny = (final_poly[v].y - tile_y1) / tile_h;
      norm_vec[2 * v]     = std::max(0.0, std::min(1.0, nx));
      norm_vec[2 * v + 1] = std::max(0.0, std::min(1.0, ny));
    }

    out_cls.push_back(static_cast<int>(class_ids[i % class_ids.size()]));
    out_coords.push_back(norm_vec);
  }

  return List::create(
    Named("class_ids") = wrap(out_cls),
    Named("polygons")  = wrap(out_coords)
  );
}
