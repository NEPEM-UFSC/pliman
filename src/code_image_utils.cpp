#include <RcppArmadillo.h>
#include <queue>
#include <stack>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// adapted from https://en.wikipedia.org/wiki/Sobel_operator#MATLAB_implementation
// [[Rcpp::export]]
NumericMatrix sobel_help(NumericMatrix A) {
  NumericMatrix Gx(3, 3);
  NumericMatrix Gy(3, 3);
  Gx(0, 0) = -1; Gx(0, 1) = 0; Gx(0, 2) = 1;
  Gx(1, 0) = -2; Gx(1, 1) = 0; Gx(1, 2) = 2;
  Gx(2, 0) = -1; Gx(2, 1) = 0; Gx(2, 2) = 1;
  Gy(0, 0) = -1; Gy(0, 1) = -2; Gy(0, 2) = -1;
  Gy(1, 0) = 0; Gy(1, 1) = 0; Gy(1, 2) = 0;
  Gy(2, 0) = 1; Gy(2, 1) = 2; Gy(2, 2) = 1;

  int rows = A.nrow();
  int columns = A.ncol();
  NumericMatrix mag(rows, columns);

  for (int i = 0; i < rows - 2; i++) {
    for (int j = 0; j < columns - 2; j++) {
      double S1 = 0;
      double S2 = 0;

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          S1 += Gx(k, l) * A(i + k, j + l);
          S2 += Gy(k, l) * A(i + k, j + l);
        }
      }
      mag(i + 1, j + 1) = sqrt(S1 * S1 + S2 * S2);
    }
  }
  return mag;
}

// [[Rcpp::export]]
NumericMatrix help_edge_thinning(NumericMatrix img) {
  int rows = img.nrow();
  int cols = img.ncol();
  NumericMatrix thinned(rows, cols);

  for (int i = 1; i < rows - 1; i++) {
    for (int j = 1; j < cols - 1; j++) {
      int p2 = img(i-1, j);
      int p3 = img(i-1, j+1);
      int p4 = img(i, j+1);
      int p5 = img(i+1, j+1);
      int p6 = img(i+1, j);
      int p7 = img(i+1, j-1);
      int p8 = img(i, j-1);
      int p9 = img(i-1, j-1);

      int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
        (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
        (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
        (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);

      int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
      int m1 = (p2 * p4 * p8);
      int m2 = (p4 * p6 * p8);

      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0) {
        thinned(i,j) = 0;
      } else {
        thinned(i,j) = img(i,j);
      }
    }
  }
  return thinned;
}

// GET THE COORDINATES OF A BOUNDING BOX OF A BINARY IMAGE
// [[Rcpp::export]]
IntegerVector bounding_box(LogicalMatrix img, int edge) {
  int nrow = img.nrow();
  int ncol = img.ncol();

  int min_row = nrow;
  int max_row = 0;
  int min_col = ncol;
  int max_col = 0;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (img(i, j)) {
        min_row = std::min(min_row, i);
        max_row = std::max(max_row, i);
        min_col = std::min(min_col, j);
        max_col = std::max(max_col, j);
      }
    }
  }
  min_row = std::max(0, min_row - edge);
  max_row = std::min(nrow - 1, max_row + edge);
  min_col = std::max(0, min_col - edge);
  max_col = std::min(ncol - 1, max_col + edge);
  return IntegerVector::create(min_row, max_row, min_col, max_col);
}

// [[Rcpp::export]]
List isolate_objects5(NumericMatrix img, IntegerMatrix labels) {
  int nrows = labels.nrow(), ncols = labels.ncol();

  IntegerVector unique_labels = sort_unique(labels);
  unique_labels = unique_labels[unique_labels != 0];

  List isolated_objects(unique_labels.length());

  for (int i = 0; i < unique_labels.length(); i++) {
    int id = unique_labels[i];
    int top = nrows, bottom = 0, left = ncols, right = 0;

    for (int j = 0; j < nrows; j++) {
      for (int k = 0; k < ncols; k++) {
        if (labels(j,k) == id) {
          top = std::min(top, j);
          bottom = std::max(bottom, j);
          left = std::min(left, k);
          right = std::max(right, k);
        }
      }
    }

    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;
    NumericMatrix cropped(crop_nrows, crop_ncols);
    for (int j = 0; j < crop_nrows; j++) {
      for (int k = 0; k < crop_ncols; k++) {
        cropped(j,k) = img(top + j, left + k);
      }
    }
    isolated_objects[i] = cropped;
  }

  return isolated_objects;
}

// HELPER FUNCTION TO ISOLATE OBJECTS BASED ON R-G-B and labels
// [[Rcpp::export]]
List help_isolate_object(SEXP R_sexp, SEXP G_sexp, SEXP B_sexp, SEXP labels_sexp, bool remove_bg, int edge) {
  int nrows = Rf_nrows(labels_sexp);
  int ncols = Rf_ncols(labels_sexp);
  if (nrows == 0 || ncols == 0) {
    SEXP dims = Rf_getAttrib(labels_sexp, R_DimSymbol);
    if (!Rf_isNull(dims) && Rf_length(dims) >= 2) {
      nrows = INTEGER(dims)[0];
      ncols = INTEGER(dims)[1];
    }
  }

  const int* p_labels = INTEGER(labels_sexp);
  int npix = nrows * ncols;

  int max_label = 0;
  for (int i = 0; i < npix; i++) {
    if (p_labels[i] > max_label) max_label = p_labels[i];
  }

  if (max_label == 0) return List();

  std::vector<int> b_top(max_label + 1, nrows);
  std::vector<int> b_bottom(max_label + 1, -1);
  std::vector<int> b_left(max_label + 1, ncols);
  std::vector<int> b_right(max_label + 1, -1);
  std::vector<bool> present(max_label + 1, false);

  for (int c = 0; c < ncols; c++) {
    int c_off = c * nrows;
    for (int r = 0; r < nrows; r++) {
      int id = p_labels[c_off + r];
      if (id > 0) {
        present[id] = true;
        if (r < b_top[id]) b_top[id] = r;
        if (r > b_bottom[id]) b_bottom[id] = r;
        if (c < b_left[id]) b_left[id] = c;
        if (c > b_right[id]) b_right[id] = c;
      }
    }
  }

  std::vector<int> valid_ids;
  for (int id = 1; id <= max_label; id++) {
    if (present[id]) valid_ids.push_back(id);
  }

  int num_objs = (int)valid_ids.size();
  List isolated_objects(num_objs);

  bool is_raw = (TYPEOF(R_sexp) == RAWSXP);
  const uint8_t* pR_raw = is_raw ? RAW(R_sexp) : nullptr;
  const uint8_t* pG_raw = is_raw ? RAW(G_sexp) : nullptr;
  const uint8_t* pB_raw = is_raw ? RAW(B_sexp) : nullptr;

  const double* pR_num = !is_raw ? REAL(R_sexp) : nullptr;
  const double* pG_num = !is_raw ? REAL(G_sexp) : nullptr;
  const double* pB_num = !is_raw ? REAL(B_sexp) : nullptr;

  for (int i = 0; i < num_objs; i++) {
    int id = valid_ids[i];
    int top = std::max(0, b_top[id] - edge);
    int bottom = std::min(nrows - 1, b_bottom[id] + edge);
    int left = std::max(0, b_left[id] - edge);
    int right = std::min(ncols - 1, b_right[id] + edge);

    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;

    if (is_raw) {
      RawMatrix croppedR(crop_nrows, crop_ncols);
      RawMatrix croppedG(crop_nrows, crop_ncols);
      RawMatrix croppedB(crop_nrows, crop_ncols);
      uint8_t* p_cR = RAW(croppedR);
      uint8_t* p_cG = RAW(croppedG);
      uint8_t* p_cB = RAW(croppedB);

      for (int c = 0; c < crop_ncols; c++) {
        int src_c = left + c;
        int src_off = src_c * nrows;
        int dst_off = c * crop_nrows;
        for (int r = 0; r < crop_nrows; r++) {
          int src_r = top + r;
          int src_idx = src_off + src_r;
          int dst_idx = dst_off + r;
          if (remove_bg && p_labels[src_idx] != id) {
            p_cR[dst_idx] = 255;
            p_cG[dst_idx] = 255;
            p_cB[dst_idx] = 255;
          } else {
            p_cR[dst_idx] = pR_raw[src_idx];
            p_cG[dst_idx] = pG_raw[src_idx];
            p_cB[dst_idx] = pB_raw[src_idx];
          }
        }
      }
      isolated_objects[i] = List::create(croppedR, croppedG, croppedB);
    } else {
      NumericMatrix croppedR(crop_nrows, crop_ncols);
      NumericMatrix croppedG(crop_nrows, crop_ncols);
      NumericMatrix croppedB(crop_nrows, crop_ncols);
      double* p_cR = REAL(croppedR);
      double* p_cG = REAL(croppedG);
      double* p_cB = REAL(croppedB);

      for (int c = 0; c < crop_ncols; c++) {
        int src_c = left + c;
        int src_off = src_c * nrows;
        int dst_off = c * crop_nrows;
        for (int r = 0; r < crop_nrows; r++) {
          int src_r = top + r;
          int src_idx = src_off + src_r;
          int dst_idx = dst_off + r;
          if (remove_bg && p_labels[src_idx] != id) {
            p_cR[dst_idx] = 1.0;
            p_cG[dst_idx] = 1.0;
            p_cB[dst_idx] = 1.0;
          } else {
            p_cR[dst_idx] = pR_num[src_idx];
            p_cG[dst_idx] = pG_num[src_idx];
            p_cB[dst_idx] = pB_num[src_idx];
          }
        }
      }
      isolated_objects[i] = List::create(croppedR, croppedG, croppedB);
    }
  }

  return isolated_objects;
}

// [[Rcpp::export]]
NumericMatrix help_shp(int rows, int cols, NumericVector dims, double buffer_x, double buffer_y) {
  double xmin = dims[0];
  double xmax = dims[1];
  double ymin = dims[2];
  double ymax = dims[3];
  double xr = xmax - xmin;
  double yr = ymax - ymin;

  double intx = xr / cols;
  double inty = yr / rows;

  NumericMatrix coords(rows * cols * 5, 2);
  int con = 0;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      con++;

      double x_start = xmin + j * intx;
      double x_end = x_start + intx;
      double y_start = ymin + i * inty;
      double y_end = y_start + inty;

      double buffered_x_start = x_start + buffer_x * intx;
      double buffered_x_end = x_end - buffer_x * intx;
      double buffered_y_start = y_start + buffer_y * inty;
      double buffered_y_end = y_end - buffer_y * inty;

      coords((con - 1) * 5, 0) = buffered_x_start;
      coords((con - 1) * 5, 1) = buffered_y_start;
      coords((con - 1) * 5 + 1, 0) = buffered_x_end;
      coords((con - 1) * 5 + 1, 1) = buffered_y_start;
      coords((con - 1) * 5 + 2, 0) = buffered_x_end;
      coords((con - 1) * 5 + 2, 1) = buffered_y_end;
      coords((con - 1) * 5 + 3, 0) = buffered_x_start;
      coords((con - 1) * 5 + 3, 1) = buffered_y_end;
      coords((con - 1) * 5 + 4, 0) = buffered_x_start;
      coords((con - 1) * 5 + 4, 1) = buffered_y_start;
    }
  }
  return coords;
}

// [[Rcpp::export]]
IntegerMatrix helper_guo_hall(IntegerMatrix image) {
  int wid = image.ncol();
  int hgt = image.nrow();
  IntegerMatrix data2 = Rcpp::clone(image);

  auto get = [&](int col, int row) { return image(row, col) != 0; };
  auto clear = [&](int col, int row) { data2(row, col) = 0; };

  IntegerMatrix stepCounter(wid, hgt);

  auto removePixel = [&](int col, int row, bool even) {
    if (!get(col, row)) return 0;
    int p2 = get(col - 1, row);
    int p3 = get(col - 1, row + 1);
    int p4 = get(col, row + 1);
    int p5 = get(col + 1, row + 1);
    int p6 = get(col + 1, row);
    int p7 = get(col + 1, row - 1);
    int p8 = get(col, row - 1);
    int p9 = get(col - 1, row - 1);
    int C = ((!p2) & (p3 | p4)) + ((!p4) & (p5 | p6)) + ((!p6) & (p7 | p8)) + ((!p8) & (p9 | p2));
    if (C != 1) return 0;
    int N1 = (p9 | p2) + (p3 | p4) + (p5 | p6) + (p7 | p8);
    int N2 = (p2 | p3) + (p4 | p5) + (p6 | p7) + (p8 | p9);
    int N = (N1 < N2) ? N1 : N2;
    if (N < 2 || N > 3) return 0;
    int m = even ? ((p6 | p7 | (!p9)) & p8) : ((p2 | p3 | (!p5)) & p4);
    if (m == 0) {
      clear(col, row);
      stepCounter(row, col) = 1;
      return 1;
    }
    return 0;
  };

  bool even = true;
  auto thinStep = [&]() {
    int result = 0;
    for (int row = 1; row < hgt - 1; row++) {
      for (int col = 1; col < wid - 1; col++) {
        result += removePixel(col, row, even);
      }
    }
    even = !even;
    image = clone(data2);
    return result;
  };

  int n = 0;
  do {
    stepCounter.fill(0);
    n = thinStep();
  } while (n > 0);

  return image;
}

// [[Rcpp::export]]
NumericVector idw_interpolation_cpp(NumericVector x, NumericVector y, NumericVector values,
                                    NumericVector new_x, NumericVector new_y, double power = 2) {
  NumericMatrix distances(new_x.size(), x.size());
  for (int i = 0; i < new_x.size(); ++i) {
    distances(i, _) = sqrt(pow(x - new_x[i], 2) + pow(y - new_y[i], 2));
  }

  NumericVector results(new_x.size(), NA_REAL);

  for (int i = 0; i < new_x.size(); ++i) {
    NumericVector weights = 1.0 / pow(distances(i, _), power);
    double weighted_sum = sum(weights * values);
    double total_weight = sum(weights);
    results[i] = total_weight > 0 ? weighted_sum / total_weight : NA_REAL;
  }
  return results;
}

NumericMatrix adjust_bbox(NumericMatrix coords, double width, double height) {
  NumericVector cent = colMeans(coords(Range(0, 3), _));
  double xmin = cent[0] - width / 2;
  double xmax = cent[0] + width / 2;
  double ymin = cent[1] - height / 2;
  double ymax = cent[1] + height / 2;

  NumericMatrix new_bbox(5, 2);
  new_bbox(0, 0) = xmin; new_bbox(0, 1) = ymin;
  new_bbox(1, 0) = xmin; new_bbox(1, 1) = ymax;
  new_bbox(2, 0) = xmax; new_bbox(2, 1) = ymax;
  new_bbox(3, 0) = xmax; new_bbox(3, 1) = ymin;
  new_bbox(4, 0) = xmin; new_bbox(4, 1) = ymin;

  return new_bbox;
}

// [[Rcpp::export]]
IntegerMatrix help_label(IntegerMatrix matrix, int max_gap = 2) {
  int rows = matrix.nrow();
  int cols = matrix.ncol();
  IntegerMatrix labels(rows, cols);
  int current_label = 0;

  auto is_within_gap = [&](int r1, int c1, int r2, int c2) {
    return abs(r1 - r2) <= max_gap && abs(c1 - c2) <= max_gap;
  };

  std::vector<int> stack_r, stack_c;

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      if (matrix(r, c) == 1 && labels(r, c) == 0) {
        ++current_label;
        stack_r.push_back(r);
        stack_c.push_back(c);

        while (!stack_r.empty()) {
          int cr = stack_r.back();
          int cc = stack_c.back();
          stack_r.pop_back();
          stack_c.pop_back();

          if (labels(cr, cc) == 0) {
            labels(cr, cc) = current_label;

            for (int dr = -max_gap; dr <= max_gap; ++dr) {
              for (int dc = -max_gap; dc <= max_gap; ++dc) {
                if (abs(dr) + abs(dc) > 0) {
                  int nr = cr + dr;
                  int nc = cc + dc;

                  if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                    if (matrix(nr, nc) == 1 && labels(nr, nc) == 0 && is_within_gap(cr, cc, nr, nc)) {
                      stack_r.push_back(nr);
                      stack_c.push_back(nc);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return labels;
}

// [[Rcpp::export]]
NumericVector rcpp_st_perimeter(List sf_coords) {
  int n = sf_coords.size();
  NumericVector perimeters(n);

  for (int i = 0; i < n; ++i) {
    List geom = sf_coords[i];
    double total_perimeter = 0.0;

    for (int j = 0; j < geom.size(); ++j) {
      NumericMatrix ring = geom[j];
      double ring_perimeter = 0.0;
      int rows = ring.nrow();

      for (int k = 0; k < rows - 1; ++k) {
        double dx = ring(k + 1, 0) - ring(k, 0);
        double dy = ring(k + 1, 1) - ring(k, 1);
        ring_perimeter += sqrt(dx * dx + dy * dy);
      }

      total_perimeter += ring_perimeter;
    }
    perimeters[i] = total_perimeter;
  }
  return perimeters;
}

static std::string generate_random_hex(int length) {
  const char hex_chars[] = "0123456789abcdef";
  std::string result(length, '0');
  GetRNGstate();
  for (int i = 0; i < length; i++) {
    result[i] = hex_chars[(int)(unif_rand() * 16)];
  }
  PutRNGstate();
  return result;
}

// [[Rcpp::export]]
std::string uuid_v7() {
  auto now = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
  long long timestamp = duration.count();

  std::stringstream ss;
  ss << std::hex << std::setw(12) << std::setfill('0') << (timestamp & 0xFFFFFFFFFFFF);
  std::string timestamp_hex = ss.str();

  std::string time_low = timestamp_hex.substr(0, 8);
  std::string time_mid = timestamp_hex.substr(8, 4);
  std::string time_high_and_version = generate_random_hex(4);
  time_high_and_version[0] = '7';

  std::string variant_and_sequence = generate_random_hex(4);
  GetRNGstate();
  variant_and_sequence[0] = "89ab"[(int)(unif_rand() * 4)];
  PutRNGstate();

  std::string node = generate_random_hex(12);
  std::string uuid = time_low + "-" + time_mid + "-" + time_high_and_version +
    "-" + variant_and_sequence + "-" + node;

  if (uuid.length() != 36) {
    throw std::runtime_error("Generated UUID has incorrect length: " + uuid);
  }

  return uuid;
}

// [[Rcpp::export]]
double helper_entropy(NumericVector values, int precision = 2) {
  std::unordered_map<double, int> freq;
  int n = values.size();

  double scale = pow(10.0, precision);
  for (int i = 0; i < n; ++i) {
    double rounded_val = round(values[i] * scale) / scale;
    freq[rounded_val]++;
  }

  double entropy = 0.0;
  for (auto& pair : freq) {
    double prob = static_cast<double>(pair.second) / n;
    entropy -= prob * log(prob);
  }

  return entropy;
}

// [[Rcpp::export]]
CharacterVector corners_to_wkt(List cornersList) {
  int nPlots = cornersList.size();
  CharacterVector out(nPlots);

  for (int k = 0; k < nPlots; ++k) {
    NumericVector v = as<NumericVector>(cornersList[k]);
    int len = v.size();
    if (len < 8 || (len % 2) != 0) {
      stop("Each element must be an even-length numeric vector of at least 8 elements");
    }
    int rows = len / 2;

    NumericMatrix m(rows, 2);
    for (int i = 0; i < rows; ++i) {
      m(i, 0) = v[i];
      m(i, 1) = v[i + rows];
    }

    NumericMatrix c4(4, 2);
    for (int i = 0; i < 4; ++i) {
      c4(i, 0) = m(i, 0);
      c4(i, 1) = m(i, 1);
    }

    double dx1 = c4(0,0) - c4(1,0);
    double dy1 = c4(0,1) - c4(1,1);
    double d1  = dx1*dx1 + dy1*dy1;
    double dx2 = c4(1,0) - c4(2,0);
    double dy2 = c4(1,1) - c4(2,1);
    double d2  = dx2*dx2 + dy2*dy2;

    double x1, y1, x2, y2;
    if (d1 < d2) {
      x1 = (c4(0,0) + c4(1,0)) * 0.5;
      y1 = (c4(0,1) + c4(1,1)) * 0.5;
      x2 = (c4(2,0) + c4(3,0)) * 0.5;
      y2 = (c4(2,1) + c4(3,1)) * 0.5;
    } else {
      x1 = (c4(1,0) + c4(2,0)) * 0.5;
      y1 = (c4(1,1) + c4(2,1)) * 0.5;
      x2 = (c4(3,0) + c4(0,0)) * 0.5;
      y2 = (c4(3,1) + c4(0,1)) * 0.5;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6)
        << "LINESTRING(" << x1 << " " << y1 << ","
        << x2 << " " << y2 << ")";

    out[k] = oss.str();
  }

  return out;
}

// [[Rcpp::export]]
arma::cube correct_image_rcpp(const arma::cube& img,
                              const arma::mat& K,
                              std::string model) {
  int n_rows = img.n_rows;
  int n_cols = img.n_cols;
  arma::cube out_img(n_rows, n_cols, 3, arma::fill::zeros);
  arma::rowvec T_pixel(3);
  double R, G, B;

  if (model == "ccm" || model == "linear") {
    if (K.n_rows != 3 || K.n_cols != 3) {
      Rcpp::stop("Erro: Para o modelo 'ccm', 'k_mat' deve ser uma matriz 3x3.");
    }
    arma::rowvec S_pixel(3);
    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        S_pixel(0) = img(i, j, 0);
        S_pixel(1) = img(i, j, 1);
        S_pixel(2) = img(i, j, 2);

        T_pixel = S_pixel * K;

        out_img(i, j, 0) = std::max(0.0, std::min(255.0, T_pixel(0)));
        out_img(i, j, 1) = std::max(0.0, std::min(255.0, T_pixel(1)));
        out_img(i, j, 2) = std::max(0.0, std::min(255.0, T_pixel(2)));
      }
    }
  } else if (model == "affine") {
    if (K.n_rows != 4 || K.n_cols != 3) {
      Rcpp::stop("Erro: Para o modelo 'affine', 'k_mat' deve ser uma matriz 4x3.");
    }
    arma::rowvec S_pixel(4);
    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        S_pixel(0) = 1.0;
        S_pixel(1) = img(i, j, 0);
        S_pixel(2) = img(i, j, 1);
        S_pixel(3) = img(i, j, 2);

        T_pixel = S_pixel * K;

        out_img(i, j, 0) = std::max(0.0, std::min(255.0, T_pixel(0)));
        out_img(i, j, 1) = std::max(0.0, std::min(255.0, T_pixel(1)));
        out_img(i, j, 2) = std::max(0.0, std::min(255.0, T_pixel(2)));
      }
    }
  } else if (model == "white_balance") {
    double g_r = 1.0, g_g = 1.0, g_b = 1.0;
    if (K.n_rows == 3 && K.n_cols == 3) {
      g_r = K(0, 0); g_g = K(1, 1); g_b = K(2, 2);
    } else if (K.n_rows == 3) {
      g_r = K(0, 0); g_g = K(1, 0); g_b = K(2, 0);
    } else if (K.n_cols == 3) {
      g_r = K(0, 0); g_g = K(0, 1); g_b = K(0, 2);
    }
    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        out_img(i, j, 0) = std::max(0.0, std::min(255.0, img(i, j, 0) * g_r));
        out_img(i, j, 1) = std::max(0.0, std::min(255.0, img(i, j, 1) * g_g));
        out_img(i, j, 2) = std::max(0.0, std::min(255.0, img(i, j, 2) * g_b));
      }
    }
  } else if (model == "cubic") {
    if (K.n_rows != 9) {
      Rcpp::stop("Erro: 'k_mat' tem %i linhas, mas o modelo 'cubic' espera 9.", K.n_rows);
    }
    arma::rowvec S_pixel(9);
    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        R = img(i, j, 0);
        G = img(i, j, 1);
        B = img(i, j, 2);

        S_pixel(0) = R;
        S_pixel(1) = G;
        S_pixel(2) = B;
        S_pixel(3) = R * R;
        S_pixel(4) = G * G;
        S_pixel(5) = B * B;
        S_pixel(6) = R * R * R;
        S_pixel(7) = G * G * G;
        S_pixel(8) = B * B * B;

        T_pixel = S_pixel * K;

        out_img(i, j, 0) = std::max(0.0, std::min(255.0, T_pixel(0)));
        out_img(i, j, 1) = std::max(0.0, std::min(255.0, T_pixel(1)));
        out_img(i, j, 2) = std::max(0.0, std::min(255.0, T_pixel(2)));
      }
    }
  } else if (model == "root_polynomial") {
    if (K.n_rows != 20) {
      Rcpp::stop("Erro: 'k_mat' tem %i linhas, mas o modelo 'root_polynomial' espera 20.", K.n_rows);
    }
    arma::rowvec S_pixel(20);
    double R2, G2, B2;

    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        R = img(i, j, 0);
        G = img(i, j, 1);
        B = img(i, j, 2);

        R2 = R * R;
        G2 = G * G;
        B2 = B * B;

        S_pixel(0) = 1.0;
        S_pixel(1) = R;
        S_pixel(2) = G;
        S_pixel(3) = B;
        S_pixel(4) = R * G;
        S_pixel(5) = R * B;
        S_pixel(6) = G * B;
        S_pixel(7) = R2;
        S_pixel(8) = G2;
        S_pixel(9) = B2;
        S_pixel(10) = R2 * R;
        S_pixel(11) = G2 * G;
        S_pixel(12) = B2 * B;
        S_pixel(13) = R2 * G;
        S_pixel(14) = R2 * B;
        S_pixel(15) = G2 * R;
        S_pixel(16) = G2 * B;
        S_pixel(17) = B2 * R;
        S_pixel(18) = B2 * G;
        S_pixel(19) = R * G * B;

        T_pixel = S_pixel * K;

        out_img(i, j, 0) = std::max(0.0, std::min(255.0, T_pixel(0)));
        out_img(i, j, 1) = std::max(0.0, std::min(255.0, T_pixel(1)));
        out_img(i, j, 2) = std::max(0.0, std::min(255.0, T_pixel(2)));
      }
    }
  } else {
    Rcpp::stop("Modelo '"+ model +"' não é reconhecido. Use 'ccm', 'affine', 'white_balance', 'cubic' ou 'root_polynomial'.");
  }

  return out_img;
}

// [[Rcpp::export]]
List transform_polygons(List geometries,
                        double shift_x,
                        double shift_y,
                        double angle_deg,
                        double scale_x,
                        double scale_y) {

  int n_polys = geometries.size();
  double angle_rad = -angle_deg * M_PI / 180.0;
  double cos_a = std::cos(angle_rad);
  double sin_a = std::sin(angle_rad);
  double min_x = 1e9, max_x = -1e9;
  double min_y = 1e9, max_y = -1e9;
  for(int i = 0; i < n_polys; i++) {
    List poly = geometries[i];
    NumericMatrix outer_ring = poly[0];
    for(int j = 0; j < outer_ring.nrow(); j++) {
      double x = outer_ring(j, 0);
      double y = outer_ring(j, 1);
      if(x < min_x) min_x = x;
      if(x > max_x) max_x = x;
      if(y < min_y) min_y = y;
      if(y > max_y) max_y = y;
    }
  }
  double center_x = (min_x + max_x) / 2.0;
  double center_y = (min_y + max_y) / 2.0;
  List out_list(n_polys);
  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");
  for(int i = 0; i < n_polys; i++) {
    List source_poly = geometries[i];
    int n_rings = source_poly.size();
    List new_poly(n_rings);
    for(int r = 0; r < n_rings; r++) {
      NumericMatrix ring = source_poly[r];
      int n_pts = ring.nrow();
      NumericMatrix new_ring(n_pts, 2);
      for(int j = 0; j < n_pts; j++) {
        double x = ring(j, 0);
        double y = ring(j, 1);
        double dx = x - center_x;
        double dy = y - center_y;
        dx *= scale_x;
        dy *= scale_y;
        double x_rot = dx * cos_a - dy * sin_a;
        double y_rot = dx * sin_a + dy * cos_a;
        new_ring(j, 0) = x_rot + center_x + shift_x;
        new_ring(j, 1) = y_rot + center_y + shift_y;
      }
      new_poly[r] = new_ring;
    }
    new_poly.attr("class") = sfg_class;
    out_list[i] = new_poly;
  }

  return out_list;
}

// [[Rcpp::export]]
IntegerMatrix cpp_as_native_raster(SEXP img_sexp, IntegerVector dims) {
  int w = dims[0];
  int h = dims[1];
  int nch = (dims.size() >= 3) ? dims[2] : 1;
  int npix = w * h;

  IntegerMatrix mat(h, w);
  int* out = INTEGER(mat);

  if (TYPEOF(img_sexp) == RAWSXP) {
    const uint8_t* ptr = RAW(img_sexp);

    if (nch >= 3) {
      const int stride_g = npix;
      const int stride_b = npix * 2;
      const bool has_alpha = (nch >= 4);
      const int stride_a = has_alpha ? npix * 3 : 0;

      for (int i = 0; i < npix; i++) {
        uint32_t rv = ptr[i];
        uint32_t gv = ptr[i + stride_g];
        uint32_t bv = ptr[i + stride_b];
        uint32_t av = has_alpha ? ptr[i + stride_a] : 255u;

        out[i] = static_cast<int>(
          (av << 24) | (bv << 16) | (gv << 8) | rv
        );
      }
    } else {
      for (int i = 0; i < npix; i++) {
        uint32_t gv = ptr[i];
        out[i] = static_cast<int>(
          (255u << 24) | (gv << 16) | (gv << 8) | gv
        );
      }
    }
  } else if (TYPEOF(img_sexp) == REALSXP) {
    const double* ptr = REAL(img_sexp);

    double max_v = 0.0;
    for (int i = 0; i < npix; i++) {
      if (std::abs(ptr[i]) > max_v) max_v = std::abs(ptr[i]);
    }
    double scale = (max_v <= 1.0) ? 255.0 : ((max_v <= 255.0) ? 1.0 : (255.0 / max_v));

    if (nch >= 3) {
      const int stride_g = npix;
      const int stride_b = npix * 2;

      for (int i = 0; i < npix; i++) {
        uint32_t rv = (uint32_t)std::max(0, std::min(255, (int)(ptr[i] * scale)));
        uint32_t gv = (uint32_t)std::max(0, std::min(255, (int)(ptr[i + stride_g] * scale)));
        uint32_t bv = (uint32_t)std::max(0, std::min(255, (int)(ptr[i + stride_b] * scale)));

        out[i] = static_cast<int>(
          (255u << 24) | (bv << 16) | (gv << 8) | rv
        );
      }
    } else {
      for (int i = 0; i < npix; i++) {
        uint32_t gv = (uint32_t)std::max(0, std::min(255, (int)(ptr[i] * scale)));
        out[i] = static_cast<int>(
          (255u << 24) | (gv << 16) | (gv << 8) | gv
        );
      }
    }
  } else if (TYPEOF(img_sexp) == LGLSXP) {
    const int* ptr = LOGICAL(img_sexp);
    const int white = static_cast<int>(0xFFFFFFFFu);
    const int black = static_cast<int>(0xFF000000u);

    for (int i = 0; i < npix; i++) {
      out[i] = ptr[i] ? white : black;
    }
  }

  mat.attr("class") = "nativeRaster";
  mat.attr("channels") = nch >= 3 ? 4 : 1;

  return mat;
}

// [[Rcpp::export]]
List image_histogram_cpp(SEXP img_sexp, int nbins = 256) {
  int nrow = Rf_nrows(img_sexp);
  int ncol = Rf_ncols(img_sexp);
  SEXP dims_sexp = Rf_getAttrib(img_sexp, R_DimSymbol);
  int nch = (TYPEOF(dims_sexp) == INTSXP && Rf_length(dims_sexp) == 3 && INTEGER(dims_sexp)[2] >= 3) ? INTEGER(dims_sexp)[2] : 1;
  int plane_size = nrow * ncol;

  List res(nch);

  if (TYPEOF(img_sexp) == RAWSXP) {
    const uint8_t* ptr = RAW(img_sexp);
    for (int ch = 0; ch < nch; ++ch) {
      const uint8_t* ch_ptr = ptr + (size_t)ch * plane_size;
      IntegerVector counts(256, 0);
      int* p_counts = INTEGER(counts);
      for (int i = 0; i < plane_size; ++i) {
        p_counts[ch_ptr[i]]++;
      }
      res[ch] = counts;
    }
  } else if (TYPEOF(img_sexp) == LGLSXP || TYPEOF(img_sexp) == INTSXP) {
    const int* ptr = INTEGER(img_sexp);
    for (int ch = 0; ch < nch; ++ch) {
      const int* ch_ptr = ptr + (size_t)ch * plane_size;
      IntegerVector counts(nbins, 0);
      int* p_counts = INTEGER(counts);
      for (int i = 0; i < plane_size; ++i) {
        int v = ch_ptr[i];
        if (v >= 0 && v < nbins) p_counts[v]++;
      }
      res[ch] = counts;
    }
  } else if (TYPEOF(img_sexp) == REALSXP) {
    const double* ptr = REAL(img_sexp);
    for (int ch = 0; ch < nch; ++ch) {
      const double* ch_ptr = ptr + (size_t)ch * plane_size;
      IntegerVector counts(nbins, 0);
      int* p_counts = INTEGER(counts);
      double scale = nbins - 1.0;
      for (int i = 0; i < plane_size; ++i) {
        double v = ch_ptr[i];
        if (std::isnan(v)) continue;
        int bin = static_cast<int>(std::round(v * scale));
        if (bin < 0) bin = 0;
        if (bin >= nbins) bin = nbins - 1;
        p_counts[bin]++;
      }
      res[ch] = counts;
    }
  }
  return res;
}

// [[Rcpp::export]]
SEXP cpp_image_transpose(SEXP img_sexp) {
  SEXP dims_sexp = Rf_getAttrib(img_sexp, R_DimSymbol);
  if (Rf_isNull(dims_sexp)) {
    Rcpp::stop("Input must be a 2D or 3D array/matrix.");
  }
  int ndim = Rf_length(dims_sexp);
  int H = INTEGER(dims_sexp)[0];
  int W = INTEGER(dims_sexp)[1];
  int C = (ndim == 3) ? INTEGER(dims_sexp)[2] : 1;

  SEXPTYPE type = TYPEOF(img_sexp);
  SEXP out_sexp = PROTECT(Rf_allocVector(type, H * W * C));

  SEXP class_attr = Rf_getAttrib(img_sexp, R_ClassSymbol);
  if (!Rf_isNull(class_attr)) Rf_setAttrib(out_sexp, R_ClassSymbol, class_attr);

  SEXP cm_attr = Rf_getAttrib(img_sexp, Rf_install("colormode"));
  if (!Rf_isNull(cm_attr)) Rf_setAttrib(out_sexp, Rf_install("colormode"), cm_attr);

  SEXP cs_attr = Rf_getAttrib(img_sexp, Rf_install("colspace"));
  if (!Rf_isNull(cs_attr)) Rf_setAttrib(out_sexp, Rf_install("colspace"), cs_attr);

  SEXP gm_attr = Rf_getAttrib(img_sexp, Rf_install("gamma"));
  if (!Rf_isNull(gm_attr)) Rf_setAttrib(out_sexp, Rf_install("gamma"), gm_attr);

  SEXP out_dims;
  if (ndim == 2) {
    out_dims = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(out_dims)[0] = W;
    INTEGER(out_dims)[1] = H;
  } else {
    out_dims = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(out_dims)[0] = W;
    INTEGER(out_dims)[1] = H;
    INTEGER(out_dims)[2] = C;
  }
  Rf_setAttrib(out_sexp, R_DimSymbol, out_dims);
  UNPROTECT(1);

  int block = 32;
  size_t HW = (size_t)H * W;
  size_t WH = (size_t)W * H;

  if (type == RAWSXP) {
    const uint8_t* src = RAW(img_sexp);
    uint8_t* dst = RAW(out_sexp);
    for (int k = 0; k < C; ++k) {
      const uint8_t* src_k = src + k * HW;
      uint8_t* dst_k = dst + k * WH;
      #pragma omp parallel for collapse(2) schedule(static) if (HW > 500000)
      for (int i0 = 0; i0 < H; i0 += block) {
        for (int j0 = 0; j0 < W; j0 += block) {
          int i_max = std::min(i0 + block, H);
          int j_max = std::min(j0 + block, W);
          for (int i = i0; i < i_max; ++i) {
            for (int j = j0; j < j_max; ++j) {
              dst_k[j + (size_t)i * W] = src_k[i + (size_t)j * H];
            }
          }
        }
      }
    }
  } else if (type == REALSXP) {
    const double* src = REAL(img_sexp);
    double* dst = REAL(out_sexp);
    for (int k = 0; k < C; ++k) {
      const double* src_k = src + k * HW;
      double* dst_k = dst + k * WH;
      #pragma omp parallel for collapse(2) schedule(static) if (HW > 500000)
      for (int i0 = 0; i0 < H; i0 += block) {
        for (int j0 = 0; j0 < W; j0 += block) {
          int i_max = std::min(i0 + block, H);
          int j_max = std::min(j0 + block, W);
          for (int i = i0; i < i_max; ++i) {
            for (int j = j0; j < j_max; ++j) {
              dst_k[j + (size_t)i * W] = src_k[i + (size_t)j * H];
            }
          }
        }
      }
    }
  } else if (type == INTSXP || type == LGLSXP) {
    const int* src = INTEGER(img_sexp);
    int* dst = INTEGER(out_sexp);
    for (int k = 0; k < C; ++k) {
      const int* src_k = src + k * HW;
      int* dst_k = dst + k * WH;
      #pragma omp parallel for collapse(2) schedule(static) if (HW > 500000)
      for (int i0 = 0; i0 < H; i0 += block) {
        for (int j0 = 0; j0 < W; j0 += block) {
          int i_max = std::min(i0 + block, H);
          int j_max = std::min(j0 + block, W);
          for (int i = i0; i < i_max; ++i) {
            for (int j = j0; j < j_max; ++j) {
              dst_k[j + (size_t)i * W] = src_k[i + (size_t)j * H];
            }
          }
        }
      }
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

// [[Rcpp::export]]
SEXP cpp_clahe(SEXP img_sexp, int nx = 8, int ny = 8, double clip_limit = 3.0, int nbins = 256) {
  SEXP dims_sexp = Rf_getAttrib(img_sexp, R_DimSymbol);
  if (Rf_isNull(dims_sexp)) {
    Rcpp::stop("Input must be a 2D or 3D array/matrix.");
  }
  int ndim = Rf_length(dims_sexp);
  int W = INTEGER(dims_sexp)[0];
  int H = INTEGER(dims_sexp)[1];
  int C = (ndim == 3) ? INTEGER(dims_sexp)[2] : 1;

  SEXPTYPE type = TYPEOF(img_sexp);
  if (type != RAWSXP && type != REALSXP) {
    Rcpp::stop("Image storage mode must be raw or double.");
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(type, (size_t)W * H * C));
  Rf_setAttrib(out_sexp, R_DimSymbol, dims_sexp);

  SEXP class_attr = Rf_getAttrib(img_sexp, R_ClassSymbol);
  if (!Rf_isNull(class_attr)) Rf_setAttrib(out_sexp, R_ClassSymbol, class_attr);

  SEXP cm_attr = Rf_getAttrib(img_sexp, Rf_install("colormode"));
  if (!Rf_isNull(cm_attr)) Rf_setAttrib(out_sexp, Rf_install("colormode"), cm_attr);

  SEXP cs_attr = Rf_getAttrib(img_sexp, Rf_install("colspace"));
  if (!Rf_isNull(cs_attr)) Rf_setAttrib(out_sexp, Rf_install("colspace"), cs_attr);

  SEXP gm_attr = Rf_getAttrib(img_sexp, Rf_install("gamma"));
  if (!Rf_isNull(gm_attr)) Rf_setAttrib(out_sexp, Rf_install("gamma"), gm_attr);

  nx = std::max(1, std::min(nx, W));
  ny = std::max(1, std::min(ny, H));
  nbins = std::max(2, std::min(nbins, 256));

  double tile_w = (double)W / nx;
  double tile_h = (double)H / ny;

  for (int ch = 0; ch < C; ++ch) {
    size_t plane_offset = (size_t)ch * W * H;

    std::vector<std::vector<double>> cdfs(nx * ny, std::vector<double>(nbins, 0.0));

    for (int ty = 0; ty < ny; ++ty) {
      for (int tx = 0; tx < nx; ++tx) {
        int x0 = (int)std::floor(tx * tile_w);
        int x1 = (int)std::floor((tx + 1) * tile_w);
        int y0 = (int)std::floor(ty * tile_h);
        int y1 = (int)std::floor((ty + 1) * tile_h);
        x1 = std::min(x1, W);
        y1 = std::min(y1, H);
        int npixels = (x1 - x0) * (y1 - y0);
        if (npixels <= 0) continue;

        std::vector<int> hist(nbins, 0);

        if (type == RAWSXP) {
          const uint8_t* ptr = RAW(img_sexp) + plane_offset;
          for (int y = y0; y < y1; ++y) {
            size_t row_off = (size_t)y * W;
            for (int x = x0; x < x1; ++x) {
              int val = ptr[x + row_off];
              int bin = (val * (nbins - 1)) / 255;
              hist[bin]++;
            }
          }
        } else {
          const double* ptr = REAL(img_sexp) + plane_offset;
          for (int y = y0; y < y1; ++y) {
            size_t row_off = (size_t)y * W;
            for (int x = x0; x < x1; ++x) {
              double val = ptr[x + row_off];
              val = std::max(0.0, std::min(1.0, val));
              int bin = (int)(val * (nbins - 1));
              bin = std::max(0, std::min(nbins - 1, bin));
              hist[bin]++;
            }
          }
        }

        int clip_val = (int)std::round(clip_limit * npixels / nbins);
        if (clip_val < 1) clip_val = 1;

        int excess = 0;
        for (int b = 0; b < nbins; ++b) {
          if (hist[b] > clip_val) {
            excess += (hist[b] - clip_val);
            hist[b] = clip_val;
          }
        }

        int add_per_bin = excess / nbins;
        int rem = excess % nbins;
        for (int b = 0; b < nbins; ++b) {
          hist[b] += add_per_bin;
        }
        for (int b = 0; b < rem; ++b) {
          hist[b]++;
        }

        int tile_idx = ty * nx + tx;
        double sum = 0.0;
        for (int b = 0; b < nbins; ++b) {
          sum += hist[b];
          cdfs[tile_idx][b] = sum / npixels;
        }
      }
    }

    if (type == RAWSXP) {
      const uint8_t* src = RAW(img_sexp) + plane_offset;
      uint8_t* dst = RAW(out_sexp) + plane_offset;

      #pragma omp parallel for collapse(2) if(W * H > 100000)
      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          double gx = (x + 0.5) / tile_w - 0.5;
          double gy = (y + 0.5) / tile_h - 0.5;

          int tx0 = (int)std::floor(gx);
          int ty0 = (int)std::floor(gy);

          tx0 = std::max(0, std::min(nx - 2, tx0));
          ty0 = std::max(0, std::min(ny - 2, ty0));

          int tx1 = tx0 + 1;
          int ty1 = ty0 + 1;

          double fx = (gx - tx0);
          double fy = (gy - ty0);
          fx = std::max(0.0, std::min(1.0, fx));
          fy = std::max(0.0, std::min(1.0, fy));

          size_t p_idx = (size_t)y * W + x;
          int val = src[p_idx];
          int bin = (val * (nbins - 1)) / 255;

          double c00 = cdfs[ty0 * nx + tx0][bin];
          double c10 = cdfs[ty0 * nx + tx1][bin];
          double c01 = cdfs[ty1 * nx + tx0][bin];
          double c11 = cdfs[ty1 * nx + tx1][bin];

          double cdf_val = (1.0 - fx) * (1.0 - fy) * c00 +
                           fx * (1.0 - fy) * c10 +
                           (1.0 - fx) * fy * c01 +
                           fx * fy * c11;

          int out_val = (int)std::round(cdf_val * 255.0);
          dst[p_idx] = (uint8_t)std::max(0, std::min(255, out_val));
        }
      }
    } else {
      const double* src = REAL(img_sexp) + plane_offset;
      double* dst = REAL(out_sexp) + plane_offset;

      #pragma omp parallel for collapse(2) if(W * H > 100000)
      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          double gx = (x + 0.5) / tile_w - 0.5;
          double gy = (y + 0.5) / tile_h - 0.5;

          int tx0 = (int)std::floor(gx);
          int ty0 = (int)std::floor(gy);

          tx0 = std::max(0, std::min(nx - 2, tx0));
          ty0 = std::max(0, std::min(ny - 2, ty0));

          int tx1 = tx0 + 1;
          int ty1 = ty0 + 1;

          double fx = (gx - tx0);
          double fy = (gy - ty0);
          fx = std::max(0.0, std::min(1.0, fx));
          fy = std::max(0.0, std::min(1.0, fy));

          size_t p_idx = (size_t)y * W + x;
          double val = src[p_idx];
          val = std::max(0.0, std::min(1.0, val));
          int bin = (int)(val * (nbins - 1));
          bin = std::max(0, std::min(nbins - 1, bin));

          double c00 = cdfs[ty0 * nx + tx0][bin];
          double c10 = cdfs[ty0 * nx + tx1][bin];
          double c01 = cdfs[ty1 * nx + tx0][bin];
          double c11 = cdfs[ty1 * nx + tx1][bin];

          double cdf_val = (1.0 - fx) * (1.0 - fy) * c00 +
                           fx * (1.0 - fy) * c10 +
                           (1.0 - fx) * fy * c01 +
                           fx * fy * c11;

          dst[p_idx] = std::max(0.0, std::min(1.0, cdf_val));
        }
      }
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

// [[Rcpp::export]]
NumericMatrix cpp_make_brush(int size, std::string shape = "disc", bool step = true, double sigma = 0.3, double angle = 45.0) {
  if (size % 2 == 0) size += 1;
  int radius = (size - 1) / 2;
  NumericMatrix mat(size, size);

  std::string sh = shape;
  for (auto &c : sh) c = std::tolower(c);
  if (sh == "circle") sh = "disc";
  if (sh == "square") sh = "box";

  if (sh == "box") {
    std::fill(mat.begin(), mat.end(), 1.0);
  } else if (sh == "disc") {
    double r2 = (double)radius * radius;
    for (int y = -radius; y <= radius; ++y) {
      for (int x = -radius; x <= radius; ++x) {
        double d2 = (double)(x * x + y * y);
        if (step) {
          mat(x + radius, y + radius) = (d2 <= r2 + 1e-7) ? 1.0 : 0.0;
        } else {
          double dist = std::sqrt(d2);
          if (dist <= radius - 0.5) {
            mat(x + radius, y + radius) = 1.0;
          } else if (dist >= radius + 0.5) {
            mat(x + radius, y + radius) = 0.0;
          } else {
            mat(x + radius, y + radius) = radius + 0.5 - dist;
          }
        }
      }
    }
  } else if (sh == "diamond") {
    for (int y = -radius; y <= radius; ++y) {
      for (int x = -radius; x <= radius; ++x) {
        int d = std::abs(x) + std::abs(y);
        mat(x + radius, y + radius) = (d <= radius) ? 1.0 : 0.0;
      }
    }
  } else if (sh == "line") {
    double rad = angle * M_PI / 180.0;
    double dx = std::cos(rad);
    double dy = std::sin(rad);
    for (int y = -radius; y <= radius; ++y) {
      for (int x = -radius; x <= radius; ++x) {
        double proj = x * dx + y * dy;
        double dist2 = (x * x + y * y) - proj * proj;
        if (dist2 <= 0.5 && std::abs(proj) <= radius + 0.5) {
          mat(x + radius, y + radius) = 1.0;
        }
      }
    }
  } else if (sh == "gaussian") {
    if (sigma <= 0) sigma = (double)radius / 3.0;
    double sig2 = 2.0 * sigma * sigma;
    double max_v = 0.0;
    for (int y = -radius; y <= radius; ++y) {
      for (int x = -radius; x <= radius; ++x) {
        double v = std::exp(-(double)(x * x + y * y) / sig2);
        mat(x + radius, y + radius) = v;
        if (v > max_v) max_v = v;
      }
    }
    if (max_v > 0) {
      for (int i = 0; i < size * size; ++i) {
        mat[i] /= max_v;
      }
    }
  } else {
    double r2 = (double)radius * radius;
    for (int y = -radius; y <= radius; ++y) {
      for (int x = -radius; x <= radius; ++x) {
        mat(x + radius, y + radius) = ((double)(x * x + y * y) <= r2 + 1e-7) ? 1.0 : 0.0;
      }
    }
  }
  return mat;
}

// [[Rcpp::export]]
RawVector float_to_raw_cpp(NumericVector input) {
  std::size_t n = input.size();
  RawVector out(n);
  const double* p_in = input.begin();
  unsigned char* p_out = reinterpret_cast<unsigned char*>(out.begin());

  #pragma omp parallel for schedule(static) if(n > 50000)
  for (std::size_t i = 0; i < n; ++i) {
    double val = p_in[i];
    if (val <= 0.0) {
      p_out[i] = 0;
    } else if (val >= 1.0) {
      p_out[i] = 255;
    } else {
      p_out[i] = static_cast<unsigned char>(val * 255.0 + 0.5);
    }
  }

  if (input.hasAttribute("dim")) {
    out.attr("dim") = input.attr("dim");
  }

  return out;
}

// [[Rcpp::export]]
NumericVector raw_to_float_cpp(RawVector input) {
  std::size_t n = input.size();
  NumericVector out(n);
  const unsigned char* p_in = reinterpret_cast<const unsigned char*>(input.begin());
  double* p_out = out.begin();

  #pragma omp parallel for schedule(static) if(n > 50000)
  for (std::size_t i = 0; i < n; ++i) {
    p_out[i] = static_cast<double>(p_in[i]) / 255.0;
  }

  if (input.hasAttribute("dim")) {
    out.attr("dim") = input.attr("dim");
  }

  return out;
}
