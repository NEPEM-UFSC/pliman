// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::RObject help_area(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::NumericMatrix coord_matrix(coord);
    double area = 0;
    int n = coord_matrix.nrow();
    for (int i = 0; i < n; ++i) {
      int j = i + 1;
      if (j == n) {
        j = 0;
      }
      area += coord_matrix(i, 0) * coord_matrix(j, 1);
      area -= coord_matrix(j, 0) * coord_matrix(i, 1);
    }
    return Rcpp::wrap(std::abs(area / 2));
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::List coord_list(coord);
    std::vector<double> result;
    for (int i = 0; i < coord_list.length(); i++) {
      Rcpp::NumericMatrix coord_matrix = coord_list[i];
      double area = 0;
      int n = coord_matrix.nrow();
      for (int i = 0; i < n; ++i) {
        int j = i + 1;
        if (j == n) {
          j = 0;
        }
        area += coord_matrix(i, 0) * coord_matrix(j, 1);
        area -= coord_matrix(j, 0) * coord_matrix(i, 1);
      }
      result.push_back(std::abs(area / 2));
    }
    return Rcpp::wrap(result);
  } else {
    Rcpp::stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}




// [[Rcpp::export]]
NumericMatrix help_slide(NumericMatrix coord, int fp = 1) {
  int n = coord.nrow();
  NumericMatrix result(n, coord.ncol());
  for(int i = 0; i < n; i++) {
    int j = (i + fp - 1) % n;
    result(i, _) = coord(j, _);
  }
  return result;
}




// [[Rcpp::export]]
NumericVector help_distpts(NumericMatrix data) {
  int n = data.nrow();
  NumericVector distances(n - 1);
  double dx, dy;
  for (int i = 0; i < n - 1; i++) {
    dx = data(i + 1, 0) - data(i, 0);
    dy = data(i + 1, 1) - data(i, 1);
    distances[i] = sqrt(dx * dx + dy * dy);
  }
  return distances;
}





// [[Rcpp::export]]
NumericVector help_centdist(NumericMatrix data) {
  int n = data.nrow();
  int m = data.ncol();
  NumericVector centroid(m);
  NumericVector distances(n);
  // calculate centroid
  for (int j = 0; j < m; j++) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += data(i, j);
    }
    centroid[j] = sum / n;
  }
  // calculate Euclidean distances
  for (int i = 0; i < n; i++) {
    double dist = 0;
    for (int j = 0; j < m; j++) {
      dist += pow(data(i, j) - centroid[j], 2);
    }
    distances(i) = sqrt(dist);
  }
  return distances;
}


arma::vec centmass(const arma::mat& coords) {
  arma::vec center(2);
  center.fill(0);
  double area = 0;
  int n = coords.n_rows;
  for (int i = 0; i < n; i++) {
    double x1 = coords(i, 0);
    double y1 = coords(i, 1);
    double x2 = coords((i + 1) % n, 0);
    double y2 = coords((i + 1) % n, 1);
    double a = x1 * y2 - x2 * y1;
    area += a;
    center(0) += (x1 + x2) * a;
    center(1) += (y1 + y2) * a;
  }
  area /= 2;
  center /= 6 * area;
  return center;
}

// [[Rcpp::export]]
NumericVector help_centdist2(NumericMatrix data) {
  int n = data.nrow();
  NumericVector distances(n);

  // Convert NumericMatrix to arma::mat
  arma::mat coords(data.begin(), data.nrow(), data.ncol(), false);

  // Calculate the center of mass using help_mc
  arma::vec center_of_mass = centmass(coords);

  // Calculate Euclidean distances from each point to the center of mass
  for (int i = 0; i < n; i++) {
    double dist = 0;
    for (int j = 0; j < data.ncol(); j++) {
      dist += pow(data(i, j) - center_of_mass[j], 2);
    }
    distances[i] = sqrt(dist);
  }

  return distances;
}



// [[Rcpp::export]]
NumericMatrix help_rotate(NumericMatrix polygon, double angle) {
  int n = polygon.nrow();
  NumericMatrix rotatedPolygon(n, 2);
  double theta = angle * datum::pi / 180;
  arma::mat R = { {cos(theta), -sin(theta)},
  {sin(theta),  cos(theta)} };

  arma::mat P = arma::mat(polygon.begin(), n, 2, false);
  arma::mat RP = P * R;
  rotatedPolygon = wrap(RP);

  return rotatedPolygon;
}



// [[Rcpp::export]]
NumericMatrix help_align(NumericMatrix coord) {
  arma::mat coord_arma(coord.begin(), coord.nrow(), coord.ncol(), false);
  arma::mat var_ = cov(coord_arma);
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, var_);
  arma::mat aligned = coord_arma * eigvec;
  return wrap(flipud(aligned));
}


// [[Rcpp::export]]
arma::mat help_lw(SEXP coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
    arma::mat var_ = cov(coord_arma);
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, var_);
    arma::mat coordinates = flipud(coord_arma * eigvec);
    double y_min = arma::min(coordinates.col(0));
    double y_max = arma::max(coordinates.col(0));
    double x_min = arma::min(coordinates.col(1));
    double x_max = arma::max(coordinates.col(1));
    arma::vec lw = {x_max - x_min, y_max - y_min};
    return lw.t();
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list(coord);
    int n = coord_list.size();
    arma::mat lw_mat(n, 2);
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = Rcpp::as<NumericMatrix>(coord_list[i]);
      arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
      arma::mat var_ = cov(coord_arma);
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, var_);
      arma::mat coordinates = flipud(coord_arma * eigvec);
      double y_min = arma::min(coordinates.col(0));
      double y_max = arma::max(coordinates.col(0));
      double x_min = arma::min(coordinates.col(1));
      double x_max = arma::max(coordinates.col(1));
      arma::vec lw = {x_max - x_min, y_max - y_min};
      lw_mat.row(i) = lw.t();
    }
    return lw_mat;
  } else {
    Rcpp::stop("Input must be either a matrix or a list of matrices");
  }
}




// [[Rcpp::export]]
Rcpp::RObject help_eigen_ratio(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    arma::mat mat_coord = Rcpp::as<arma::mat>(coord);
    arma::mat covmat = cov(mat_coord);
    arma::vec eigenval;
    arma::mat eigenvec;
    eig_sym(eigenval, eigenvec, covmat);
    return Rcpp::wrap(eigenval[0]/eigenval[1]);
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::List list_coord = Rcpp::as<Rcpp::List>(coord);
    int n = list_coord.size();
    arma::vec result(n);
    for (int i = 0; i < n; i++) {
      arma::mat mat_coord = Rcpp::as<arma::mat>(list_coord[i]);
      arma::mat covmat = cov(mat_coord);
      arma::vec eigenval;
      arma::mat eigenvec;
      eig_sym(eigenval, eigenvec, covmat);
      result[i] = eigenval[0]/eigenval[1];
    }
    return Rcpp::wrap(result);
  } else {
    stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}








// Function calculates distance
// between two points


long dist(pair<long, long> p1,
          pair<long, long> p2)
{
  long x0 = p1.first - p2.first;
  long y0 = p1.second - p2.second;
  return sqrt(x0 * x0 + y0 * y0);
}


// [[Rcpp::export]]
Rcpp::RObject help_calliper(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    int n = coord_matrix.nrow();
    double Max = 0;
    for(int i = 0; i < n; i++)
    {
      for(int j = i + 1; j < n; j++)
      {
        pair<long, long> point1 = make_pair(coord_matrix(i, 0), coord_matrix(i, 1));
        pair<long, long> point2 = make_pair(coord_matrix(j, 0), coord_matrix(j, 1));
        Max = max(Max, (double)dist(point1, point2));
      }
    }
    return Rcpp::wrap(sqrt(Max));
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list = as<List>(coord);
    std::vector<double> result;
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = coord_list[i];
      int n = coord_matrix.nrow();
      double Max = 0;
      for(int i = 0; i < n; i++)
      {
        for(int j = i + 1; j < n; j++)
        {
          pair<long, long> point1 = make_pair(coord_matrix(i, 0), coord_matrix(i, 1));
          pair<long, long> point2 = make_pair(coord_matrix(j, 0), coord_matrix(j, 1));
          Max = max(Max, (double)dist(point1, point2));
        }
      }
      result.push_back(sqrt(Max));
    }
    return Rcpp::wrap(result);
  } else {
    stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}





// [[Rcpp::export]]
Rcpp::RObject help_elongation(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
    arma::mat var_ = cov(coord_arma);
    arma::vec eigval;
    arma::vec lw;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, var_);
    arma::mat coordinates = fliplr(coord_arma * eigvec);
    double x_min = arma::min(coordinates.col(0));
    double x_max = arma::max(coordinates.col(0));
    double y_min = arma::min(coordinates.col(1));
    double y_max = arma::max(coordinates.col(1));
    lw = 1 - (y_max - y_min) / (x_max - x_min);
    return Rcpp::wrap(lw);
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list(coord);
    int n = coord_list.size();
    arma::vec lw_vec(n);
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = Rcpp::as<NumericMatrix>(coord_list[i]);
      arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
      arma::mat var_ = cov(coord_arma);
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, var_);
      arma::mat coordinates = fliplr(coord_arma * eigvec);
      double x_min = arma::min(coordinates.col(0));
      double x_max = arma::max(coordinates.col(0));
      double y_min = arma::min(coordinates.col(1));
      double y_max = arma::max(coordinates.col(1));
      lw_vec(i) = 1 - (y_max - y_min) / (x_max - x_min);
    }
    return Rcpp::wrap(lw_vec);
  } else {
    Rcpp::stop("Input must be either a matrix or a list of matrices");
  }
}


// [[Rcpp::export]]
arma::mat help_flip_y(arma::mat shape) {
  shape.col(1) = -shape.col(1);
  return shape;
}


// [[Rcpp::export]]
arma::mat help_flip_x(arma::mat shape) {
  shape.col(0) = -shape.col(0);
  return shape;
}


// [[Rcpp::export]]
arma::vec help_mc(const arma::mat& coords) {
  arma::vec center(2);
  center.fill(0);
  double area = 0;
  int n = coords.n_rows;
  for (int i = 0; i < n; i++) {
    double x1 = coords(i, 0);
    double y1 = coords(i, 1);
    double x2 = coords((i + 1) % n, 0);
    double y2 = coords((i + 1) % n, 1);
    double a = x1 * y2 - x2 * y1;
    area += a;
    center(0) += (x1 + x2) * a;
    center(1) += (y1 + y2) * a;
  }
  area /= 2;
  center /= 6 * area;
  return center;
}


// [[Rcpp::export]]
NumericVector help_limits(NumericMatrix mat) {
  int nrow = mat.nrow();
  double minx = mat(0, 0), miny = mat(0, 1), maxx = mat(0, 0), maxy = mat(0, 1);
  for (int i = 0; i < nrow; i++) {
    if (mat(i, 0) < minx) minx = mat(i, 0);
    if (mat(i, 0) > maxx) maxx = mat(i, 0);
    if (mat(i, 1) < miny) miny = mat(i, 1);
    if (mat(i, 1) > maxy) maxy = mat(i, 1);
  }
  NumericVector result(4);
  result[0] = minx;
  result[1] = maxx;
  result[2] = miny;
  result[3] = maxy;
  return result;
}



// [[Rcpp::export]]
NumericVector help_moments(const NumericMatrix & data) {
  int n = data.nrow();
  double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;

  for (int i = 0; i < n; i++) {
    double x = data(i, 0);
    double y = data(i, 1);
    sum_x += x;
    sum_y += y;
    sum_x2 += x * x;
    sum_y2 += y * y;
    sum_xy += x * y;
  }

  double cov_xy = sum_xy / n - sum_x * sum_y / n / n;
  double var_x = sum_x2 / n - sum_x * sum_x / n / n;
  double var_y = sum_y2 / n - sum_y * sum_y / n / n;

  double theta = 0.5 * atan2(2 * cov_xy, var_x - var_y);
  double a = sqrt(0.5 * (var_x + var_y + sqrt(pow(var_x - var_y, 2) + 4 * pow(cov_xy, 2))));
  double b = sqrt(0.5 * (var_x + var_y - sqrt(pow(var_x - var_y, 2) + 4 * pow(cov_xy, 2))));

  NumericVector result(4);
  result[0] = fmax(a, b);
  result[1] = fmin(a, b);
  result[2] = sqrt(1 - pow(result[1]/result[0], 2));
  result[3] = theta;
  return result;
}


// Rcpp function to count the number of pixels in each connected object
// [[Rcpp::export]]
NumericVector get_area_mask(SEXP mask_sexp) {
  size_t n = Rf_xlength(mask_sexp);
  int type = TYPEOF(mask_sexp);

  if (n == 0) {
    return NumericVector(0);
  }

  // Fast-path 1: Logical matrix (LGLSXP - 0/1 booleans)
  if (type == LGLSXP) {
    const int* p = LOGICAL(mask_sexp);
    size_t count = 0;
    #if defined(_OPENMP)
    #pragma omp parallel for reduction(+:count) schedule(static)
    #endif
    for (size_t i = 0; i < n; i++) {
      if (p[i] != 0 && p[i] != NA_LOGICAL) count++;
    }
    NumericVector area(1);
    area[0] = static_cast<double>(count);
    return area;
  }

  // Fast-path 2: Integer matrix (INTSXP - label matrix from bwlabel_cpp / watershed_cpp)
  if (type == INTSXP) {
    const int* p = INTEGER(mask_sexp);

    #if defined(_OPENMP)
    int nthreads = omp_get_max_threads();
    std::vector<std::vector<uint32_t>> thread_hists(nthreads);
    std::vector<int> thread_max(nthreads, 0);

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::vector<uint32_t>& local_h = thread_hists[tid];
      local_h.assign(4096, 0);
      uint32_t* ptr = local_h.data();
      int local_max = 0;

      #pragma omp for schedule(static)
      for (size_t i = 0; i < n; i++) {
        int v = p[i];
        if (v > 0) {
          if (v > local_max) local_max = v;
          if (v < 4096) {
            ptr[v]++;
          } else {
            size_t uv = static_cast<size_t>(v);
            if (uv >= local_h.size()) local_h.resize(uv + 512, 0);
            local_h[uv]++;
            ptr = local_h.data();
          }
        }
      }
      thread_max[tid] = local_max;
    }

    int global_max = 0;
    for (int t = 0; t < nthreads; t++) {
      if (thread_max[t] > global_max) {
        global_max = thread_max[t];
      }
    }

    if (global_max <= 0) return NumericVector(0);

    int max_id = global_max;
    NumericVector area(max_id, 0.0);
    double* p_area = REAL(area);

    for (int t = 0; t < nthreads; t++) {
      const auto& h = thread_hists[t];
      size_t h_sz = std::min(h.size(), static_cast<size_t>(max_id + 1));
      for (size_t l = 1; l < h_sz; l++) {
        if (h[l] > 0) {
          p_area[l - 1] += h[l];
        }
      }
    }
    return area;

    #else
    int max_val = 0;
    for (size_t i = 0; i < n; i++) {
      if (p[i] > max_val) max_val = p[i];
    }
    if (max_val <= 0) return NumericVector(0);
    NumericVector area(max_val, 0.0);
    double* p_area = REAL(area);
    for (size_t i = 0; i < n; i++) {
      int val = p[i];
      if (val > 0 && val <= max_val) p_area[val - 1] += 1.0;
    }
    return area;
    #endif
  }

  // Fast-path 3: Double matrix (REALSXP - float labels or 0.0/1.0)
  if (type == REALSXP) {
    const double* p = REAL(mask_sexp);

    #if defined(_OPENMP)
    int nthreads = omp_get_max_threads();
    std::vector<std::vector<uint32_t>> thread_hists(nthreads);
    std::vector<int> thread_max(nthreads, 0);

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::vector<uint32_t>& local_h = thread_hists[tid];
      local_h.assign(4096, 0);
      uint32_t* ptr = local_h.data();
      int local_max = 0;

      #pragma omp for schedule(static)
      for (size_t i = 0; i < n; i++) {
        int v = static_cast<int>(p[i]);
        if (v > 0) {
          if (v > local_max) local_max = v;
          if (v < 4096) {
            ptr[v]++;
          } else {
            size_t uv = static_cast<size_t>(v);
            if (uv >= local_h.size()) local_h.resize(uv + 512, 0);
            local_h[uv]++;
            ptr = local_h.data();
          }
        }
      }
      thread_max[tid] = local_max;
    }

    int global_max = 0;
    for (int t = 0; t < nthreads; t++) {
      if (thread_max[t] > global_max) {
        global_max = thread_max[t];
      }
    }

    if (global_max <= 0) return NumericVector(0);

    int max_id = global_max;
    NumericVector area(max_id, 0.0);
    double* p_area = REAL(area);

    for (int t = 0; t < nthreads; t++) {
      const auto& h = thread_hists[t];
      size_t h_sz = std::min(h.size(), static_cast<size_t>(max_id + 1));
      for (size_t l = 1; l < h_sz; l++) {
        if (h[l] > 0) {
          p_area[l - 1] += h[l];
        }
      }
    }
    return area;

    #else
    int max_val = 0;
    for (size_t i = 0; i < n; i++) {
      int val = static_cast<int>(p[i]);
      if (val > max_val) max_val = val;
    }
    if (max_val <= 0) return NumericVector(0);
    NumericVector area(max_val, 0.0);
    double* p_area = REAL(area);
    for (size_t i = 0; i < n; i++) {
      int val = static_cast<int>(p[i]);
      if (val > 0 && val <= max_val) p_area[val - 1] += 1.0;
    }
    return area;
    #endif
  }

  stop("Input must be a numeric, integer, or logical matrix/array.");
}


// helper for convex hull using Andrew's Monotone Chain
arma::mat get_convex_hull(const arma::mat& P) {
    int n = P.n_rows, k = 0;
    if (n <= 3) return P;
    arma::mat H(2*n, 2);

    std::vector<std::pair<double, double>> pts(n);
    for(int i=0; i<n; i++) pts[i] = {P(i, 0), P(i, 1)};

    std::sort(pts.begin(), pts.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first || (a.first == b.first && a.second < b.second);
    });

    auto cross = [](const std::pair<double, double>& O, const std::pair<double, double>& A, const std::pair<double, double>& B) {
        return (A.first - O.first) * (B.second - O.second) - (A.second - O.second) * (B.first - O.first);
    };

    std::vector<std::pair<double, double>> hull(2*n);

    // Lower hull
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0) k--;
        hull[k++] = pts[i];
    }

    // Upper hull
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0) k--;
        hull[k++] = pts[i];
    }

    arma::mat res(k - 1, 2);
    for (int i = 0; i < k - 1; i++) {
        res(i, 0) = hull[i].first;
        res(i, 1) = hull[i].second;
    }
    return res;
}

// Rcpp function to find the 4 outer corners of a card/rectangle via Maximum Area Quadrilateral on Convex Hull
// [[Rcpp::export]]
NumericMatrix find_card_corners_cpp(NumericMatrix contour) {
  int n_pts = contour.nrow();
  if (n_pts < 4) {
    stop("Contour must have at least 4 points.");
  }

  // Convert contour to Armadillo matrix to compute convex hull
  arma::mat P(n_pts, 2);
  for (int i = 0; i < n_pts; i++) {
    P(i, 0) = contour(i, 0);
    P(i, 1) = contour(i, 1);
  }

  arma::mat H = get_convex_hull(P);
  int n_hull = H.n_rows;

  if (n_hull < 4) {
    NumericMatrix res(4, 2);
    for (int i = 0; i < std::min(4, n_hull); i++) {
      res(i, 0) = H(i, 0);
      res(i, 1) = H(i, 1);
    }
    return res;
  }

  // Find 4 points in convex hull that maximize the quadrilateral area
  double max_area = -1.0;
  int best_i = 0, best_j = 1, best_k = 2, best_l = 3;

  for (int i = 0; i < n_hull - 3; i++) {
    for (int j = i + 1; j < n_hull - 2; j++) {
      for (int k = j + 1; k < n_hull - 1; k++) {
        for (int l = k + 1; l < n_hull; l++) {
          double x1 = H(i, 0), y1 = H(i, 1);
          double x2 = H(j, 0), y2 = H(j, 1);
          double x3 = H(k, 0), y3 = H(k, 1);
          double x4 = H(l, 0), y4 = H(l, 1);

          double area = 0.5 * std::abs(
            (x1 * y2 - x2 * y1) +
            (x2 * y3 - x3 * y2) +
            (x3 * y4 - x4 * y3) +
            (x4 * y1 - x1 * y4)
          );

          if (area > max_area) {
            max_area = area;
            best_i = i; best_j = j; best_k = k; best_l = l;
          }
        }
      }
    }
  }

  double pts[4][2] = {
    {H(best_i, 0), H(best_i, 1)},
    {H(best_j, 0), H(best_j, 1)},
    {H(best_k, 0), H(best_k, 1)},
    {H(best_l, 0), H(best_l, 1)}
  };

  // Classify corners: BL (min x+y), TR (max x+y), BR (max x-y), TL (min x-y)
  double min_sum = 1e18, max_sum = -1e18;
  double min_diff = 1e18, max_diff = -1e18;
  int idx_bl = 0, idx_tr = 0, idx_br = 0, idx_tl = 0;

  for (int idx = 0; idx < 4; idx++) {
    double s = pts[idx][0] + pts[idx][1];
    double d = pts[idx][0] - pts[idx][1];
    if (s < min_sum) { min_sum = s; idx_bl = idx; }
    if (s > max_sum) { max_sum = s; idx_tr = idx; }
    if (d > max_diff) { max_diff = d; idx_br = idx; }
    if (d < min_diff) { min_diff = d; idx_tl = idx; }
  }

  NumericMatrix corners(4, 2);
  rownames(corners) = CharacterVector::create("bottom_left", "top_right", "bottom_right", "top_left");
  colnames(corners) = CharacterVector::create("x", "y");

  corners(0, 0) = pts[idx_bl][0]; corners(0, 1) = pts[idx_bl][1];
  corners(1, 0) = pts[idx_tr][0]; corners(1, 1) = pts[idx_tr][1];
  corners(2, 0) = pts[idx_br][0]; corners(2, 1) = pts[idx_br][1];
  corners(3, 0) = pts[idx_tl][0]; corners(3, 1) = pts[idx_tl][1];

  return corners;
}

// Forward declaration
NumericMatrix help_smoth(NumericMatrix coords, int niter);

// [[Rcpp::export]]
DataFrame poly_measures_cpp(List contours, bool calc_pcv = false) {
    int n_contours = contours.size();

    std::vector<double> mass_x(n_contours, NA_REAL);
    std::vector<double> mass_y(n_contours, NA_REAL);
    std::vector<double> area(n_contours, NA_REAL);
    std::vector<double> area_ch(n_contours, NA_REAL);
    std::vector<double> perimeter(n_contours, NA_REAL);
    std::vector<double> radius_mean(n_contours, NA_REAL);
    std::vector<double> radius_min(n_contours, NA_REAL);
    std::vector<double> radius_max(n_contours, NA_REAL);
    std::vector<double> radius_sd(n_contours, NA_REAL);
    std::vector<double> radius_ratio(n_contours, NA_REAL);
    std::vector<double> diam_mean(n_contours, NA_REAL);
    std::vector<double> diam_min(n_contours, NA_REAL);
    std::vector<double> diam_max(n_contours, NA_REAL);
    std::vector<double> caliper(n_contours, NA_REAL);
    std::vector<double> length_m(n_contours, NA_REAL);
    std::vector<double> width_m(n_contours, NA_REAL);
    std::vector<double> solidity(n_contours, NA_REAL);
    std::vector<double> convexity(n_contours, NA_REAL);
    std::vector<double> elongation(n_contours, NA_REAL);
    std::vector<double> circularity(n_contours, NA_REAL);
    std::vector<double> circularity_haralick(n_contours, NA_REAL);
    std::vector<double> circularity_norm(n_contours, NA_REAL);
    std::vector<double> eccentricity(n_contours, NA_REAL);
    std::vector<double> maj_axis(n_contours, NA_REAL);
    std::vector<double> min_axis(n_contours, NA_REAL);
    std::vector<double> theta(n_contours, NA_REAL);
    std::vector<double> pcv(n_contours, NA_REAL);

    struct ContourData {
        std::vector<double> x;
        std::vector<double> y;
        bool valid = false;
    };
    std::vector<ContourData> cdata(n_contours);
    for (int i = 0; i < n_contours; i++) {
        SEXP curr = contours[i];
        if (Rf_isNull(curr)) continue;
        NumericMatrix mat = as<NumericMatrix>(curr);
        int n_pts = mat.nrow();
        if (n_pts < 3) continue;
        cdata[i].x.resize(n_pts);
        cdata[i].y.resize(n_pts);
        const double* ptr = REAL(mat);
        for (int j = 0; j < n_pts; j++) {
            cdata[i].x[j] = ptr[j];
            cdata[i].y[j] = ptr[j + n_pts];
        }
        cdata[i].valid = true;
    }

    #pragma omp parallel for schedule(dynamic) if(n_contours > 10)
    for (int i = 0; i < n_contours; i++) {
        if (!cdata[i].valid) continue;
        const auto& px = cdata[i].x;
        const auto& py = cdata[i].y;
        int n_pts = px.size();

        // 1. Centroid / Mass via Shoelace area weighting
        double a = 0.0;
        double cm_x = 0.0, cm_y = 0.0;
        for (int j = 0; j < n_pts; ++j) {
            int next_j = (j + 1) % n_pts;
            double x1 = px[j], y1 = py[j];
            double x2 = px[next_j], y2 = py[next_j];
            double cross = x1 * y2 - x2 * y1;
            a += cross;
            cm_x += (x1 + x2) * cross;
            cm_y += (y1 + y2) * cross;
        }
        double abs_area = std::abs(a / 2.0);
        area[i] = abs_area;
        if (abs_area > 0.0 && a != 0.0) {
            mass_x[i] = cm_x / (6.0 * (a / 2.0));
            mass_y[i] = cm_y / (6.0 * (a / 2.0));
        } else {
            mass_x[i] = 0.0; mass_y[i] = 0.0;
        }

        // 2. Convex Hull (Monotone Chain)
        std::vector<std::pair<double, double>> pts(n_pts);
        for (int j = 0; j < n_pts; j++) pts[j] = {px[j], py[j]};
        std::sort(pts.begin(), pts.end(), [](const std::pair<double, double>& u, const std::pair<double, double>& v) {
            return u.first < v.first || (u.first == v.first && u.second < v.second);
        });

        std::vector<std::pair<double, double>> hull(2 * n_pts);
        int k = 0;
        for (int j = 0; j < n_pts; ++j) {
            while (k >= 2) {
                double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                               (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
                if (cross <= 0) k--; else break;
            }
            hull[k++] = pts[j];
        }
        for (int j = n_pts - 2, t = k + 1; j >= 0; j--) {
            while (k >= t) {
                double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                               (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
                if (cross <= 0) k--; else break;
            }
            hull[k++] = pts[j];
        }

        int n_ch = k - 1;
        double a_ch = 0.0;
        double p_ch = 0.0;
        double max_d_sq = 0.0;

        for (int j = 0; j < n_ch; ++j) {
            int next_j = (j + 1) % n_ch;
            a_ch += hull[j].first * hull[next_j].second - hull[next_j].first * hull[j].second;
            double dx = hull[next_j].first - hull[j].first;
            double dy = hull[next_j].second - hull[j].second;
            p_ch += std::sqrt(dx * dx + dy * dy);

            for (int m = j + 1; m < n_ch; m++) {
                double dxx = hull[j].first - hull[m].first;
                double dyy = hull[j].second - hull[m].second;
                double d_sq = dxx * dxx + dyy * dyy;
                if (d_sq > max_d_sq) max_d_sq = d_sq;
            }
        }
        area_ch[i] = std::abs(a_ch / 2.0);
        caliper[i] = std::sqrt(max_d_sq);

        // 3. Perimeter
        double p = 0.0;
        for (int j = 0; j < n_pts - 1; j++) {
            double dx = px[j+1] - px[j];
            double dy = py[j+1] - py[j];
            p += std::sqrt(dx * dx + dy * dy);
        }
        perimeter[i] = p;

        // 4. Centdist (radius metrics)
        double cent_x = 0.0, cent_y = 0.0;
        for (int j = 0; j < n_pts; j++) {
            cent_x += px[j];
            cent_y += py[j];
        }
        cent_x /= n_pts;
        cent_y /= n_pts;

        double c_mean = 0.0, c_min = 1e15, c_max = -1e15;
        std::vector<double> cdists(n_pts);
        for (int j = 0; j < n_pts; j++) {
            double dx = px[j] - cent_x;
            double dy = py[j] - cent_y;
            double d = std::sqrt(dx * dx + dy * dy);
            cdists[j] = d;
            c_mean += d;
            if (d < c_min) c_min = d;
            if (d > c_max) c_max = d;
        }
        c_mean /= n_pts;

        double c_var = 0.0;
        for (int j = 0; j < n_pts; j++) {
            c_var += (cdists[j] - c_mean) * (cdists[j] - c_mean);
        }
        double c_sd = std::sqrt(c_var / (n_pts > 1 ? (n_pts - 1) : 1));

        radius_mean[i] = c_mean;
        radius_min[i] = c_min;
        radius_max[i] = c_max;
        radius_sd[i] = c_sd;
        radius_ratio[i] = (c_min > 0) ? (c_max / c_min) : 0.0;
        diam_mean[i] = c_mean * 2.0;
        diam_min[i] = c_min * 2.0;
        diam_max[i] = c_max * 2.0;

        // 5. Moments, theta, maj_axis, min_axis, eccentricity, length & width
        double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
        for (int j = 0; j < n_pts; j++) {
            double x = px[j], y = py[j];
            sum_x += x; sum_y += y;
            sum_x2 += x * x; sum_y2 += y * y;
            sum_xy += x * y;
        }
        double cov_xy = sum_xy / n_pts - sum_x * sum_y / n_pts / n_pts;
        double var_x = sum_x2 / n_pts - sum_x * sum_x / n_pts / n_pts;
        double var_y = sum_y2 / n_pts - sum_y * sum_y / n_pts / n_pts;
        double t = 0.5 * std::atan2(2.0 * cov_xy, var_x - var_y);

        double cos_t = std::cos(t);
        double sin_t = std::sin(t);
        double min_u = 1e15, max_u = -1e15;
        double min_v = 1e15, max_v = -1e15;
        for (int j = 0; j < n_pts; j++) {
            double u = px[j] * cos_t + py[j] * sin_t;
            double v = -px[j] * sin_t + py[j] * cos_t;
            if (u < min_u) min_u = u;
            if (u > max_u) max_u = u;
            if (v < min_v) min_v = v;
            if (v > max_v) max_v = v;
        }
        double dim1 = max_u - min_u;
        double dim2 = max_v - min_v;
        double l = std::max(dim1, dim2);
        double w = std::min(dim1, dim2);

        length_m[i] = l;
        width_m[i] = w;

        double a_maj = std::sqrt(std::max(0.0, 0.5 * (var_x + var_y + std::sqrt(std::pow(var_x - var_y, 2) + 4.0 * std::pow(cov_xy, 2)))));
        double b_min = std::sqrt(std::max(0.0, 0.5 * (var_x + var_y - std::sqrt(std::pow(var_x - var_y, 2) + 4.0 * std::pow(cov_xy, 2)))));

        maj_axis[i] = std::fmax(a_maj, b_min);
        min_axis[i] = std::fmin(a_maj, b_min);
        eccentricity[i] = (maj_axis[i] > 0) ? std::sqrt(std::max(0.0, 1.0 - std::pow(min_axis[i] / maj_axis[i], 2))) : 0.0;
        theta[i] = t;

        // 6. Shape factors
        solidity[i] = (area_ch[i] > 0) ? (abs_area / area_ch[i]) : 0.0;
        convexity[i] = (p > 0) ? (p_ch / p) : 0.0;
        elongation[i] = (l > 0) ? (1.0 - (w / l)) : 0.0;
        circularity[i] = (abs_area > 0) ? ((p * p) / abs_area) : 0.0;
        circularity_haralick[i] = (c_sd > 0) ? (c_mean / c_sd) : 0.0;
        circularity_norm[i] = (p > 0) ? ((abs_area * 4.0 * M_PI) / (p * p)) : 0.0;

        // 7. PCV (Only if requested)
        if (calc_pcv && p > 0) {
            std::vector<double> sx = px, sy = py;
            std::vector<double> smx(n_pts), smy(n_pts);
            for (int a = 0; a < 100; a++) {
                for (int k = 0; k < n_pts; k++) {
                    int prev = (k == 0) ? (n_pts - 1) : (k - 1);
                    int next = (k == n_pts - 1) ? 0 : (k + 1);
                    smx[k] = (sx[k] + sx[prev] + sx[next]) / 3.0;
                    smy[k] = (sy[k] + sy[prev] + sy[next]) / 3.0;
                }
                sx = smx; sy = smy;
            }
            double sum_d = 0.0;
            std::vector<double> sdists(n_pts);
            for (int k = 0; k < n_pts; k++) {
                double dx = px[k] - sx[k];
                double dy = py[k] - sy[k];
                double d = std::sqrt(dx * dx + dy * dy);
                sdists[k] = d;
                sum_d += d;
            }
            double m_sd = sum_d / n_pts;
            double v_sd = 0.0;
            for (int k = 0; k < n_pts; k++) v_sd += (sdists[k] - m_sd) * (sdists[k] - m_sd);
            double sd_s = std::sqrt(v_sd / (n_pts > 1 ? (n_pts - 1) : 1));
            pcv[i] = (sum_d * sd_s) / p;
        }
    }

    return DataFrame::create(
        Named("x") = wrap(mass_x),
        Named("y") = wrap(mass_y),
        Named("area") = wrap(area),
        Named("area_ch") = wrap(area_ch),
        Named("perimeter") = wrap(perimeter),
        Named("radius_mean") = wrap(radius_mean),
        Named("radius_min") = wrap(radius_min),
        Named("radius_max") = wrap(radius_max),
        Named("radius_sd") = wrap(radius_sd),
        Named("radius_ratio") = wrap(radius_ratio),
        Named("diam_mean") = wrap(diam_mean),
        Named("diam_min") = wrap(diam_min),
        Named("diam_max") = wrap(diam_max),
        Named("caliper") = wrap(caliper),
        Named("length") = wrap(length_m),
        Named("width") = wrap(width_m),
        Named("solidity") = wrap(solidity),
        Named("convexity") = wrap(convexity),
        Named("elongation") = wrap(elongation),
        Named("circularity") = wrap(circularity),
        Named("circularity_haralick") = wrap(circularity_haralick),
        Named("circularity_norm") = wrap(circularity_norm),
        Named("eccentricity") = wrap(eccentricity),
        Named("maj_axis") = wrap(maj_axis),
        Named("min_axis") = wrap(min_axis),
        Named("theta") = wrap(theta),
        Named("pcv") = wrap(pcv)
    );
}

// [[Rcpp::export]]
List compute_chulls_cpp(List contours) {
  int n_contours = contours.size();
  List out_list(n_contours);
  CharacterVector col_names = CharacterVector::create("x", "y");

  struct HullData {
    std::vector<double> x;
    std::vector<double> y;
    bool valid = false;
  };
  std::vector<HullData> hulls(n_contours);

  struct ContourInput {
    std::vector<double> x;
    std::vector<double> y;
    bool valid = false;
  };
  std::vector<ContourInput> inputs(n_contours);

  for (int i = 0; i < n_contours; i++) {
    SEXP curr = contours[i];
    if (Rf_isNull(curr)) continue;
    NumericMatrix mat = as<NumericMatrix>(curr);
    int n_pts = mat.nrow();
    if (n_pts < 3) continue;
    inputs[i].x.resize(n_pts);
    inputs[i].y.resize(n_pts);
    const double* ptr = REAL(mat);
    for (int j = 0; j < n_pts; j++) {
      inputs[i].x[j] = ptr[j];
      inputs[i].y[j] = ptr[j + n_pts];
    }
    inputs[i].valid = true;
  }

  #pragma omp parallel for schedule(dynamic) if(n_contours > 10)
  for (int i = 0; i < n_contours; i++) {
    if (!inputs[i].valid) continue;
    const auto& px = inputs[i].x;
    const auto& py = inputs[i].y;
    int n_pts = px.size();

    std::vector<std::pair<double, double>> pts(n_pts);
    for (int j = 0; j < n_pts; j++) pts[j] = {px[j], py[j]};

    std::sort(pts.begin(), pts.end(), [](const std::pair<double, double>& u, const std::pair<double, double>& v) {
      return u.first < v.first || (u.first == v.first && u.second < v.second);
    });

    std::vector<std::pair<double, double>> hull(2 * n_pts);
    int k = 0;
    for (int j = 0; j < n_pts; ++j) {
      while (k >= 2) {
        double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                       (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
        if (cross <= 0) k--; else break;
      }
      hull[k++] = pts[j];
    }
    for (int j = n_pts - 2, t = k + 1; j >= 0; j--) {
      while (k >= t) {
        double cross = (hull[k-1].first - hull[k-2].first) * (pts[j].second - hull[k-2].second) -
                       (hull[k-1].second - hull[k-2].second) * (pts[j].first - hull[k-2].first);
        if (cross <= 0) k--; else break;
      }
      hull[k++] = pts[j];
    }

    int n_ch = k - 1;
    if (n_ch > 0) {
      hulls[i].x.resize(n_ch + 1);
      hulls[i].y.resize(n_ch + 1);
      for (int j = 0; j < n_ch; j++) {
        hulls[i].x[j] = hull[j].first;
        hulls[i].y[j] = hull[j].second;
      }
      hulls[i].x[n_ch] = hull[0].first;
      hulls[i].y[n_ch] = hull[0].second;
      hulls[i].valid = true;
    }
  }

  for (int i = 0; i < n_contours; i++) {
    if (hulls[i].valid) {
      int n_ch = hulls[i].x.size();
      NumericMatrix res(n_ch, 2);
      double* rptr = REAL(res);
      for (int j = 0; j < n_ch; j++) {
        rptr[j] = hulls[i].x[j];
        rptr[j + n_ch] = hulls[i].y[j];
      }
      colnames(res) = col_names;
      out_list[i] = res;
    }
  }
  return out_list;
}

// [[Rcpp::export]]
DataFrame poly_measures_minimal_cpp(List contours) {
    int n_contours = contours.size();
    
    NumericVector mass_x(n_contours, NA_REAL);
    NumericVector mass_y(n_contours);
    NumericVector area(n_contours);
    NumericVector perimeter(n_contours);
    NumericVector length_m(n_contours);
    NumericVector width_m(n_contours);
    NumericVector circularity_norm(n_contours);
    NumericVector eccentricity(n_contours);
    NumericVector maj_axis(n_contours);
    NumericVector min_axis(n_contours);

    for (int i = 0; i < n_contours; i++) {
        SEXP curr_contour = contours[i];
        if (Rf_isNull(curr_contour)) continue;
        
        NumericMatrix coord = as<NumericMatrix>(curr_contour);
        arma::mat C = as<arma::mat>(coord);
        int n_pts = C.n_rows;
        if (n_pts < 3) continue;

        // 1. mass
        arma::vec cm = centmass(C);
        mass_x[i] = cm(0);
        mass_y[i] = cm(1);

        // 2. Area
        double a = 0;
        for (int j = 0; j < n_pts; ++j) {
            int next_j = (j + 1) % n_pts;
            a += C(j, 0) * C(next_j, 1) - C(next_j, 0) * C(j, 1);
        }
        area[i] = std::abs(a / 2.0);

        // 5. Perimeter (distpts)
        double p = 0;
        for (int j = 0; j < n_pts - 1; j++) {
            double dx = C(j+1, 0) - C(j, 0);
            double dy = C(j+1, 1) - C(j, 1);
            p += std::sqrt(dx*dx + dy*dy);
        }
        perimeter[i] = p;

        // 7. Length and Width (from help_lw)
        arma::mat covmat = arma::cov(C);
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covmat);
        arma::mat rotated = arma::flipud(C * eigvec);
        double l = arma::max(rotated.col(1)) - arma::min(rotated.col(1));
        double w = arma::max(rotated.col(0)) - arma::min(rotated.col(0));
        length_m[i] = l;
        width_m[i] = w;

        // 8. Shape factors
        circularity_norm[i] = (area[i] * 4.0 * datum::pi) / (p * p);
        
        // Eigenvalues from population covariance for exact match with help_moments
        double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
        for (int j = 0; j < n_pts; j++) {
            double x = C(j, 0);
            double y = C(j, 1);
            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_y2 += y * y;
            sum_xy += x * y;
        }
        double cov_xy = sum_xy / n_pts - sum_x * sum_y / n_pts / n_pts;
        double var_x = sum_x2 / n_pts - sum_x * sum_x / n_pts / n_pts;
        double var_y = sum_y2 / n_pts - sum_y * sum_y / n_pts / n_pts;
        double a_maj = std::sqrt(0.5 * (var_x + var_y + std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        double b_min = std::sqrt(0.5 * (var_x + var_y - std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        
        maj_axis[i] = std::fmax(a_maj, b_min);
        min_axis[i] = std::fmin(a_maj, b_min);
        eccentricity[i] = std::sqrt(1 - std::pow(min_axis[i] / maj_axis[i], 2));
    }

    return DataFrame::create(
        Named("x") = mass_x,
        Named("y") = mass_y,
        Named("area") = area,
        Named("perimeter") = perimeter,
        Named("length") = length_m,
        Named("width") = width_m,
        Named("circularity_norm") = circularity_norm,
        Named("eccentricity") = eccentricity,
        Named("maj_axis") = maj_axis,
        Named("min_axis") = min_axis
    );
}

// [[Rcpp::export]]
DataFrame poly_measures_disease_cpp(List contours) {
    int n_contours = contours.size();
    
    NumericVector mass_x(n_contours, NA_REAL);
    NumericVector mass_y(n_contours);
    NumericVector area(n_contours);
    NumericVector perimeter(n_contours);
    NumericVector radius_mean(n_contours);
    NumericVector radius_min(n_contours);
    NumericVector radius_max(n_contours);
    NumericVector radius_sd(n_contours);
    NumericVector diam_mean(n_contours);
    NumericVector diam_min(n_contours);
    NumericVector diam_max(n_contours);
    NumericVector length_m(n_contours);
    NumericVector width_m(n_contours);
    NumericVector maj_axis(n_contours);
    NumericVector min_axis(n_contours);

    for (int i = 0; i < n_contours; i++) {
        SEXP curr_contour = contours[i];
        if (Rf_isNull(curr_contour)) continue;
        
        NumericMatrix coord = as<NumericMatrix>(curr_contour);
        arma::mat C = as<arma::mat>(coord);
        int n_pts = C.n_rows;
        if (n_pts < 3) continue;

        // 1. mass
        arma::vec cm = centmass(C);
        mass_x[i] = cm(0);
        mass_y[i] = cm(1);

        // 2. Area
        double a = 0;
        for (int j = 0; j < n_pts; ++j) {
            int next_j = (j + 1) % n_pts;
            a += C(j, 0) * C(next_j, 1) - C(next_j, 0) * C(j, 1);
        }
        area[i] = std::abs(a / 2.0);

        // 5. Perimeter (distpts)
        double p = 0;
        for (int j = 0; j < n_pts - 1; j++) {
            double dx = C(j+1, 0) - C(j, 0);
            double dy = C(j+1, 1) - C(j, 1);
            p += std::sqrt(dx*dx + dy*dy);
        }
        perimeter[i] = p;

        // 6. Centdist
        double cent_x = 0, cent_y = 0;
        for (int j = 0; j < n_pts; j++) {
            cent_x += C(j, 0);
            cent_y += C(j, 1);
        }
        cent_x /= n_pts;
        cent_y /= n_pts;

        arma::vec cdists(n_pts);
        double c_mean = 0, c_min = 1e9, c_max = -1e9;
        for (int j = 0; j < n_pts; j++) {
            double dx = C(j, 0) - cent_x;
            double dy = C(j, 1) - cent_y;
            double d = std::sqrt(dx*dx + dy*dy);
            cdists[j] = d;
            c_mean += d;
            if (d < c_min) c_min = d;
            if (d > c_max) c_max = d;
        }
        c_mean /= n_pts;

        double c_var = 0;
        for (int j = 0; j < n_pts; j++) {
            c_var += (cdists[j] - c_mean) * (cdists[j] - c_mean);
        }
        double c_sd = std::sqrt(c_var / (n_pts - 1));

        radius_mean[i] = c_mean;
        radius_min[i] = c_min;
        radius_max[i] = c_max;
        radius_sd[i] = c_sd;
        diam_mean[i] = c_mean * 2.0;
        diam_min[i] = c_min * 2.0;
        diam_max[i] = c_max * 2.0;

        // 7. Length and Width (from help_lw)
        arma::mat covmat = arma::cov(C);
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covmat);
        arma::mat rotated = arma::flipud(C * eigvec);
        double l = arma::max(rotated.col(1)) - arma::min(rotated.col(1));
        double w = arma::max(rotated.col(0)) - arma::min(rotated.col(0));
        length_m[i] = l;
        width_m[i] = w;
        
        // Eigenvalues from population covariance for exact match with help_moments
        double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
        for (int j = 0; j < n_pts; j++) {
            double x = C(j, 0);
            double y = C(j, 1);
            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_y2 += y * y;
            sum_xy += x * y;
        }
        double cov_xy = sum_xy / n_pts - sum_x * sum_y / n_pts / n_pts;
        double var_x = sum_x2 / n_pts - sum_x * sum_x / n_pts / n_pts;
        double var_y = sum_y2 / n_pts - sum_y * sum_y / n_pts / n_pts;
        double a_maj = std::sqrt(0.5 * (var_x + var_y + std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        double b_min = std::sqrt(0.5 * (var_x + var_y - std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        
        maj_axis[i] = std::fmax(a_maj, b_min);
        min_axis[i] = std::fmin(a_maj, b_min);
    }

    return DataFrame::create(
        Named("x") = mass_x,
        Named("y") = mass_y,
        Named("area") = area,
        Named("perimeter") = perimeter,
        Named("radius_mean") = radius_mean,
        Named("radius_min") = radius_min,
        Named("radius_max") = radius_max,
        Named("radius_sd") = radius_sd,
        Named("diam_mean") = diam_mean,
        Named("diam_min") = diam_min,
        Named("diam_max") = diam_max,
        Named("length") = length_m,
        Named("width") = width_m,
        Named("maj_axis") = maj_axis,
        Named("min_axis") = min_axis
    );
}

// Function to check if a point is inside a polygon
bool pointInPolygon(NumericMatrix polygon, double x, double y) {
  int n = polygon.nrow();
  bool inside = false;
  for (int i = 0, j = n-1; i < n; j = i++) {
    if (((polygon(i, 1) > y) != (polygon(j, 1) > y)) &&
        (x < (polygon(j, 0) - polygon(i, 0)) * (y - polygon(i, 1)) / (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))) {
      inside = !inside;
    }
  }
  return inside;
}

// Function to create a binary image from a polygon
// Input: a matrix containing the coordinates of the polygon vertices (each row is a vertex)
// Output: a binary matrix where the polygon is marked as 1 and the remaining area is marked as 0
// Note: the function assumes that the polygon is closed (i.e., the first and last vertices are the same)
// [[Rcpp::export]]
LogicalMatrix polygon_to_binary(NumericMatrix polygon) {
  int xmin = floor(min(polygon(_, 0)));  // Compute minimum x-coordinate
  int ymin = floor(min(polygon(_, 1)));  // Compute minimum y-coordinate
  int xmax = ceil(max(polygon(_, 0)));   // Compute maximum x-coordinate
  int ymax = ceil(max(polygon(_, 1)));   // Compute maximum y-coordinate
  int width = xmax - xmin + 1;           // Compute width
  int height = ymax - ymin + 1;          // Compute height

  LogicalMatrix binaryImage(width, height);

  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      if (pointInPolygon(polygon, i + xmin, j + ymin)) {
        binaryImage(i, j) = true;
      } else {
        binaryImage(i, j) = false;
      }
    }
  }
  return binaryImage;
}



// [[Rcpp::export]]
IntegerVector sum_true_cols(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  IntegerVector col_sums(ncol);

  for (int j = 0; j < ncol; j++) {
    int col_sum = 0;
    for (int i = 0; i < nrow; i++) {
      if (x(i, j) == 1) {
        col_sum++;
      }
    }
    col_sums[j] = col_sum;
  }
  return col_sums;
}


// [[Rcpp::export]]
NumericVector help_poly_angles(NumericMatrix coords) {
  int num_vertices = coords.nrow();
  NumericMatrix sides(num_vertices, num_vertices);
  NumericVector angles(num_vertices);

  // Calculate the sides of the polygon using the distance formula
  for (int i = 0; i < num_vertices; i++) {
    for (int j = 0; j < num_vertices; j++) {
      sides(i,j) = sqrt(pow(coords(i,0) - coords(j,0), 2) + pow(coords(i,1) - coords(j,1), 2));
    }
  }
  // Calculate the internal angles of the polygon using the law of cosines
  for (int i = 0; i < num_vertices; i++) {
    int prev_vertex = (i == 0) ? num_vertices - 1 : i - 1;
    int next_vertex = (i == num_vertices - 1) ? 0 : i + 1;
    angles(i) = acos((pow(sides(prev_vertex,i), 2) + pow(sides(next_vertex,i), 2) - pow(sides(prev_vertex,next_vertex), 2)) /
      (2 * sides(prev_vertex,i) * sides(next_vertex,i)));
  }

  // Convert the angles from radians to degrees
  angles = angles * 180 / M_PI;
  return angles;
}


// [[Rcpp::export]]
NumericMatrix help_smoth(NumericMatrix coords, int niter) {
  int p = coords.nrow();
  NumericMatrix current = clone(coords);
  NumericMatrix next(p, 2);

  for (int a = 0; a < niter; a++) {
    for (int i = 0; i < p; i++) {
      int prevIndex = (i == 0) ? (p - 1) : (i - 1);
      int nextIndex = (i == p - 1) ? 0 : (i + 1);

      next(i, 0) = (current(i, 0) + current(prevIndex, 0) + current(nextIndex, 0)) / 3.0;
      next(i, 1) = (current(i, 1) + current(prevIndex, 1) + current(nextIndex, 1)) / 3.0;
    }

    for (int i = 0; i < p; i++) {
      current(i, 0) = next(i, 0);
      current(i, 1) = next(i, 1);
    }
  }

  return current;
}

// [[Rcpp::export]]
List smoothContours(List contours, int window_size = 3) {
  if (window_size < 1) {
    stop("Window size must be at least 1.");
  }

  int half_window = window_size / 2;
  List smoothedContours(contours.size());

  for (int i = 0; i < contours.size(); ++i) {
    NumericMatrix contour = contours[i];
    int n = contour.nrow();
    NumericMatrix smoothedContour(n, 2);

    for (int j = 0; j < n; ++j) {
      // Smoothing X and Y separately
      double x_sum = 0.0, y_sum = 0.0;
      int count = 0;

      for (int k = -half_window; k <= half_window; ++k) {
        int idx = j + k;
        if (idx >= 0 && idx < n) {
          x_sum += contour(idx, 0);
          y_sum += contour(idx, 1);
          count++;
        }
      }

      smoothedContour(j, 0) = x_sum / count;
      smoothedContour(j, 1) = y_sum / count;
    }

    smoothedContours[i] = smoothedContour;
  }

  return smoothedContours;
}


// [[Rcpp::export]]
List efourier_cpp(List coords_list, int nharm) {
  int n_contours = coords_list.size();
  List results(n_contours);

  // Convert Rcpp::List of matrices to pure C++ types for thread-safe parallel processing
  std::vector<std::vector<double>> coords_x(n_contours);
  std::vector<std::vector<double>> coords_y(n_contours);
  std::vector<int> nrows(n_contours);
  std::vector<int> nharm_vec(n_contours);

  for (int i = 0; i < n_contours; ++i) {
    NumericMatrix coord = coords_list[i];
    int nr = coord.nrow();
    nrows[i] = nr;
    
    // Adjust nharm for this specific contour if necessary
    int nh = nharm;
    if (nh * 2 > nr) {
      nh = nr / 2;
    }
    if (nh == -1) {
      nh = nr / 2;
    }
    nharm_vec[i] = nh;

    coords_x[i].resize(nr);
    coords_y[i].resize(nr);
    for (int r = 0; r < nr; ++r) {
      coords_x[i][r] = coord(r, 0);
      coords_y[i][r] = coord(r, 1);
    }
  }

  // Pre-allocate pure C++ result structures
  std::vector<std::vector<double>> an_results(n_contours);
  std::vector<std::vector<double>> bn_results(n_contours);
  std::vector<std::vector<double>> cn_results(n_contours);
  std::vector<std::vector<double>> dn_results(n_contours);
  std::vector<double> a0_results(n_contours);
  std::vector<double> c0_results(n_contours);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_contours; ++i) {
    int nr = nrows[i];
    int nh = nharm_vec[i];

    an_results[i].resize(nh, 0.0);
    bn_results[i].resize(nh, 0.0);
    cn_results[i].resize(nh, 0.0);
    dn_results[i].resize(nh, 0.0);

    if (nr < 3 || nh <= 0) {
      continue;
    }

    const std::vector<double>& cx = coords_x[i];
    const std::vector<double>& cy = coords_y[i];

    std::vector<double> Dx(nr);
    std::vector<double> Dy(nr);
    std::vector<double> Dt(nr);
    std::vector<double> t1(nr);
    std::vector<double> t1m1(nr);

    double cum_sum = 0.0;
    for (int r = 0; r < nr; ++r) {
      int prev_r = (r == 0) ? (nr - 1) : (r - 1);
      Dx[r] = cx[r] - cx[prev_r];
      Dy[r] = cy[r] - cy[prev_r];
      double dt_val = std::sqrt(Dx[r]*Dx[r] + Dy[r]*Dy[r]);
      if (dt_val < 1e-10) {
        dt_val = 1e-10;
      }
      Dt[r] = dt_val;
      cum_sum += dt_val;
      t1[r] = cum_sum;
      t1m1[r] = (r == 0) ? 0.0 : t1[r-1];
    }
    double t = cum_sum;

    for (int h = 1; h <= nh; ++h) {
      double ti = t / (2.0 * datum::pi * datum::pi * h * h);
      double r_val = 2.0 * h * datum::pi;

      double an_sum = 0.0;
      double bn_sum = 0.0;
      double cn_sum = 0.0;
      double dn_sum = 0.0;

      for (int r = 0; r < nr; ++r) {
        double dx_dt = Dx[r] / Dt[r];
        double dy_dt = Dy[r] / Dt[r];

        double angle_t1 = r_val * t1[r] / t;
        double angle_t1m1 = r_val * t1m1[r] / t;

        double diff_cos = std::cos(angle_t1) - std::cos(angle_t1m1);
        double diff_sin = std::sin(angle_t1) - std::sin(angle_t1m1);

        an_sum += dx_dt * diff_cos;
        bn_sum += dx_dt * diff_sin;
        cn_sum += dy_dt * diff_cos;
        dn_sum += dy_dt * diff_sin;
      }

      an_results[i][h - 1] = ti * an_sum;
      bn_results[i][h - 1] = ti * bn_sum;
      cn_results[i][h - 1] = ti * cn_sum;
      dn_results[i][h - 1] = ti * dn_sum;
    }

    double a0_sum = 0.0;
    double c0_sum = 0.0;
    for (int r = 0; r < nr; ++r) {
      a0_sum += cx[r] * Dt[r];
      c0_sum += cy[r] * Dt[r];
    }
    a0_results[i] = 2.0 * a0_sum / t;
    c0_results[i] = 2.0 * c0_sum / t;
  }

  // Convert pure C++ results back to Rcpp structures (single-threaded)
  for (int i = 0; i < n_contours; ++i) {
    int nh = nharm_vec[i];
    int nr = nrows[i];

    NumericVector an(nh);
    NumericVector bn(nh);
    NumericVector cn(nh);
    NumericVector dn(nh);
    for (int h = 0; h < nh; ++h) {
      an[h] = an_results[i][h];
      bn[h] = bn_results[i][h];
      cn[h] = cn_results[i][h];
      dn[h] = dn_results[i][h];
    }

    // Prepare coords matrix to return
    NumericMatrix coords_mat(nr, 2);
    for (int r = 0; r < nr; ++r) {
      coords_mat(r, 0) = coords_x[i][r];
      coords_mat(r, 1) = coords_y[i][r];
    }

    List coefs = List::create(
      Named("an") = an,
      Named("bn") = bn,
      Named("cn") = cn,
      Named("dn") = dn,
      Named("a0") = a0_results[i],
      Named("c0") = c0_results[i],
      Named("nr") = nr,
      Named("nharm") = nh,
      Named("coords") = coords_mat
    );
    coefs.attr("class") = "efourier";

    results[i] = coefs;
  }

  return results;
}


// Helper function to replicate R's %% operator behavior
inline double r_mod(double x, double y) {
  double r = std::fmod(x, y);
  return r < 0.0 ? r + y : r;
}

// [[Rcpp::export]]
List efourier_norm_cpp(List efourier_list, bool start) {
  int n_objects = efourier_list.size();
  List results(n_objects);

  // Convert Rcpp structures to pure C++ types for thread safety
  std::vector<std::vector<double>> an_vec(n_objects);
  std::vector<std::vector<double>> bn_vec(n_objects);
  std::vector<std::vector<double>> cn_vec(n_objects);
  std::vector<std::vector<double>> dn_vec(n_objects);
  std::vector<double> a0_vec(n_objects);
  std::vector<double> c0_vec(n_objects);
  std::vector<int> nharm_vec(n_objects);

  for (int i = 0; i < n_objects; ++i) {
    List obj = efourier_list[i];
    NumericVector an = obj["an"];
    NumericVector bn = obj["bn"];
    NumericVector cn = obj["cn"];
    NumericVector dn = obj["dn"];
    
    int nh = an.size();
    nharm_vec[i] = nh;
    a0_vec[i] = obj["a0"];
    c0_vec[i] = obj["c0"];

    an_vec[i].assign(an.begin(), an.end());
    bn_vec[i].assign(bn.begin(), bn.end());
    cn_vec[i].assign(cn.begin(), cn.end());
    dn_vec[i].assign(dn.begin(), dn.end());
  }

  // Pre-allocate C++ results
  std::vector<std::vector<double>> A_results(n_objects);
  std::vector<std::vector<double>> B_results(n_objects);
  std::vector<std::vector<double>> C_results(n_objects);
  std::vector<std::vector<double>> D_results(n_objects);
  std::vector<double> size_results(n_objects);
  std::vector<double> theta_results(n_objects);
  std::vector<double> psi_results(n_objects);
  std::vector<std::vector<double>> lnef_results(n_objects);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_objects; ++i) {
    int nh = nharm_vec[i];
    if (nh <= 0) continue;

    A_results[i].resize(nh, 0.0);
    B_results[i].resize(nh, 0.0);
    C_results[i].resize(nh, 0.0);
    D_results[i].resize(nh, 0.0);
    lnef_results[i].resize(4, 0.0);

    const std::vector<double>& an = an_vec[i];
    const std::vector<double>& bn = bn_vec[i];
    const std::vector<double>& cn = cn_vec[i];
    const std::vector<double>& dn = dn_vec[i];

    double A1 = an[0];
    double B1 = bn[0];
    double C1 = cn[0];
    double D1 = dn[0];

    // Compute theta
    double num = 2.0 * (A1 * B1 + C1 * D1);
    double den = A1 * A1 + C1 * C1 - B1 * B1 - D1 * D1;
    double theta_val = 0.0;
    if (std::abs(den) < 1e-15) {
      theta_val = (num >= 0.0) ? (datum::pi / 2.0) : (-datum::pi / 2.0);
    } else {
      theta_val = std::atan(num / den);
    }
    double theta_r = 0.5 * theta_val;
    theta_r = r_mod(theta_r, datum::pi);

    // Compute phaseshift and M2
    double cos_t = std::cos(theta_r);
    double sin_t = std::sin(theta_r);

    double M2_00 = A1 * cos_t + B1 * sin_t;
    double M2_01 = -A1 * sin_t + B1 * cos_t;
    double M2_10 = C1 * cos_t + D1 * sin_t;
    double M2_11 = -C1 * sin_t + D1 * cos_t;

    double v0 = M2_00*M2_00 + M2_10*M2_10;
    double v1 = M2_01*M2_01 + M2_11*M2_11;

    if (v0 < v1) {
      theta_r += datum::pi / 2.0;
    }
    theta_r = r_mod(theta_r + datum::pi / 2.0, datum::pi) - datum::pi / 2.0;

    double Aa = A1 * std::cos(theta_r) + B1 * std::sin(theta_r);
    double Cc = C1 * std::cos(theta_r) + D1 * std::sin(theta_r);
    double scale = std::sqrt(Aa*Aa + Cc*Cc);
    double psi = r_mod(std::atan(Cc / Aa), datum::pi);
    if (Aa < 0.0) {
      psi += datum::pi;
    }
    double size = 1.0 / scale;

    double cos_p = std::cos(psi);
    double sin_p = std::sin(psi);

    if (start) {
      theta_r = 0.0;
    }

    for (int h = 1; h <= nh; ++h) {
      double cos_ht = std::cos(h * theta_r);
      double sin_ht = std::sin(h * theta_r);

      // Temp matrix multiplication:
      // Temp = [[an[h-1], bn[h-1]], [cn[h-1], dn[h-1]]] %*% [[cos_ht, -sin_ht], [sin_ht, cos_ht]]
      double T00 = an[h-1] * cos_ht + bn[h-1] * sin_ht;
      double T01 = -an[h-1] * sin_ht + bn[h-1] * cos_ht;
      double T10 = cn[h-1] * cos_ht + dn[h-1] * sin_ht;
      double T11 = -cn[h-1] * sin_ht + dn[h-1] * cos_ht;

      // mat = size * rotation %*% Temp
      // rotation = [[cos_p, sin_p], [-sin_p, cos_p]]
      double A_val = size * (cos_p * T00 + sin_p * T10);
      double B_val = size * (cos_p * T01 + sin_p * T11);
      double C_val = size * (-sin_p * T00 + cos_p * T10);
      double D_val = size * (-sin_p * T01 + cos_p * T11);

      A_results[i][h - 1] = A_val;
      B_results[i][h - 1] = B_val;
      C_results[i][h - 1] = C_val;
      D_results[i][h - 1] = D_val;

      lnef_results[i][0] = A_val;
      lnef_results[i][1] = B_val;
      lnef_results[i][2] = C_val;
      lnef_results[i][3] = D_val;
    }

    size_results[i] = scale;
    theta_results[i] = theta_r;
    psi_results[i] = psi;
  }

  // Convert results to Rcpp (single-threaded)
  for (int i = 0; i < n_objects; ++i) {
    int nh = nharm_vec[i];
    NumericVector A(nh);
    NumericVector B(nh);
    NumericVector C(nh);
    NumericVector D(nh);

    for (int h = 0; h < nh; ++h) {
      A[h] = A_results[i][h];
      B[h] = B_results[i][h];
      C[h] = C_results[i][h];
      D[h] = D_results[i][h];
    }

    NumericVector lnef(4);
    lnef[0] = lnef_results[i][0];
    lnef[1] = lnef_results[i][1];
    lnef[2] = lnef_results[i][2];
    lnef[3] = lnef_results[i][3];

    List coefs = List::create(
      Named("A") = A,
      Named("B") = B,
      Named("C") = C,
      Named("D") = D,
      Named("size") = size_results[i],
      Named("theta") = theta_results[i],
      Named("psi") = psi_results[i],
      Named("a0") = a0_vec[i],
      Named("c0") = c0_vec[i],
      Named("lnef") = lnef,
      Named("nharm") = nh
    );
    coefs.attr("class") = "nefourier";

    results[i] = coefs;
  }

  return results;
}

// [[Rcpp::export]]
List rfourier_cpp(List contours, IntegerVector nharm_vec) {
  int n_objects = contours.size();
  List results(n_objects);

  std::vector<std::vector<double>> a_results(n_objects);
  std::vector<std::vector<double>> b_results(n_objects);
  std::vector<double> a0_results(n_objects);

  std::vector<arma::mat> mats(n_objects);
  for (int i = 0; i < n_objects; ++i) {
    mats[i] = as<arma::mat>(contours[i]);
    a_results[i].resize(nharm_vec[i]);
    b_results[i].resize(nharm_vec[i]);
  }

  #pragma omp parallel for
  for (int i = 0; i < n_objects; ++i) {
    const arma::mat& C = mats[i];
    int N = C.n_rows;
    int nh = nharm_vec[i];
    
    double xc = 0, yc = 0;
    for (int j = 0; j < N; ++j) {
      xc += C(j, 0);
      yc += C(j, 1);
    }
    xc /= N;
    yc /= N;

    std::vector<double> r(N);
    double a0 = 0;
    for (int j = 0; j < N; ++j) {
      double dx = C(j, 0) - xc;
      double dy = C(j, 1) - yc;
      r[j] = std::sqrt(dx*dx + dy*dy);
      a0 += r[j];
    }
    a0 /= N;
    a0_results[i] = a0;

    for (int h = 1; h <= nh; ++h) {
      double an = 0;
      double bn = 0;
      for (int j = 0; j < N; ++j) {
        double theta = (2.0 * datum::pi * j) / N;
        an += r[j] * std::cos(h * theta);
        bn += r[j] * std::sin(h * theta);
      }
      an *= (2.0 / N);
      bn *= (2.0 / N);
      a_results[i][h-1] = an;
      b_results[i][h-1] = bn;
    }
  }

  for (int i = 0; i < n_objects; ++i) {
    int nh = nharm_vec[i];
    NumericVector a(nh);
    NumericVector b(nh);
    for (int h = 0; h < nh; ++h) {
      a[h] = a_results[i][h];
      b[h] = b_results[i][h];
    }
    List coefs = List::create(
      Named("an") = a,
      Named("bn") = b,
      Named("a0") = a0_results[i]
    );
    results[i] = coefs;
  }
  return results;
}

// [[Rcpp::export]]
List tfourier_cpp(List contours, IntegerVector nharm_vec) {
  int n_objects = contours.size();
  List results(n_objects);

  std::vector<std::vector<double>> a_results(n_objects);
  std::vector<std::vector<double>> b_results(n_objects);

  std::vector<arma::mat> mats(n_objects);
  for (int i = 0; i < n_objects; ++i) {
    mats[i] = as<arma::mat>(contours[i]);
    a_results[i].resize(nharm_vec[i]);
    b_results[i].resize(nharm_vec[i]);
  }

  #pragma omp parallel for
  for (int i = 0; i < n_objects; ++i) {
    const arma::mat& C = mats[i];
    int N = C.n_rows;
    int nh = nharm_vec[i];
    
    std::vector<double> phi(N);
    for (int j = 0; j < N; ++j) {
      int next_j = (j + 1) % N;
      double dx = C(next_j, 0) - C(j, 0);
      double dy = C(next_j, 1) - C(j, 1);
      phi[j] = std::atan2(dy, dx);
    }
    
    for (int j = 1; j < N; ++j) {
      double diff = phi[j] - phi[j-1];
      while (diff > datum::pi) {
        phi[j] -= 2.0 * datum::pi;
        diff = phi[j] - phi[j-1];
      }
      while (diff < -datum::pi) {
        phi[j] += 2.0 * datum::pi;
        diff = phi[j] - phi[j-1];
      }
    }
    
    double phi0 = phi[0];
    for (int h = 1; h <= nh; ++h) {
      double an = 0;
      double bn = 0;
      for (int j = 0; j < N; ++j) {
        double t = (2.0 * datum::pi * j) / N;
        double Phi_t = phi[j] - phi0 - t;
        an += Phi_t * std::cos(h * t);
        bn += Phi_t * std::sin(h * t);
      }
      an *= (2.0 / N);
      bn *= (2.0 / N);
      a_results[i][h-1] = an;
      b_results[i][h-1] = bn;
    }
  }

  for (int i = 0; i < n_objects; ++i) {
    int nh = nharm_vec[i];
    NumericVector a(nh);
    NumericVector b(nh);
    for (int h = 0; h < nh; ++h) {
      a[h] = a_results[i][h];
      b[h] = b_results[i][h];
    }
    List coefs = List::create(
      Named("an") = a,
      Named("bn") = b
    );
    results[i] = coefs;
  }
  return results;
}

// [[Rcpp::export]]
List gpa_cpp(List contours, double tol = 1e-5, int max_iter = 100) {
  int n_objects = contours.size();
  std::vector<arma::mat> shapes(n_objects);
  std::vector<double> centroid_sizes(n_objects);
  
  for (int i = 0; i < n_objects; ++i) {
    arma::mat C = as<arma::mat>(contours[i]);
    arma::rowvec centroid = arma::mean(C, 0);
    C.each_row() -= centroid;
    
    double cs = arma::norm(C, "fro");
    centroid_sizes[i] = cs;
    if (cs > 0) {
      C /= cs;
    }
    shapes[i] = C;
  }
  
  arma::mat ref_shape = shapes[0];
  arma::mat mean_shape = ref_shape;
  
  int iter = 0;
  double diff = 1e9;
  
  while (iter < max_iter && diff > tol) {
    arma::mat sum_shapes = arma::zeros<arma::mat>(ref_shape.n_rows, ref_shape.n_cols);
    
    for (int i = 0; i < n_objects; ++i) {
      arma::mat Y = shapes[i];
      arma::mat U, V;
      arma::vec s;
      arma::svd(U, s, V, Y.t() * ref_shape);
      
      arma::mat R = U * V.t();
      
      if (arma::det(R) < 0) {
        arma::mat U_mod = U;
        U_mod.col(U.n_cols - 1) *= -1;
        R = U_mod * V.t();
      }
      
      arma::mat Y_rot = Y * R;
      shapes[i] = Y_rot;
      sum_shapes += Y_rot;
    }
    
    arma::mat new_mean = sum_shapes / n_objects;
    new_mean.each_row() -= arma::mean(new_mean, 0);
    new_mean /= arma::norm(new_mean, "fro");
    
    diff = arma::norm(new_mean - ref_shape, "fro");
    ref_shape = new_mean;
    iter++;
  }
  
  List aligned_shapes(n_objects);
  for (int i = 0; i < n_objects; ++i) {
    aligned_shapes[i] = wrap(shapes[i]);
  }
  
  return List::create(
    Named("aligned") = aligned_shapes,
    Named("consensus") = wrap(ref_shape),
    Named("centroid_sizes") = wrap(centroid_sizes),
    Named("iterations") = iter
  );
}

// [[Rcpp::export]]
List object_bbox_cpp(List contours) {
  int n = contours.length();
  List bbox_list(n);

  for (int k = 0; k < n; ++k) {
    SEXP item = contours[k];
    if (Rf_isNull(item)) {
      List item_box;
      item_box["x_min"] = NA_REAL;
      item_box["y_min"] = NA_REAL;
      item_box["x_max"] = NA_REAL;
      item_box["y_max"] = NA_REAL;
      bbox_list[k] = item_box;
      continue;
    }

    NumericMatrix coords(item);
    int npts = coords.nrow();
    if (npts == 0) {
      List item_box;
      item_box["x_min"] = NA_REAL;
      item_box["y_min"] = NA_REAL;
      item_box["x_max"] = NA_REAL;
      item_box["y_max"] = NA_REAL;
      bbox_list[k] = item_box;
      continue;
    }

    double x_min = coords(0, 0);
    double x_max = coords(0, 0);
    double y_min = coords(0, 1);
    double y_max = coords(0, 1);

    for (int i = 1; i < npts; ++i) {
      double x = coords(i, 0);
      double y = coords(i, 1);
      if (x < x_min) x_min = x;
      if (x > x_max) x_max = x;
      if (y < y_min) y_min = y;
      if (y > y_max) y_max = y;
    }

    List item_box;
    item_box["x_min"] = x_min;
    item_box["y_min"] = y_min;
    item_box["x_max"] = x_max;
    item_box["y_max"] = y_max;
    bbox_list[k] = item_box;
  }

  return bbox_list;
}
