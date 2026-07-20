#include <RcppArmadillo.h>
#include <math.h>
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
NumericVector get_area_mask(IntegerVector mask) {
  int n = mask.length();
  int max_value = max(mask);
  NumericVector area(max_value);
  std::fill(area.begin(), area.end(), 0);
  for (int i = 0; i < n; i++) {
    int x = mask[i];
    if (x > 0 && x <= max_value) {
      area[x - 1]++;
    }
  }
  return area;
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

// Forward declaration
NumericMatrix help_smoth(NumericMatrix coords, int niter);

// [[Rcpp::export]]
DataFrame poly_measures_cpp(List contours) {
    int n_contours = contours.size();

    NumericVector mass_x(n_contours);
    NumericVector mass_y(n_contours);
    NumericVector area(n_contours);
    NumericVector area_ch(n_contours);
    NumericVector perimeter(n_contours);
    NumericVector radius_mean(n_contours);
    NumericVector radius_min(n_contours);
    NumericVector radius_max(n_contours);
    NumericVector radius_sd(n_contours);
    NumericVector radius_ratio(n_contours);
    NumericVector diam_mean(n_contours);
    NumericVector diam_min(n_contours);
    NumericVector diam_max(n_contours);
    NumericVector caliper(n_contours);
    NumericVector length_m(n_contours);
    NumericVector width_m(n_contours);
    NumericVector solidity(n_contours);
    NumericVector convexity(n_contours);
    NumericVector elongation(n_contours);
    NumericVector circularity(n_contours);
    NumericVector circularity_haralick(n_contours);
    NumericVector circularity_norm(n_contours);
    NumericVector eccentricity(n_contours);
    NumericVector maj_axis(n_contours);
    NumericVector min_axis(n_contours);
    NumericVector theta(n_contours);
    NumericVector pcv(n_contours, NA_REAL);
    
    // Initialize mass_x to NA_REAL to easily filter invalid rows in R
    mass_x.fill(NA_REAL);

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

        // 3. Convex hull
        arma::mat ch = get_convex_hull(C);

        // 4. Area CH
        double a_ch = 0;
        int n_ch = ch.n_rows;
        for (int j = 0; j < n_ch; ++j) {
            int next_j = (j + 1) % n_ch;
            a_ch += ch(j, 0) * ch(next_j, 1) - ch(next_j, 0) * ch(j, 1);
        }
        area_ch[i] = std::abs(a_ch / 2.0);

        // 5. Perimeter (distpts)
        double p = 0;
        for (int j = 0; j < n_pts - 1; j++) {
            double dx = C(j+1, 0) - C(j, 0);
            double dy = C(j+1, 1) - C(j, 1);
            p += std::sqrt(dx*dx + dy*dy);
        }
        perimeter[i] = p;

        // 5.1 Perimeter of Convex Hull (for convexity)
        double p_ch = 0;
        if (n_ch > 0) {
            // Convex hull might not close the loop automatically like contours do
            // Wait, get_convex_hull returns k-1 points, which means the last point is NOT the first point.
            // Monotone chain: hull[0] to hull[k-1]. hull[k-1] == hull[0].
            // I returned k-1 points, so res(0) to res(k-2).
            // To compute perimeter, we need to close the loop!
            for (int j = 0; j < n_ch; j++) {
                int next_j = (j + 1) % n_ch;
                double dx = ch(next_j, 0) - ch(j, 0);
                double dy = ch(next_j, 1) - ch(j, 1);
                p_ch += std::sqrt(dx*dx + dy*dy);
            }
        }

        // 6. Centdist (mean/sd of distances to centroid calculated as simple average)
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
        radius_ratio[i] = c_max / c_min;
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

        caliper[i] = l;
        length_m[i] = l;
        width_m[i] = w;

        // 8. Shape factors
        solidity[i] = area[i] / area_ch[i];
        convexity[i] = p_ch / p;
        elongation[i] = 1.0 - (w / l);
        circularity[i] = (p * p) / area[i];
        circularity_haralick[i] = c_mean / c_sd;
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
        double t = 0.5 * std::atan2(2 * cov_xy, var_x - var_y);
        double a_maj = std::sqrt(0.5 * (var_x + var_y + std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        double b_min = std::sqrt(0.5 * (var_x + var_y - std::sqrt(std::pow(var_x - var_y, 2) + 4 * std::pow(cov_xy, 2))));
        
        maj_axis[i] = std::fmax(a_maj, b_min);
        min_axis[i] = std::fmin(a_maj, b_min);
        eccentricity[i] = std::sqrt(1 - std::pow(min_axis[i] / maj_axis[i], 2));
        theta[i] = t;

        // 9. PCV (poly_pcv)
        // Mathematically correct Jacobi update (symmetric moving average)
        arma::mat smoth = C;
        int niter = 100;
        arma::mat smoothed(n_pts, 2);
        for (int a = 0; a < niter; a++) {
            for (int k = 0; k < n_pts; k++) {
                int prev = (k == 0) ? (n_pts - 1) : (k - 1);
                int next = (k == n_pts - 1) ? 0 : (k + 1);
                smoothed(k, 0) = (smoth(k, 0) + smoth(prev, 0) + smoth(next, 0)) / 3.0;
                smoothed(k, 1) = (smoth(k, 1) + smoth(prev, 1) + smoth(next, 1)) / 3.0;
            }
            smoth = smoothed;
        }
        
        double sum_dists = 0;
        arma::vec smooth_dists(n_pts);
        for (int k = 0; k < n_pts; k++) {
            double dx = C(k, 0) - smoth(k, 0);
            double dy = C(k, 1) - smoth(k, 1);
            double d = std::sqrt(dx*dx + dy*dy);
            smooth_dists[k] = d;
            sum_dists += d;
        }
        double mean_sdists = sum_dists / n_pts;
        double var_sdists = 0;
        for (int k = 0; k < n_pts; k++) {
            var_sdists += (smooth_dists[k] - mean_sdists) * (smooth_dists[k] - mean_sdists);
        }
        double sd_sdists = std::sqrt(var_sdists / (n_pts - 1));

        pcv[i] = (sum_dists * sd_sdists) / p;
    }

    return DataFrame::create(
        Named("x") = mass_x,
        Named("y") = mass_y,
        Named("area") = area,
        Named("area_ch") = area_ch,
        Named("perimeter") = perimeter,
        Named("radius_mean") = radius_mean,
        Named("radius_min") = radius_min,
        Named("radius_max") = radius_max,
        Named("radius_sd") = radius_sd,
        Named("radius_ratio") = radius_ratio,
        Named("diam_mean") = diam_mean,
        Named("diam_min") = diam_min,
        Named("diam_max") = diam_max,
        Named("caliper") = caliper,
        Named("length") = length_m,
        Named("width") = width_m,
        Named("solidity") = solidity,
        Named("convexity") = convexity,
        Named("elongation") = elongation,
        Named("circularity") = circularity,
        Named("circularity_haralick") = circularity_haralick,
        Named("circularity_norm") = circularity_norm,
        Named("eccentricity") = eccentricity,
        Named("maj_axis") = maj_axis,
        Named("min_axis") = min_axis,
        Named("theta") = theta,
        Named("pcv") = pcv
    );
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
