#include <RcppArmadillo.h>
#include <queue>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <random>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// adapted from imagerExtra  https://bit.ly/3HtxumB
Rcpp::NumericMatrix int_sum(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix res(nrow, ncol);

  res(0,0) = mat(0,0);
  for (int i = 1; i < nrow; ++i) {
    res(i,0) = mat(i,0) + res(i-1,0);
  }
  for (int j = 1; j < ncol; ++j) {
    res(0,j) = mat(0,j) + res(0,j-1);
  }
  for (int i = 1; i < nrow; ++i) {
    for (int j = 1; j < ncol; ++j) {
      res(i,j) = mat(i,j) + res(i-1,j) + res(i,j-1) - res(i-1,j-1);
    }
  }
  return res;
}

Rcpp::NumericMatrix int_sum_squared(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix mat_squared(nrow, ncol);
  Rcpp::NumericMatrix res(nrow, ncol);

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      mat_squared(i,j) = mat(i,j) * mat(i,j);
    }
  }

  res(0,0) = mat_squared(0,0);
  for (int i = 1; i < nrow; ++i) {
    res(i,0) = mat_squared(i,0) + res(i-1,0);
  }
  for (int j = 1; j < ncol; ++j) {
    res(0,j) = mat_squared(0,j) + res(0,j-1);
  }
  for (int i = 1; i < nrow; ++i) {
    for (int j = 1; j < ncol; ++j) {
      res(i,j) = mat_squared(i,j) + res(i-1,j) + res(i,j-1) - res(i-1,j-1);
    }
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix threshold_adaptive(Rcpp::NumericMatrix mat, double k, int windowsize, double maxsd) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix res(nrow, ncol);
  Rcpp::NumericMatrix integ_sum = int_sum(mat);
  Rcpp::NumericMatrix int_sum_sqr = int_sum_squared(mat);
  int winhalf = windowsize / 2;
  int winsize_squared = windowsize * windowsize;
  int nrow_center = nrow - windowsize;
  int ncol_center = ncol - windowsize;


  for (int i = 0; i < winhalf; ++i) {
    for (int j = 0; j < winhalf; ++j) {
      int temp_winsize = (winhalf + i + 1) * (winhalf + j + 1);
      double mean_local = integ_sum(i+winhalf,j+winhalf) / temp_winsize;
      double sd_local = sqrt(int_sum_sqr(i+winhalf,j+winhalf) / temp_winsize - mean_local * mean_local);
      double threshold_local  = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j =0; j < winhalf; ++j) {
      int temp_winsize = windowsize * (winhalf + j + 1);
      double mean_local = (integ_sum(i+winhalf,j+winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local  = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = 0; j < winhalf; ++j) {
      int temp_winsize = (winhalf + nrow - i) * (winhalf + j + 1);
      double mean_local = (integ_sum(nrow-1,j+winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,j+winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = 0; i < winhalf; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      int temp_winsize = (winhalf + i + 1) * windowsize;
      double mean_local = (integ_sum(i+winhalf,j+winhalf) - integ_sum(i+winhalf,j-winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) - int_sum_sqr(i+winhalf,j-winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      double mean_local = (integ_sum(i+winhalf,j+winhalf) + integ_sum(i-winhalf,j-winhalf) - integ_sum(i+winhalf,j-winhalf) - integ_sum(i-winhalf,j+winhalf)) / winsize_squared;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(i+winhalf,j-winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / winsize_squared - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      int temp_winsize = (winhalf + nrow - i) * windowsize;
      double mean_local = (integ_sum(nrow-1,j+winhalf) + integ_sum(i-winhalf,j-winhalf) - integ_sum(nrow-1,j-winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,j+winhalf) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(nrow-1,j-winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = 0; i < winhalf; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = (winhalf + i + 1) * (winhalf + ncol - j);
      double mean_local = (integ_sum(i+winhalf,ncol-1) - integ_sum(i+winhalf,j-winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,ncol-1) - int_sum_sqr(i+winhalf,j-winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = windowsize * (winhalf + ncol - j);
      double mean_local = (integ_sum(i+winhalf,ncol-1) + integ_sum(i-winhalf,j-winhalf) - integ_sum(i+winhalf,j-winhalf) - integ_sum(i-winhalf,ncol-1)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,ncol-1) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(i+winhalf,j-winhalf) - int_sum_sqr(i-winhalf,ncol-1)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = (winhalf + nrow - i) * (winhalf + ncol - j);
      double mean_local = (integ_sum(nrow-1,ncol-1) + integ_sum(i-winhalf,j-winhalf) - integ_sum(nrow-1,j-winhalf) - integ_sum(i-winhalf,ncol-1)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,ncol-1) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(nrow-1,j-winhalf) - int_sum_sqr(i-winhalf,ncol-1)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }
  return res;
}


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
NumericMatrix rgb_to_hsb_help(NumericVector r, NumericVector g, NumericVector b) {
  NumericMatrix hsb(r.size(), 3);
  for (int i = 0; i < r.size(); i++) {
    double max_val = std::max(std::max(r[i], g[i]), b[i]);
    double min_val = std::min(std::min(r[i], g[i]), b[i]);
    double diff = max_val - min_val;
    if (max_val == r[i]) {
      hsb(i, 0) = 60 * ((g[i] - b[i]) / diff);
    } else if (max_val == g[i]) {
      hsb(i, 0) = 60 * (2 + (b[i] - r[i]) / diff);
    } else {
      hsb(i, 0) = 60 * (4 + (r[i] - g[i]) / diff);
    }
    hsb(i, 1) = (max_val - min_val) / max_val * 100;
    hsb(i, 2) = max_val * 100;
  }
  return hsb;
}


// [[Rcpp::export]]
arma::mat rgb_to_srgb_help(const arma::mat& rgb) {
  double gamma = 2.2;
  arma::mat rgb_gamma = pow(rgb, gamma);

  arma::mat matrix = {
    { 3.2406, -1.5372, -0.4986 },
    { -0.9689, 1.8758, 0.0415 },
    { 0.0557, -0.2040, 1.0570 }
  };
  arma::mat rgb_srgb = rgb_gamma * matrix;

  rgb_srgb.elem(find(rgb_srgb < 0)).zeros();
  rgb_srgb.elem(find(rgb_srgb > 1)).ones();
  return rgb_srgb;
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


// DISTANCE TRANSFORM
// [[Rcpp::export]]
NumericMatrix help_dist_transform(const LogicalMatrix &bin) {
  int nrow = bin.nrow(), ncol = bin.ncol();
  NumericMatrix dist(nrow, ncol);
  std::fill(dist.begin(), dist.end(), -1.0);
  std::queue<std::pair<int, int> > q;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (bin(i, j)) {
        dist(i, j) = 0.0;
        q.push(std::make_pair(i, j));
      }
    }
  }
  int dx[] = {-1, 0, 1, 0};
  int dy[] = {0, 1, 0, -1};
  while (!q.empty()) {
    int x = q.front().first, y = q.front().second;
    q.pop();
    for (int d = 0; d < 4; d++) {
      int nx = x + dx[d], ny = y + dy[d];
      if (nx >= 0 && nx < nrow && ny >= 0 && ny < ncol && dist(nx, ny) < 0) {
        dist(nx, ny) = dist(x, y) + 1.0;
        q.push(std::make_pair(nx, ny));
      }
    }
  }
  return dist;
}


// WATERSHED
// [[Rcpp::export]]
IntegerMatrix help_watershed(IntegerMatrix binary, IntegerMatrix labels, IntegerMatrix distances) {
  int rows = binary.nrow();
  int cols = binary.ncol();

  // create the output matrix with zeros
  IntegerMatrix output(rows, cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      output(i, j) = 0;
    }
  }
  // find the location of all markers and store them in a queue
  std::queue<std::pair<int, int>> markers;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (labels(i, j) != 0) {
        markers.push(std::make_pair(i, j));
      }
    }
  }
  // loop through the markers and perform the watershed segmentation
  while (!markers.empty()) {
    std::pair<int, int> current = markers.front();
    markers.pop();
    int x = current.first;
    int y = current.second;
    int currentLabel = labels(x, y);
    if (x > 0 && binary(x-1, y) == 1 && labels(x-1, y) == 0) {
      labels(x-1, y) = currentLabel;
      markers.push(std::make_pair(x-1, y));
    }
    if (x < rows - 1 && binary(x+1, y) == 1 && labels(x+1, y) == 0) {
      labels(x+1, y) = currentLabel;
      markers.push(std::make_pair(x+1, y));
    }
    if (y > 0 && binary(x, y-1) == 1 && labels(x, y-1) == 0) {
      labels(x, y-1) = currentLabel;markers.push(std::make_pair(x, y-1));
    }
    if (y < cols - 1 && binary(x, y+1) == 1 && labels(x, y+1) == 0) {
      labels(x, y+1) = currentLabel;
      markers.push(std::make_pair(x, y+1));
    }
  }
  // perform the transform distance calculation
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (binary(i, j) == 1) {
        int minDistance = INT_MAX;
        int minLabel = -1;
        if (i > 0 && labels(i-1, j) != 0 && distances(i-1, j) < minDistance) {
          minDistance = distances(i-1, j);
          minLabel = labels(i-1, j);
        }
        if (i < rows - 1 && labels(i+1, j) != 0 && distances(i+1, j) < minDistance) {
          minDistance = distances(i+1, j);
          minLabel = labels(i+1, j);
        }
        if (j > 0 && labels(i, j-1) != 0 && distances(i, j-1) < minDistance) {
          minDistance = distances(i, j-1);
          minLabel = labels(i, j-1);
        }
        if (j < cols - 1 && labels(i, j+1) != 0 && distances(i, j+1) < minDistance) {
          minDistance = distances(i, j+1);
          minLabel = labels(i, j+1);
        }
        if (minLabel != -1) output(i, j) = minLabel;
      }
    }
  }
  return output;
}


// EXTRACT PIXELS
// [[Rcpp::export]]
std::vector<std::vector<double>> help_get_rgb(const NumericMatrix &R, const NumericMatrix &G, const NumericMatrix &B, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = R.nrow();
  int ncol = R.ncol();
  // get the number of labels
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  // create a list to store the RGB values for each label
  std::vector<std::vector<double>> result(labelsCount);

  // loop through each label
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int label = labels(i, j);
      if (label > 0) {
        result[label].push_back(label);
        result[label].push_back(R(i, j));
        result[label].push_back(G(i, j));
        result[label].push_back(B(i, j));
      }
    }
  }
  return result;
}

// EXTRACT RE and NIR
// [[Rcpp::export]]
std::vector<std::vector<double>> help_get_renir(const NumericMatrix &RE, const NumericMatrix &NIR, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = RE.nrow();
  int ncol = RE.ncol();
  // get the number of labels
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  // create a list to store the RGB values for each label
  std::vector<std::vector<double>> result(labelsCount);

  // loop through each label
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int label = labels(i, j);
      if (label > 0) {
        result[label].push_back(label);
        result[label].push_back(RE(i, j));
        result[label].push_back(NIR(i, j));
      }
    }
  }
  return result;
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

  // Get the unique labels in the mask, starting from 1
  IntegerVector unique_labels = sort_unique(labels);
  unique_labels = unique_labels[unique_labels != 0];

  // Create a list to store the isolated objects
  List isolated_objects(unique_labels.length());

  // Loop over the unique labels in the mask
  for (int i = 0; i < unique_labels.length(); i++) {
    int id = unique_labels[i];

    int top = nrows, bottom = 0, left = ncols, right = 0;

    // Loop over the rows and columns of the labels matrix
    for (int j = 0; j < nrows; j++) {
      for (int k = 0; k < ncols; k++) {
        if (labels(j,k) == id) {
          // Update the bounding box
          top = std::min(top, j);
          bottom = std::max(bottom, j);
          left = std::min(left, k);
          right = std::max(right, k);
        }
      }
    }

    // Crop the object
    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;
    NumericMatrix cropped(crop_nrows, crop_ncols);
    for (int j = 0; j < crop_nrows; j++) {
      for (int k = 0; k < crop_ncols; k++) {
        cropped(j,k) = img(top + j, left + k);
      }
    }

    // Add the isolated object to the list
    isolated_objects[i] = cropped;
  }

  return isolated_objects;
}




// HELPER FUNCTION TO ISOLATE OBJECTS BASED ON R-G-B and labels
// [[Rcpp::export]]
List help_isolate_object(NumericMatrix R, NumericMatrix G, NumericMatrix B, IntegerMatrix labels, bool remove_bg, int edge) {

  int nrows = labels.nrow(), ncols = labels.ncol();

  // Get the unique labels in the mask, starting from 1
  IntegerVector unique_labels = sort_unique(labels);
  unique_labels = unique_labels[unique_labels != 0];

  // Create a list to store the isolated objects
  List isolated_objects(unique_labels.length());

  // Loop over the unique labels in the mask
  for (int i = 0; i < unique_labels.length(); i++) {
    int id = unique_labels[i];

    int top = nrows, bottom = 0, left = ncols, right = 0;

    // Loop over the rows and columns of the labels matrix
    for (int j = 0; j < nrows; j++) {
      for (int k = 0; k < ncols; k++) {
        if (labels(j,k) == id) {
          // Update the bounding box
          top = std::min(top, j);
          bottom = std::max(bottom, j);
          left = std::min(left, k);
          right = std::max(right, k);

        }
      }
    }
    // Expand the bounding box by edge pixels
    top = std::max(0, top - edge);
    bottom = std::min(nrows - 1, bottom + edge);
    left = std::max(0, left - edge);
    right = std::min(ncols - 1, right + edge);

    // Crop the objects
    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;
    NumericMatrix croppedR(crop_nrows, crop_ncols);
    NumericMatrix croppedG(crop_nrows, crop_ncols);
    NumericMatrix croppedB(crop_nrows, crop_ncols);

    for (int j = 0; j < crop_nrows; j++) {
      for (int k = 0; k < crop_ncols; k++) {
        croppedR(j,k) = R(top + j, left + k);
        croppedG(j,k) = G(top + j, left + k);
        croppedB(j,k) = B(top + j, left + k);
      }
    }

    if(remove_bg){
      // Fill the pixels that are not part of the object with white
      for (int j = 0; j < crop_nrows; j++) {
        for (int k = 0; k < crop_ncols; k++) {
          if (labels(top + j, left + k) != id) {
            croppedR(j,k) = 1;
            croppedG(j,k) = 1;
            croppedB(j,k) = 1;
          }
        }
      }
    }

    // Store the isolated object in the list
    isolated_objects[i] = List::create(croppedR, croppedG, croppedB);
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





// Function to compute Otsu's threshold
// [[Rcpp::export]]
double help_otsu(const NumericVector& img) {
  int n = img.size();

  double x_max = max(img);
  double x_min = min(img);

  // Compute histogram
  std::vector<int> histogram(256, 0);
  for (int i = 0; i < n; i++) {
    int intensity = (int)(((1 - 0) / (x_max - x_min) * (img[i] - x_max) + 1) * 255);
    histogram[intensity]++;
  }

  // Compute total number of pixels
  int totalPixels = n;

  // Compute sum of intensities
  double sum = 0;
  for (int i = 0; i < 256; i++) {
    sum += i * histogram[i];
  }

  // Compute sum of background intensities
  double sumBackground = 0;
  int backgroundPixels = 0;

  // Initialize variables for storing optimal threshold and maximum between-class variance
  double maxVariance = 0;
  double threshold = 0;


  // Iterate through all possible thresholds
  for (int i = 0; i < 256; i++) {
    // Update background sum and number of background pixels
    backgroundPixels += histogram[i];
    sumBackground += i * histogram[i];

    // Calculate foreground and background weights
    double weightBackground = (double)backgroundPixels / totalPixels;
    double weightForeground = 1 - weightBackground;

    // Calculate mean intensities
    double meanBackground = sumBackground / backgroundPixels;
    double meanForeground = (sum - sumBackground) / (totalPixels - backgroundPixels);

    // Calculate between-class variance
    double variance = weightBackground * weightForeground * pow((meanBackground - meanForeground), 2);

    // Update maximum variance and threshold
    if (variance > maxVariance) {
      maxVariance = variance;
      threshold = i;
    }
  }

  // Scale the threshold value back to the range of 0-1

  return threshold * (x_max - x_min) / 255 + x_min;
}





// Function to apply Guo-Hall thinning algorithm to a binary image
// Adapted from https://observablehq.com/@esperanc/thinning#guoHall
// [[Rcpp::export]]
IntegerMatrix helper_guo_hall(IntegerMatrix image) {
  int wid = image.ncol();
  int hgt = image.nrow();
  IntegerMatrix data2 = Rcpp::clone(image);

  auto get = [&](int col, int row) { return image(row, col) != 0; };
  auto clear = [&](int col, int row) { data2(row, col) = 0; };

  IntegerMatrix stepCounter(wid, hgt);

  // Performs the conditional removal of one pixel. Even is true
  // if this is an even iteration.
  // Returns 1 if pixel was removed and 0 if not
  auto removePixel = [&](int col, int row, bool even) {
    if (!get(col, row)) return 0; // Not a 1-pixel
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

  // Performs one thinning step.
  // Returns the number of removed pixels
  auto thinStep = [&]() {
    int result = 0;
    for (int row = 1; row < hgt - 1; row++) {
      for (int col = 1; col < wid - 1; col++) {
        result += removePixel(col, row, even);
      }
    }
    even = !even;
    image = clone(data2); // Copy data2 back to image
    return result;
  };

  // Performs the thinning algorithm
  int n = 0;
  do {
    stepCounter.fill(0);
    n = thinStep();
  } while (n > 0);

  return image;
}

// IDW interpolation function using Rcpp
// [[Rcpp::export]]
NumericVector idw_interpolation_cpp(NumericVector x, NumericVector y, NumericVector values,
                                    NumericVector new_x, NumericVector new_y, double power = 2) {
  // Calculate distances between new points and existing points
  NumericMatrix distances(new_x.size(), x.size());
  for (int i = 0; i < new_x.size(); ++i) {
    distances(i, _) = sqrt(pow(x - new_x[i], 2) + pow(y - new_y[i], 2));
  }

  // Initialize a NumericVector to store results
  NumericVector results(new_x.size(), NA_REAL);

  for (int i = 0; i < new_x.size(); ++i) {
    // Inverse distance weighting formula
    NumericVector weights = 1.0 / pow(distances(i, _), power);
    double weighted_sum = sum(weights * values);
    double total_weight = sum(weights);
    results[i] = total_weight > 0 ? weighted_sum / total_weight : NA_REAL;

  }
  return results;
}
// Function to adjust the bounding box around the centroid
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
  new_bbox(4, 0) = xmin; new_bbox(4, 1) = ymin; // Closing the polygon

  return new_bbox;
}


// [[Rcpp::export]]
IntegerMatrix help_label(IntegerMatrix matrix, int max_gap = 2) {
  int rows = matrix.nrow();
  int cols = matrix.ncol();
  IntegerMatrix labels(rows, cols);
  int current_label = 0;

  // Função auxiliar para verificar se dois pixels estão dentro do gap permitido
  auto is_within_gap = [&](int r1, int c1, int r2, int c2) {
    return abs(r1 - r2) <= max_gap && abs(c1 - c2) <= max_gap;
  };

  // Pilhas para busca em profundidade
  std::vector<int> stack_r, stack_c;

  // Iterar sobre cada pixel da matriz
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      if (matrix(r, c) == 1 && labels(r, c) == 0) { // Novo objeto encontrado
        ++current_label;
        stack_r.push_back(r);
        stack_c.push_back(c);

        // Rotular os pixels conectados
        while (!stack_r.empty()) {
          int cr = stack_r.back();
          int cc = stack_c.back();
          stack_r.pop_back();
          stack_c.pop_back();

          if (labels(cr, cc) == 0) {
            labels(cr, cc) = current_label;

            // Verificar vizinhos
            for (int dr = -max_gap; dr <= max_gap; ++dr) {
              for (int dc = -max_gap; dc <= max_gap; ++dc) {
                if (abs(dr) + abs(dc) > 0) { // Ignorar o próprio pixel
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
    List geom = sf_coords[i]; // Get each geometry (may consist of multiple rings)
    double total_perimeter = 0.0;

    for (int j = 0; j < geom.size(); ++j) {
      NumericMatrix ring = geom[j]; // Single ring (matrix of coordinates)
      double ring_perimeter = 0.0;
      int rows = ring.nrow();

      for (int k = 0; k < rows - 1; ++k) {
        // Calculate Euclidean distance between consecutive points
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

// Helper function to generate random hexadecimal characters
std::string generate_random_hex(int length) {
  const char hex_chars[] = "0123456789abcdef";
  std::string result(length, '0');
  GetRNGstate(); // Inicia o RNG do R
  for (int i = 0; i < length; i++) {
    result[i] = hex_chars[(int)(unif_rand() * 16)]; // Gera um índice entre 0 e 15
  }
  PutRNGstate(); // Finaliza o RNG do R

  return result;
}

// [[Rcpp::export]]
std::string  uuid_v7() {
  // Step 1: Get current timestamp in milliseconds since Unix epoch
  auto now = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
  long long timestamp = duration.count();

  // Convert timestamp to a hexadecimal string with exactly 12 characters
  std::stringstream ss;
  ss << std::hex << std::setw(12) << std::setfill('0') << (timestamp & 0xFFFFFFFFFFFF); // Ensure 12 hex digits
  std::string timestamp_hex = ss.str();

  // Step 2: Extract components from the timestamp
  std::string time_low = timestamp_hex.substr(0, 8);  // First 8 hex digits (32 bits)
  std::string time_mid = timestamp_hex.substr(8, 4);  // Next 4 hex digits (16 bits)
  std::string time_high_and_version = generate_random_hex(4); // Generate 4 random hex chars
  time_high_and_version[0] = '7'; // Set the version to 7 (Version 7 UUID)

  // Step 3: Generate the variant and clock sequence
  std::string variant_and_sequence = generate_random_hex(4);
  GetRNGstate(); // Inicia o RNG do R para a variante
  variant_and_sequence[0] = "89ab"[(int)(unif_rand() * 4)]; // Define a variante (8, 9, a, b)
  PutRNGstate(); // Finaliza o RNG do R

  // Step 4: Generate the node (random bits for uniqueness)
  std::string node = generate_random_hex(12);

  // Step 5: Combine all components into the UUID structure
  std::string uuid = time_low + "-" + time_mid + "-" + time_high_and_version +
    "-" + variant_and_sequence + "-" + node;

  // Ensure UUID is exactly 36 characters long (with 4 hyphens)
  if (uuid.length() != 36) {
    throw std::runtime_error("Generated UUID has incorrect length: " + uuid);
  }

  return uuid;
}

// [[Rcpp::export]]
double helper_entropy(NumericVector values, int precision = 2) {
  std::unordered_map<double, int> freq;
  int n = values.size();

  // Compute frequencies with rounding
  double scale = pow(10.0, precision);
  for (int i = 0; i < n; ++i) {
    double rounded_val = round(values[i] * scale) / scale;
    freq[rounded_val]++;
  }

  // Compute entropy
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

    // build rows×2 matrix of coordinates
    NumericMatrix m(rows, 2);
    for (int i = 0; i < rows; ++i) {
      m(i, 0) = v[i];
      m(i, 1) = v[i + rows];
    }

    // take first four unique corners (drop closing point)
    NumericMatrix c4(4, 2);
    for (int i = 0; i < 4; ++i) {
      c4(i, 0) = m(i, 0);
      c4(i, 1) = m(i, 1);
    }

    // squared lengths of edges 1–2 and 2–3
    double dx1 = c4(0,0) - c4(1,0);
    double dy1 = c4(0,1) - c4(1,1);
    double d1  = dx1*dx1 + dy1*dy1;
    double dx2 = c4(1,0) - c4(2,0);
    double dy2 = c4(1,1) - c4(2,1);
    double d2  = dx2*dx2 + dy2*dy2;

    // choose longer-opposite edges midpoints
    double x1, y1, x2, y2;
    if (d1 < d2) {
      // edges 1–2 & 3–4 shorter => mids of those
      x1 = (c4(0,0) + c4(1,0)) * 0.5;
      y1 = (c4(0,1) + c4(1,1)) * 0.5;
      x2 = (c4(2,0) + c4(3,0)) * 0.5;
      y2 = (c4(2,1) + c4(3,1)) * 0.5;
    } else {
      // edges 2–3 & 4–1 shorter
      x1 = (c4(1,0) + c4(2,0)) * 0.5;
      y1 = (c4(1,1) + c4(2,1)) * 0.5;
      x2 = (c4(3,0) + c4(0,0)) * 0.5;
      y2 = (c4(3,1) + c4(0,1)) * 0.5;
    }

    // format WKT with fixed decimals
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
  arma::rowvec T_pixel(3); // A saída é sempre 3 (R,G,B)
  double R, G, B;

  // --- CAMINHO 1: Modelo "cubic" (9 termos) ---
  if (model == "cubic") {
    if (K.n_rows != 9) {
      Rcpp::stop("Erro: 'k_mat' tem %i linhas, mas o modelo 'cubic' espera 9.", K.n_rows);
    }

    arma::rowvec S_pixel(9); // Vetor de 9 termos

    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        R = img(i, j, 0);
        G = img(i, j, 1);
        B = img(i, j, 2);

        // Preenche os 9 termos
        S_pixel(0) = R;
        S_pixel(1) = G;
        S_pixel(2) = B;
        S_pixel(3) = R * R;
        S_pixel(4) = G * G;
        S_pixel(5) = B * B;
        S_pixel(6) = R * R * R;
        S_pixel(7) = G * G * G;
        S_pixel(8) = B * B * B;

        // (1x9) * (9x3) = (1x3)
        T_pixel = S_pixel * K;

        out_img(i, j, 0) = T_pixel(0);
        out_img(i, j, 1) = T_pixel(1);
        out_img(i, j, 2) = T_pixel(2);
      }
    }
  } else if (model == "root_polynomial") {
    if (K.n_rows != 20) {
      Rcpp::stop("Erro: 'k_mat' tem %i linhas, mas o modelo 'root_polynomial' espera 20.", K.n_rows);
    }

    arma::rowvec S_pixel(20); // Vetor de 20 termos
    double R2, G2, B2; // Variáveis intermediárias

    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        R = img(i, j, 0);
        G = img(i, j, 1);
        B = img(i, j, 2);

        R2 = R * R;
        G2 = G * G;
        B2 = B * B;

        // Preenche os 20 termos
        S_pixel(0) = 1.0; // Intercept
        S_pixel(1) = R;
        S_pixel(2) = G;
        S_pixel(3) = B;
        S_pixel(4) = R * G;
        S_pixel(5) = R * B;
        S_pixel(6) = G * B;
        S_pixel(7) = R2;
        S_pixel(8) = G2;
        S_pixel(9) = B2;
        S_pixel(10) = R2 * R; // R3
        S_pixel(11) = G2 * G; // G3
        S_pixel(12) = B2 * B; // B3
        S_pixel(13) = R2 * G;
        S_pixel(14) = R2 * B;
        S_pixel(15) = G2 * R;
        S_pixel(16) = G2 * B;
        S_pixel(17) = B2 * R;
        S_pixel(18) = B2 * G;
        S_pixel(19) = R * G * B;

        // (1x20) * (20x3) = (1x3)
        T_pixel = S_pixel * K;

        out_img(i, j, 0) = T_pixel(0);
        out_img(i, j, 1) = T_pixel(1);
        out_img(i, j, 2) = T_pixel(2);
      }
    }

    // --- CAMINHO 3: Erro ---
  } else {
    Rcpp::stop("Modelo '"+ model +"' não é reconhecido. Use 'cubic' ou 'root_polynomial'.");
  }

  return out_img;
}

double dist_eucl(double x1, double y1, double x2, double y2) {
  return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}
double dist_sq(double x1, double y1, double x2, double y2) {
  return std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2);
}
NumericVector get_point_at_dist(NumericMatrix coords, NumericVector cum_dist, double target_dist) {
  int n = coords.nrow();
  if (target_dist <= 0) return coords(0, _);
  if (target_dist >= cum_dist[n-1]) return coords(n-1, _);
  int i = 0;
  while(i < n - 1 && cum_dist[i+1] < target_dist) i++;
  double segment_len = cum_dist[i+1] - cum_dist[i];
  if (segment_len == 0) return coords(i, _);
  double t = (target_dist - cum_dist[i]) / segment_len;
  NumericVector out(2);
  out[0] = coords(i, 0) * (1 - t) + coords(i+1, 0) * t;
  out[1] = coords(i, 1) * (1 - t) + coords(i+1, 1) * t;
  return out;
}
int get_closest_idx_forward(NumericMatrix rail, double tx, double ty, int start_idx) {
  int n = rail.nrow();
  int best_idx = start_idx;
  double min_d2 = dist_sq(rail(start_idx, 0), rail(start_idx, 1), tx, ty);
  int search_limit = std::min(n, start_idx + 1000);

  for(int i = start_idx + 1; i < search_limit; i++) {
    double d2 = dist_sq(rail(i, 0), rail(i, 1), tx, ty);
    if (d2 < min_d2) {
      min_d2 = d2;
      best_idx = i;
    }
  }
  return best_idx;
}
// [[Rcpp::export]]
List make_grid_structure(NumericMatrix rail1,
                         NumericMatrix rail2,
                         int nrow,
                         int ncol,
                         double buffer_col,
                         double buffer_row,
                         Nullable<double> plot_width_opt,
                         Nullable<double> plot_height_opt) {
  int n1 = rail1.nrow();
  NumericVector cum_dist1(n1); cum_dist1[0] = 0;
  for(int i=1; i<n1; i++) cum_dist1[i] = cum_dist1[i-1] + dist_eucl(rail1(i-1,0), rail1(i-1,1), rail1(i,0), rail1(i,1));
  double total_len1 = cum_dist1[n1-1];

  int n2 = rail2.nrow();
  NumericVector cum_dist2(n2); cum_dist2[0] = 0;
  for(int i=1; i<n2; i++) cum_dist2[i] = cum_dist2[i-1] + dist_eucl(rail2(i-1,0), rail2(i-1,1), rail2(i,0), rail2(i,1));
  double total_len2 = cum_dist2[n2-1];

  double cell_along_avg = (total_len1 / ncol + total_len2 / ncol) / 2.0;
  double margin_along = 0.0;
  if (plot_width_opt.isNotNull()) {
    double pw = as<double>(plot_width_opt);
    margin_along = std::max(0.0, (cell_along_avg - pw) / 2.0);
  } else if (buffer_col > 0) {
    margin_along = buffer_col / 2.0;
  }

  List out_list(nrow * ncol);
  int idx = 0;

  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");

  for (int i = 0; i < ncol; i++) {
    double t_base_start = (double)i / ncol;
    double t_base_end   = (double)(i + 1) / ncol;
    double dist_start_l1 = t_base_start * total_len1 + margin_along;
    double dist_end_l1   = t_base_end * total_len1   - margin_along;
    double dist_start_l2 = t_base_start * total_len2 + margin_along;
    double dist_end_l2   = t_base_end * total_len2   - margin_along;

    NumericVector p_r1_start = get_point_at_dist(rail1, cum_dist1, dist_start_l1);
    NumericVector p_r1_end   = get_point_at_dist(rail1, cum_dist1, dist_end_l1);
    NumericVector p_r2_start = get_point_at_dist(rail2, cum_dist2, dist_start_l2);
    NumericVector p_r2_end   = get_point_at_dist(rail2, cum_dist2, dist_end_l2);

    double w_top = dist_eucl(p_r1_start[0], p_r1_start[1], p_r2_start[0], p_r2_start[1]);
    double w_bot = dist_eucl(p_r1_end[0],   p_r1_end[1],   p_r2_end[0],   p_r2_end[1]);
    double cell_cross_avg = (w_top + w_bot) / 2.0 / nrow;

    double margin_cross = 0.0;
    if (plot_height_opt.isNotNull()) {
      double ph = as<double>(plot_height_opt);
      margin_cross = std::max(0.0, (cell_cross_avg - ph) / 2.0);
    } else if (buffer_row > 0) {
      margin_cross = buffer_row / 2.0;
    }
    for (int j = 0; j < nrow; j++) {
      double u_base_start = (double)j / nrow;
      double u_base_end   = (double)(j + 1) / nrow;
      double u_margin_top = margin_cross / w_top;
      double u_margin_bot = margin_cross / w_bot;
      double u_start_top = u_base_start + u_margin_top;
      double u_end_top   = u_base_end   - u_margin_top;
      double u_start_bot = u_base_start + u_margin_bot;
      double u_end_bot   = u_base_end   - u_margin_bot;
      double x1 = p_r1_start[0] * (1-u_start_top) + p_r2_start[0] * u_start_top;
      double y1 = p_r1_start[1] * (1-u_start_top) + p_r2_start[1] * u_start_top;
      double x2 = p_r1_start[0] * (1-u_end_top)   + p_r2_start[0] * u_end_top;
      double y2 = p_r1_start[1] * (1-u_end_top)   + p_r2_start[1] * u_end_top;
      double x3 = p_r1_end[0] * (1-u_end_bot)   + p_r2_end[0] * u_end_bot;
      double y3 = p_r1_end[1] * (1-u_end_bot)   + p_r2_end[1] * u_end_bot;
      double x4 = p_r1_end[0] * (1-u_start_bot) + p_r2_end[0] * u_start_bot;
      double y4 = p_r1_end[1] * (1-u_start_bot) + p_r2_end[1] * u_start_bot;
      NumericMatrix ring(5, 2);
      ring(0,0) = x1; ring(0,1) = y1;
      ring(1,0) = x2; ring(1,1) = y2;
      ring(2,0) = x3; ring(2,1) = y3;
      ring(3,0) = x4; ring(3,1) = y4;
      ring(4,0) = x1; ring(4,1) = y1;
      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;

      out_list[idx] = polygon_sfg;
      idx++;
    }
  }
  return out_list;
}
// [[Rcpp::export]]
List make_grid_curved(NumericMatrix rail1,
                      NumericMatrix rail2,
                      int nrow,
                      int ncol,
                      bool curved = true,
                      int density = 20) {

  int n_dense = rail1.nrow();
  NumericMatrix centerline(n_dense, 2);
  NumericVector cl_cum_dist(n_dense);
  cl_cum_dist[0] = 0;

  centerline(0, _) = (rail1(0, _) + rail2(0, _)) / 2.0;

  for(int i=1; i<n_dense; i++) {
    centerline(i, 0) = (rail1(i, 0) + rail2(i, 0)) / 2.0;
    centerline(i, 1) = (rail1(i, 1) + rail2(i, 1)) / 2.0;
    cl_cum_dist[i] = cl_cum_dist[i-1] + dist_eucl(centerline(i-1,0), centerline(i-1,1), centerline(i,0), centerline(i,1));
  }
  double total_cl_len = cl_cum_dist[n_dense-1];

  List out_list(nrow * ncol);
  int idx = 0;
  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");

  int last_idx_r1 = 0;
  int last_idx_r2 = 0;
  int n_steps = curved ? density : 1;

  for (int i = 0; i < ncol; i++) {
    int idx_r1_s, idx_r2_s, idx_r1_e, idx_r2_e;
    if (i == 0) {
      idx_r1_s = 0; idx_r2_s = 0;
    } else {
      double dist_s = ((double)i / ncol) * total_cl_len;
      NumericVector p_cl_s = get_point_at_dist(centerline, cl_cum_dist, dist_s);
      idx_r1_s = get_closest_idx_forward(rail1, p_cl_s[0], p_cl_s[1], last_idx_r1);
      idx_r2_s = get_closest_idx_forward(rail2, p_cl_s[0], p_cl_s[1], last_idx_r2);
    }
    if (i == ncol - 1) {
      idx_r1_e = n_dense - 1; idx_r2_e = n_dense - 1;
    } else {
      double dist_e = ((double)(i + 1) / ncol) * total_cl_len;
      NumericVector p_cl_e = get_point_at_dist(centerline, cl_cum_dist, dist_e);
      idx_r1_e = get_closest_idx_forward(rail1, p_cl_e[0], p_cl_e[1], idx_r1_s);
      idx_r2_e = get_closest_idx_forward(rail2, p_cl_e[0], p_cl_e[1], idx_r2_s);
    }
    last_idx_r1 = idx_r1_s;
    last_idx_r2 = idx_r2_s;

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      NumericMatrix ring(2 * (n_steps + 1) + 1, 2);
      int pt_idx = 0;

      // Top Edge
      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;
        double curr_idx_r1 = idx_r1_s + t * (idx_r1_e - idx_r1_s);
        double curr_idx_r2 = idx_r2_s + t * (idx_r2_e - idx_r2_s);

        int i1 = (int)curr_idx_r1;
        int i2 = (int)curr_idx_r2;

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_top) + x_r2 * u_top;
        ring(pt_idx, 1) = y_r1 * (1 - u_top) + y_r2 * u_top;
        pt_idx++;
      }
      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        double curr_idx_r1 = idx_r1_s + t * (idx_r1_e - idx_r1_s);
        double curr_idx_r2 = idx_r2_s + t * (idx_r2_e - idx_r2_s);

        int i1 = (int)curr_idx_r1;
        int i2 = (int)curr_idx_r2;

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_bot) + x_r2 * u_bot;
        ring(pt_idx, 1) = y_r1 * (1 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }
      ring(pt_idx, 0) = ring(0, 0);
      ring(pt_idx, 1) = ring(0, 1);

      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;
      out_list[idx++] = polygon_sfg;
    }
  }
  return out_list;
}

// [[Rcpp::export]]
List make_grid_landmarks(NumericMatrix rail1,
                         NumericMatrix rail2,
                         IntegerVector anchors1,
                         IntegerVector anchors2,
                         int nrow,
                         bool curved = true,
                         int density = 30) {

  int n_cols = anchors1.size() - 1;
  if (anchors2.size() != anchors1.size()) {
    stop("Rail 1 and Rail 2 must have the same number of control points for manual mode.");
  }

  List out_list(nrow * n_cols);
  int idx = 0;
  CharacterVector sfg_class = CharacterVector::create("XY", "POLYGON", "sfg");
  int n_steps = curved ? density : 1;

  for (int i = 0; i < n_cols; i++) {
    int idx_r1_start = anchors1[i];
    int idx_r1_end   = anchors1[i+1];

    int idx_r2_start = anchors2[i];
    int idx_r2_end   = anchors2[i+1];

    double diff_r1 = (double)(idx_r1_end - idx_r1_start);
    double diff_r2 = (double)(idx_r2_end - idx_r2_start);

    for (int j = 0; j < nrow; j++) {
      double u_top = (double)j / nrow;
      double u_bot = (double)(j + 1) / nrow;

      NumericMatrix ring(2 * (n_steps + 1) + 1, 2);
      int pt_idx = 0;
      for (int k = 0; k <= n_steps; k++) {
        double t = (double)k / n_steps;

        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);

        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }

        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);

        ring(pt_idx, 0) = x_r1 * (1 - u_top) + x_r2 * u_top;
        ring(pt_idx, 1) = y_r1 * (1 - u_top) + y_r2 * u_top;
        pt_idx++;
      }

      for (int k = n_steps; k >= 0; k--) {
        double t = (double)k / n_steps;
        int i1 = idx_r1_start + (int)(t * diff_r1);
        int i2 = idx_r2_start + (int)(t * diff_r2);
        if (k == n_steps) { i1 = idx_r1_end; i2 = idx_r2_end; }
        double x_r1 = rail1(i1, 0); double y_r1 = rail1(i1, 1);
        double x_r2 = rail2(i2, 0); double y_r2 = rail2(i2, 1);
        ring(pt_idx, 0) = x_r1 * (1 - u_bot) + x_r2 * u_bot;
        ring(pt_idx, 1) = y_r1 * (1 - u_bot) + y_r2 * u_bot;
        pt_idx++;
      }
      ring(pt_idx, 0) = ring(0, 0);
      ring(pt_idx, 1) = ring(0, 1);

      List polygon_sfg(1);
      polygon_sfg[0] = ring;
      polygon_sfg.attr("class") = sfg_class;
      out_list[idx++] = polygon_sfg;
    }
  }
  return out_list;
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
