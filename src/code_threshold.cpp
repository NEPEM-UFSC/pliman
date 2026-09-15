#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// adapted from imagerExtra https://bit.ly/3HtxumB
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
      double sd_local = sqrt((int_sum_sqr(nrow-1,j+winhalf) + int_sum_sqr(i-winhalf,j-winhalf) - integ_sum(nrow-1,j-winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
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

// Function to compute Otsu's threshold
// [[Rcpp::export]]
double help_otsu(SEXP img_sexp) {
  int npix = Rf_length(img_sexp);
  if (npix == 0) return 0.0;

  std::vector<int> histogram(256, 0);
  double x_min = 0.0, x_max = 1.0;
  int valid_pixels = 0;

  if (TYPEOF(img_sexp) == RAWSXP) {
    const uint8_t* ptr = RAW(img_sexp);
    for (int i = 0; i < npix; i++) {
      histogram[ptr[i]]++;
    }
    x_min = 0.0;
    x_max = 255.0;
    valid_pixels = npix;
  } else if (TYPEOF(img_sexp) == REALSXP) {
    const double* ptr = REAL(img_sexp);
    bool first = true;
    for (int i = 0; i < npix; i++) {
      double v = ptr[i];
      if (!std::isnan(v)) {
        if (first) {
          x_min = v;
          x_max = v;
          first = false;
        } else {
          if (v < x_min) x_min = v;
          if (v > x_max) x_max = v;
        }
        valid_pixels++;
      }
    }
    if (valid_pixels == 0) return 0.0;
    double range = x_max - x_min;
    if (range <= 0.0) return x_min;

    for (int i = 0; i < npix; i++) {
      double v = ptr[i];
      if (!std::isnan(v)) {
        int idx = (int)std::round(((v - x_min) / range) * 255.0);
        if (idx < 0) idx = 0;
        if (idx > 255) idx = 255;
        histogram[idx]++;
      }
    }
  } else {
    Rcpp::NumericVector img(img_sexp);
    bool first = true;
    for (int i = 0; i < npix; i++) {
      double v = img[i];
      if (!Rcpp::NumericVector::is_na(v) && !std::isnan(v)) {
        if (first) {
          x_min = v;
          x_max = v;
          first = false;
        } else {
          if (v < x_min) x_min = v;
          if (v > x_max) x_max = v;
        }
        valid_pixels++;
      }
    }
    if (valid_pixels == 0) return 0.0;
    double range = x_max - x_min;
    if (range <= 0.0) return x_min;

    for (int i = 0; i < npix; i++) {
      double v = img[i];
      if (!Rcpp::NumericVector::is_na(v) && !std::isnan(v)) {
        int idx = (int)std::round(((v - x_min) / range) * 255.0);
        if (idx < 0) idx = 0;
        if (idx > 255) idx = 255;
        histogram[idx]++;
      }
    }
  }

  int totalPixels = valid_pixels;
  if (totalPixels == 0) return 0.0;

  double sum = 0;
  for (int i = 0; i < 256; i++) {
    sum += (double)i * histogram[i];
  }

  double sumBackground = 0;
  int backgroundPixels = 0;
  double maxVariance = 0;
  double threshold = 0;

  for (int i = 0; i < 256; i++) {
    backgroundPixels += histogram[i];
    if (backgroundPixels == 0) continue;
    int foregroundPixels = totalPixels - backgroundPixels;
    if (foregroundPixels == 0) break;

    sumBackground += (double)i * histogram[i];

    double meanBackground = sumBackground / backgroundPixels;
    double meanForeground = (sum - sumBackground) / foregroundPixels;

    double variance = (double)backgroundPixels * (double)foregroundPixels * (meanBackground - meanForeground) * (meanBackground - meanForeground);

    if (variance > maxVariance) {
      maxVariance = variance;
      threshold = i;
    }
  }

  if (TYPEOF(img_sexp) == RAWSXP) {
    return threshold;
  } else {
    return threshold * (x_max - x_min) / 255.0 + x_min;
  }
}

// [[Rcpp::export]]
SEXP cpp_binary_threshold(SEXP img_sexp, double threshold, int op = 1, bool return_raw = true) {
  SEXP dims_sexp = Rf_getAttrib(img_sexp, R_DimSymbol);
  int npix = Rf_length(img_sexp);
  if (npix == 0) return R_NilValue;

  int thresh_raw = (threshold <= 1.0 && threshold >= 0) ? (int)std::round(threshold * 255.0) : (int)std::round(threshold);
  thresh_raw = std::max(0, std::min(255, thresh_raw));

  int r_h = (!Rf_isNull(dims_sexp) && Rf_length(dims_sexp) >= 1) ? INTEGER(dims_sexp)[0] : npix;
  int r_w = (!Rf_isNull(dims_sexp) && Rf_length(dims_sexp) >= 2) ? INTEGER(dims_sexp)[1] : 1;

  if (return_raw) {
    SEXP res = PROTECT(Rf_allocMatrix(RAWSXP, r_h, r_w));
    if (!Rf_isNull(dims_sexp)) {
      Rf_setAttrib(res, R_DimSymbol, dims_sexp);
    }
    uint8_t* out = RAW(res);

    if (TYPEOF(img_sexp) == RAWSXP) {
      const uint8_t* ptr = RAW(img_sexp);
      const uint8_t tr = (uint8_t)thresh_raw;
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < tr) ? 255 : 0;
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= tr) ? 255 : 0;
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > tr) ? 255 : 0;
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= tr) ? 255 : 0;
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == tr) ? 255 : 0;
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != tr) ? 255 : 0;
          break;
      }
    } else if (TYPEOF(img_sexp) == REALSXP) {
      const double* ptr = REAL(img_sexp);
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < threshold) ? 255 : 0;
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= threshold) ? 255 : 0;
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > threshold) ? 255 : 0;
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= threshold) ? 255 : 0;
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == threshold) ? 255 : 0;
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != threshold) ? 255 : 0;
          break;
      }
    } else if (TYPEOF(img_sexp) == INTSXP) {
      const int* ptr = INTEGER(img_sexp);
      int tr_i = (int)threshold;
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < tr_i) ? 255 : 0;
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= tr_i) ? 255 : 0;
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > tr_i) ? 255 : 0;
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= tr_i) ? 255 : 0;
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == tr_i) ? 255 : 0;
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != tr_i) ? 255 : 0;
          break;
      }
    }
    SEXP cls = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(cls, 0, Rf_mkChar("image"));
    SET_STRING_ELT(cls, 1, Rf_mkChar("array"));
    Rf_setAttrib(res, R_ClassSymbol, cls);
    Rf_setAttrib(res, Rf_install("colormode"), Rf_mkString("Grayscale"));
    UNPROTECT(1);

    UNPROTECT(1);
    return res;
  } else {
    SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, r_h, r_w));
    if (!Rf_isNull(dims_sexp)) {
      Rf_setAttrib(res, R_DimSymbol, dims_sexp);
    }
    int* out = LOGICAL(res);

    if (TYPEOF(img_sexp) == RAWSXP) {
      const uint8_t* ptr = RAW(img_sexp);
      const uint8_t tr = (uint8_t)thresh_raw;
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < tr);
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= tr);
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > tr);
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= tr);
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == tr);
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != tr);
          break;
      }
    } else if (TYPEOF(img_sexp) == REALSXP) {
      const double* ptr = REAL(img_sexp);
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < threshold);
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= threshold);
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > threshold);
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= threshold);
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == threshold);
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != threshold);
          break;
      }
    } else if (TYPEOF(img_sexp) == INTSXP) {
      const int* ptr = INTEGER(img_sexp);
      int tr_i = (int)threshold;
      switch (op) {
        case 1:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] < tr_i);
          break;
        case 2:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] <= tr_i);
          break;
        case 3:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] > tr_i);
          break;
        case 4:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] >= tr_i);
          break;
        case 5:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] == tr_i);
          break;
        case 6:
          #pragma omp parallel for if(npix > 250000)
          for (int i = 0; i < npix; i++) out[i] = (ptr[i] != tr_i);
          break;
      }
    }
    SEXP cls = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(cls, 0, Rf_mkChar("image"));
    SET_STRING_ELT(cls, 1, Rf_mkChar("array"));
    Rf_setAttrib(res, R_ClassSymbol, cls);
    Rf_setAttrib(res, Rf_install("colormode"), Rf_mkString("Grayscale"));
    UNPROTECT(1);

    UNPROTECT(1);
    return res;
  }
}
