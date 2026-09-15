#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include "lsd.h"
}

// [[Rcpp::export]]
List detect_line_segments(SEXP image_sexp, int X, int Y,
                          double scale = 0.8,
                          double sigma_scale = 0.6,
                          double quant = 2.0,
                          double ang_th = 22.5,
                          double log_eps = 0.0,
                          double density_th = 0.7,
                          int n_bins = 1024,
                          int need_to_union = 0,
                          double union_ang_th = 7,
                          int union_use_NFA = 0,
                          double union_log_eps = 0.0,
                          double length_threshold = 5,
                          double dist_threshold = 5)
{
  int point_num = Rf_length(image_sexp);
  if (point_num != X * Y) {
    stop("Size of image not the same as X*Y");
  }

  std::vector<double> img_buffer;
  double* img_ptr = nullptr;

  if (TYPEOF(image_sexp) == RAWSXP) {
    img_buffer.resize(point_num);
    Rbyte* rp = RAW(image_sexp);
    for(int i = 0; i < point_num; i++) img_buffer[i] = rp[i] / 255.0;
    img_ptr = img_buffer.data();
  } else if (TYPEOF(image_sexp) == REALSXP) {
    img_ptr = REAL(image_sexp);
  } else if (TYPEOF(image_sexp) == INTSXP) {
    img_buffer.resize(point_num);
    int* ip = INTEGER(image_sexp);
    for(int i = 0; i < point_num; i++) img_buffer[i] = ip[i] / 255.0;
    img_ptr = img_buffer.data();
  } else {
    stop("Unsupported image type");
  }
  int n_out, reg_x, reg_y;
  double * out;
  int * reg_img;
  out = LineSegmentDetection( &n_out, img_ptr, X, Y, scale, sigma_scale, quant,
                              ang_th, log_eps, density_th, union_ang_th,
                              union_use_NFA, union_log_eps, n_bins,
                              need_to_union, &reg_img, &reg_x, &reg_y,
                              length_threshold, dist_threshold );

  NumericMatrix y(n_out, 7);
  for(int i = 0; i < n_out; i++)
  {
    // x1, y1, x2, y2, width, p, -log_nfa
    y(i, 0) = out[7*i+1]; 
    y(i, 1) = X - out[7*i+0]; 
    y(i, 2) = out[7*i+3]; 
    y(i, 3) = X - out[7*i+2]; 
    y(i, 4) = out[7*i+4]; 
    y(i, 5) = out[7*i+5];
    y(i, 6) = out[7*i+6]; 
  }
  NumericMatrix p(reg_x, reg_y);
  for(int i = 0; i < reg_x * reg_y; i++) p[i] = reg_img[i];
  List z = List::create(y, p);
  return z;
}



