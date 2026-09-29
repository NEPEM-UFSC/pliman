#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix haralick_features_cpp(IntegerMatrix labels, SEXP ref_sexp, int nc = 32) {
  int nrow = labels.nrow(), ncol = labels.ncol();
  const int* p_labels = labels.begin();
  
  // Find maximum label
  int max_label = 0;
  int total_cells = labels.size();
  for (int i = 0; i < total_cells; i++) {
    if (p_labels[i] > max_label) max_label = p_labels[i];
  }
  
  // Create output matrix: max_label rows, 13 columns
  NumericMatrix out(max_label, 13);
  colnames(out) = CharacterVector::create(
    "h.asm", "h.con", "h.cor", "h.var", "h.idm",
    "h.sav", "h.sva", "h.sen", "h.ent", "h.dva", "h.den", "h.f12", "h.f13"
  );
  
  if (max_label == 0) return out;
  
  // Group coordinates by label for fast O(N) lookup.
  // Pack r and c into a single 32-bit integer: (r << 16) | c
  std::vector<std::vector<int>> label_pixels(max_label + 1);
  for (int c = 0; c < ncol; c++) {
    int c_offset = c * nrow;
    for (int r = 0; r < nrow; r++) {
      int lbl = p_labels[c_offset + r];
      if (lbl > 0 && lbl <= max_label) {
        label_pixels[lbl].push_back((r << 16) | c);
      }
    }
  }
  
  // Pre-quantize reference image to range [0, nc - 1]
  std::vector<int> quant(total_cells, 0);
  if (TYPEOF(ref_sexp) == RAWSXP) {
    Rbyte* p_ref = RAW(ref_sexp);
    for (int i = 0; i < total_cells; i++) {
      double val = p_ref[i] / 255.0;
      if (val > 1.0) val = 1.0;
      if (val < 0.0) val = 0.0;
      int q = static_cast<int>(std::floor(val * (nc - 1)));
      if (q < 0) q = 0;
      if (q >= nc) q = nc - 1;
      quant[i] = q;
    }
  } else if (TYPEOF(ref_sexp) == REALSXP) {
    double* p_ref = REAL(ref_sexp);
    double max_val = 0.0;
    for (int i = 0; i < total_cells; i++) {
      if (p_ref[i] > max_val) max_val = p_ref[i];
    }
    double denom = (max_val > 1.0) ? 255.0 : 1.0;
    for (int i = 0; i < total_cells; i++) {
      double val = p_ref[i] / denom;
      if (val > 1.0) val = 1.0;
      if (val < 0.0) val = 0.0;
      int q = static_cast<int>(std::floor(val * (nc - 1)));
      if (q < 0) q = 0;
      if (q >= nc) q = nc - 1;
      quant[i] = q;
    }
  } else if (TYPEOF(ref_sexp) == INTSXP) {
    int* p_ref = INTEGER(ref_sexp);
    for (int i = 0; i < total_cells; i++) {
      double val = p_ref[i] / 255.0;
      if (val > 1.0) val = 1.0;
      if (val < 0.0) val = 0.0;
      int q = static_cast<int>(std::floor(val * (nc - 1)));
      if (q < 0) q = 0;
      if (q >= nc) q = nc - 1;
      quant[i] = q;
    }
  } else {
    Rcpp::stop("haralick_features_cpp: unsupported ref type.");
  }
  
  // Helper log10 function
  auto log10_safe = [](double x) {
    return (x > 1e-15) ? std::log10(x) : 0.0;
  };
  
  // For each label, compute GLCM and Haralick features
  for (int lbl = 1; lbl <= max_label; lbl++) {
    const auto& pixels = label_pixels[lbl];
    if (pixels.empty()) continue;
    
    // Allocate GLCM matrix (symmetric)
    std::vector<double> glcm(nc * nc, 0.0);
    double total_pairs = 0.0;
    
    // 4 directions: (0, 1), (1, 0), (1, 1), (1, -1)
    int drs[4] = {0, 1, 1, 1};
    int dcs[4] = {1, 0, 1, -1};
    
    for (int val : pixels) {
      int r = val >> 16;
      int c = val & 0xFFFF;
      int idx = c * nrow + r;
      int val_i = quant[idx];
      
      for (int d = 0; d < 4; d++) {
        int nr = r + drs[d];
        int nc_val = c + dcs[d];
        
        if (static_cast<unsigned>(nr) < static_cast<unsigned>(nrow) &&
            static_cast<unsigned>(nc_val) < static_cast<unsigned>(ncol)) {
          
          int n_idx = nc_val * nrow + nr;
          if (p_labels[n_idx] == lbl) {
            int val_j = quant[n_idx];
            glcm[val_i * nc + val_j] += 1.0;
            glcm[val_j * nc + val_i] += 1.0;
            total_pairs += 2.0;
          }
        }
      }
    }
    
    if (total_pairs == 0.0) continue;
    
    // Normalize GLCM
    for (int i = 0; i < nc * nc; i++) {
      glcm[i] /= total_pairs;
    }
    
    // Compute Haralick features
    double asm_val = 0.0;
    double con = 0.0;
    double idm = 0.0;
    double ent = 0.0;
    
    std::vector<double> px(nc, 0.0);
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nc; j++) {
        double p = glcm[i * nc + j];
        asm_val += p * p;
        con += (i - j) * (i - j) * p;
        idm += p / (1.0 + (i - j) * (i - j));
        ent -= p * log10_safe(p);
        px[i] += p;
      }
    }
    
    // Mean and variance
    double mu = 0.0;
    for (int i = 0; i < nc; i++) {
      mu += i * px[i];
    }
    double var = 0.0;
    for (int i = 0; i < nc; i++) {
      var += (i + 1.0 - mu) * (i + 1.0 - mu) * px[i];
    }
    
    double cor = 0.0;
    double var_0 = var - 1.0;
    if (var_0 > 1e-10) {
      for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
          double p = glcm[i * nc + j];
          cor += (i - mu) * (j - mu) * p;
        }
      }
      cor /= var_0;
    } else {
      cor = 1.0;
    }
    
    // Sum average, Sum variance, Sum entropy
    int num_sums = 2 * nc - 1;
    std::vector<double> p_xplusy(num_sums, 0.0);
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nc; j++) {
        p_xplusy[i + j] += glcm[i * nc + j];
      }
    }
    
    double sav = 0.0;
    double sen = 0.0;
    for (int k = 0; k < num_sums; k++) {
      sav += (k + 2.0) * p_xplusy[k];
      double p = p_xplusy[k];
      sen -= p * log10_safe(p);
    }
    
    double sva = 0.0;
    for (int k = 0; k < num_sums; k++) {
      sva += (k + 2.0 - sen) * (k + 2.0 - sen) * p_xplusy[k];
    }
    
    // Difference variance, Difference entropy
    std::vector<double> p_xminusy(nc, 0.0);
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nc; j++) {
        p_xminusy[std::abs(i - j)] += glcm[i * nc + j];
      }
    }
    
    double den = 0.0;
    for (int k = 0; k < nc; k++) {
      double p = p_xminusy[k];
      den -= p * log10_safe(p);
    }
    
    double dva = con;
    
    double hx = 0.0;
    for (int i = 0; i < nc; i++) {
      hx -= px[i] * log10_safe(px[i]);
    }
    
    double hxy1 = 0.0;
    double hxy2 = 0.0;
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nc; j++) {
        double p_x = px[i];
        double p_y = px[j];
        if (p_x * p_y > 1e-15) {
          hxy1 -= glcm[i * nc + j] * log10_safe(p_x * p_y);
          hxy2 -= p_x * p_y * log10_safe(p_x * p_y);
        }
      }
    }
    
    double f12 = 0.0;
    if (hx > 1e-10) {
      f12 = (hxy1 - ent) / hx;
    }
    
    double f13 = 0.0;
    double term = hxy2 - ent;
    if (term > 0.0) {
      f13 = std::sqrt(1.0 - std::exp(-2.0 * term));
    }
    
    out(lbl - 1, 0) = asm_val;
    out(lbl - 1, 1) = con;
    out(lbl - 1, 2) = cor;
    out(lbl - 1, 3) = var;
    out(lbl - 1, 4) = idm;
    out(lbl - 1, 5) = sav;
    out(lbl - 1, 6) = sva;
    out(lbl - 1, 7) = sen;
    out(lbl - 1, 8) = ent;
    out(lbl - 1, 9) = dva;
    out(lbl - 1, 10) = den;
    out(lbl - 1, 11) = f12;
    out(lbl - 1, 12) = f13;
  }
  
  return out;
}
