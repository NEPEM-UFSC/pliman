#include <RcppArmadillo.h>
#include <stack>
#include <cmath>
#include <cctype>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

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

// EXTRACT PIXELS
std::vector<std::vector<double>> help_get_rgb(const NumericMatrix &R, const NumericMatrix &G, const NumericMatrix &B, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = R.nrow();
  int ncol = R.ncol();
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  std::vector<std::vector<double>> result(labelsCount);
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
std::vector<std::vector<double>> help_get_renir(const NumericMatrix &RE, const NumericMatrix &NIR, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = RE.nrow();
  int ncol = RE.ncol();
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  std::vector<std::vector<double>> result(labelsCount);
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

enum TokenType {
  TOK_NUM,
  TOK_VAR_R, TOK_VAR_G, TOK_VAR_B, TOK_VAR_RE, TOK_VAR_NIR, TOK_VAR_SWIR,
  TOK_ADD, TOK_SUB, TOK_MUL, TOK_DIV, TOK_POW, TOK_NEG,
  TOK_ABS, TOK_SQRT, TOK_EXP, TOK_LOG, TOK_SIN, TOK_COS, TOK_TAN
};

struct RPNToken {
  TokenType type;
  double val;
};

inline std::vector<RPNToken> parse_to_rpn(const std::string& expr_str) {
  std::vector<RPNToken> rpn;
  std::stack<std::string> op_stack;

  auto precedence = [](const std::string& op) -> int {
    if (op == "neg") return 4;
    if (op == "^") return 3;
    if (op == "*" || op == "/") return 2;
    if (op == "+" || op == "-") return 1;
    return 0;
  };

  auto is_right_assoc = [](const std::string& op) -> bool {
    return (op == "^" || op == "neg");
  };

  auto is_func = [](const std::string& s) -> bool {
    return (s == "abs" || s == "sqrt" || s == "exp" || s == "log" || s == "sin" || s == "cos" || s == "tan");
  };

  size_t n = expr_str.size();
  size_t i = 0;
  bool expect_unary = true;

  while (i < n) {
    if (std::isspace(expr_str[i])) {
      i++;
      continue;
    }

    if (std::isdigit(expr_str[i]) || (expr_str[i] == '.' && i + 1 < n && std::isdigit(expr_str[i + 1]))) {
      size_t len = 0;
      double val = std::stod(expr_str.substr(i), &len);
      i += len;
      RPNToken tok;
      tok.type = TOK_NUM;
      tok.val = val;
      rpn.push_back(tok);
      expect_unary = false;
      continue;
    }

    if (std::isalpha(expr_str[i])) {
      std::string id = "";
      while (i < n && (std::isalnum(expr_str[i]) || expr_str[i] == '_')) {
        id += (char)std::toupper(expr_str[i]);
        i++;
      }

      std::string id_lower = id;
      std::transform(id_lower.begin(), id_lower.end(), id_lower.begin(), ::tolower);

      if (is_func(id_lower)) {
        op_stack.push(id_lower);
        expect_unary = true;
      } else {
        RPNToken tok;
        if (id == "R") tok.type = TOK_VAR_R;
        else if (id == "G") tok.type = TOK_VAR_G;
        else if (id == "B") tok.type = TOK_VAR_B;
        else if (id == "RE") tok.type = TOK_VAR_RE;
        else if (id == "NIR") tok.type = TOK_VAR_NIR;
        else if (id == "SWIR") tok.type = TOK_VAR_SWIR;
        else {
          throw std::runtime_error("Unknown variable: " + id);
        }
        rpn.push_back(tok);
        expect_unary = false;
      }
      continue;
    }

    char c = expr_str[i];
    if (c == '(') {
      op_stack.push("(");
      i++;
      expect_unary = true;
    } else if (c == ')') {
      while (!op_stack.empty() && op_stack.top() != "(") {
        std::string top_op = op_stack.top();
        op_stack.pop();
        RPNToken tok;
        if (top_op == "+") tok.type = TOK_ADD;
        else if (top_op == "-") tok.type = TOK_SUB;
        else if (top_op == "*") tok.type = TOK_MUL;
        else if (top_op == "/") tok.type = TOK_DIV;
        else if (top_op == "^") tok.type = TOK_POW;
        else if (top_op == "neg") tok.type = TOK_NEG;
        else if (top_op == "abs") tok.type = TOK_ABS;
        else if (top_op == "sqrt") tok.type = TOK_SQRT;
        else if (top_op == "exp") tok.type = TOK_EXP;
        else if (top_op == "log") tok.type = TOK_LOG;
        else if (top_op == "sin") tok.type = TOK_SIN;
        else if (top_op == "cos") tok.type = TOK_COS;
        else if (top_op == "tan") tok.type = TOK_TAN;
        rpn.push_back(tok);
      }
      if (!op_stack.empty() && op_stack.top() == "(") {
        op_stack.pop();
      }
      if (!op_stack.empty() && is_func(op_stack.top())) {
        std::string fn = op_stack.top();
        op_stack.pop();
        RPNToken tok;
        if (fn == "abs") tok.type = TOK_ABS;
        else if (fn == "sqrt") tok.type = TOK_SQRT;
        else if (fn == "exp") tok.type = TOK_EXP;
        else if (fn == "log") tok.type = TOK_LOG;
        else if (fn == "sin") tok.type = TOK_SIN;
        else if (fn == "cos") tok.type = TOK_COS;
        else if (fn == "tan") tok.type = TOK_TAN;
        rpn.push_back(tok);
      }
      i++;
      expect_unary = false;
    } else if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
      std::string op(1, c);
      if (c == '-' && expect_unary) {
        op = "neg";
      }

      int p_curr = precedence(op);
      while (!op_stack.empty() && op_stack.top() != "(") {
        std::string top_op = op_stack.top();
        int p_top = precedence(top_op);
        if (p_top > p_curr || (p_top == p_curr && !is_right_assoc(op))) {
          op_stack.pop();
          RPNToken tok;
          if (top_op == "+") tok.type = TOK_ADD;
          else if (top_op == "-") tok.type = TOK_SUB;
          else if (top_op == "*") tok.type = TOK_MUL;
          else if (top_op == "/") tok.type = TOK_DIV;
          else if (top_op == "^") tok.type = TOK_POW;
          else if (top_op == "neg") tok.type = TOK_NEG;
          else if (top_op == "abs") tok.type = TOK_ABS;
          else if (top_op == "sqrt") tok.type = TOK_SQRT;
          else if (top_op == "exp") tok.type = TOK_EXP;
          else if (top_op == "log") tok.type = TOK_LOG;
          else if (top_op == "sin") tok.type = TOK_SIN;
          else if (top_op == "cos") tok.type = TOK_COS;
          else if (top_op == "tan") tok.type = TOK_TAN;
          rpn.push_back(tok);
        } else {
          break;
        }
      }
      op_stack.push(op);
      i++;
      expect_unary = true;
    } else {
      i++;
    }
  }

  while (!op_stack.empty()) {
    std::string top_op = op_stack.top();
    op_stack.pop();
    if (top_op != "(" && top_op != ")") {
      RPNToken tok;
      if (top_op == "+") tok.type = TOK_ADD;
      else if (top_op == "-") tok.type = TOK_SUB;
      else if (top_op == "*") tok.type = TOK_MUL;
      else if (top_op == "/") tok.type = TOK_DIV;
      else if (top_op == "^") tok.type = TOK_POW;
      else if (top_op == "neg") tok.type = TOK_NEG;
      else if (top_op == "abs") tok.type = TOK_ABS;
      else if (top_op == "sqrt") tok.type = TOK_SQRT;
      else if (top_op == "exp") tok.type = TOK_EXP;
      else if (top_op == "log") tok.type = TOK_LOG;
      else if (top_op == "sin") tok.type = TOK_SIN;
      else if (top_op == "cos") tok.type = TOK_COS;
      else if (top_op == "tan") tok.type = TOK_TAN;
      rpn.push_back(tok);
    }
  }

  return rpn;
}

inline double eval_rpn(const std::vector<RPNToken>& rpn, double r, double g, double b, double re, double nir, double swir) {
  double st[32];
  int top = 0;

  for (size_t k = 0; k < rpn.size(); ++k) {
    const RPNToken& tok = rpn[k];
    switch (tok.type) {
      case TOK_NUM: st[top++] = tok.val; break;
      case TOK_VAR_R: st[top++] = r; break;
      case TOK_VAR_G: st[top++] = g; break;
      case TOK_VAR_B: st[top++] = b; break;
      case TOK_VAR_RE: st[top++] = re; break;
      case TOK_VAR_NIR: st[top++] = nir; break;
      case TOK_VAR_SWIR: st[top++] = swir; break;
      case TOK_ADD: { double val2 = st[--top]; double val1 = st[--top]; st[top++] = val1 + val2; break; }
      case TOK_SUB: { double val2 = st[--top]; double val1 = st[--top]; st[top++] = val1 - val2; break; }
      case TOK_MUL: { double val2 = st[--top]; double val1 = st[--top]; st[top++] = val1 * val2; break; }
      case TOK_DIV: { double val2 = st[--top]; double val1 = st[--top]; st[top++] = (val2 != 0.0) ? (val1 / val2) : 0.0; break; }
      case TOK_POW: { double val2 = st[--top]; double val1 = st[--top]; st[top++] = std::pow(val1, val2); break; }
      case TOK_NEG: { st[top - 1] = -st[top - 1]; break; }
      case TOK_ABS: { st[top - 1] = std::abs(st[top - 1]); break; }
      case TOK_SQRT: { double v = st[top - 1]; st[top - 1] = (v >= 0.0) ? std::sqrt(v) : 0.0; break; }
      case TOK_EXP: { st[top - 1] = std::exp(st[top - 1]); break; }
      case TOK_LOG: { double v = st[top - 1]; st[top - 1] = (v > 0.0) ? std::log(v) : 0.0; break; }
      case TOK_SIN: { st[top - 1] = std::sin(st[top - 1]); break; }
      case TOK_COS: { st[top - 1] = std::cos(st[top - 1]); break; }
      case TOK_TAN: { st[top - 1] = std::tan(st[top - 1]); break; }
    }
  }
  return (top > 0) ? st[0] : 0.0;
}

enum IndexOp {
  OP_R, OP_G, OP_B, OP_RE, OP_NIR, OP_SWIR, OP_NR, OP_NG, OP_NB, OP_GRAY, OP_GRAY2,
  OP_IPCA, OP_GLI, OP_VARI, OP_NGBDI, OP_BGI, OP_BI, OP_CI, OP_CIVE, OP_EGVI, OP_ERVI,
  OP_GB, OP_GD, OP_GLAI, OP_GR, OP_HI, OP_HUE, OP_HUE2, OP_I, OP_L, OP_MGVRI, OP_RB,
  OP_RI, OP_S, OP_SAVI, OP_SCI, OP_SHP, OP_SI, OP_BRVI, OP_MVARI, OP_NDI, OP_RGBVI,
  OP_TGI, OP_VEG, OP_VNDVI, OP_WI, OP_LAB_L, OP_LAB_A, OP_LAB_B, OP_LAB_BA, OP_LAB_LA,
  OP_LAB_LB, OP_DGCI, OP_ARI, OP_ARVI, OP_BAI, OP_BWDRVI, OP_CCCI, OP_CIG, OP_CIRE,
  OP_CVI, OP_EVI, OP_GARI, OP_GDVI, OP_GEMI, OP_GNDVI, OP_GOSAVI, OP_GRVI, OP_GSAVI,
  OP_IPVI, OP_LAI, OP_MCARI1, OP_MCARI2, OP_MSAVI, OP_MSR, OP_NDMI, OP_NDRE, OP_NDVI,
  OP_NDWI, OP_NLI, OP_OSAVI, OP_PNDVI, OP_PSRI, OP_RDVI, OP_REDVI, OP_RESR, OP_RVI,
  OP_SAVI2, OP_TCARI, OP_TDVI, OP_TSAVI, OP_TVI, OP_VARIRE, OP_VIG, OP_VIN, OP_VIRE,
  OP_WDRVI, OP_WSI, OP_RPN, OP_UNKNOWN
};

inline IndexOp get_index_op(const std::string& ind, bool rpn_valid) {
  if (ind == "R") return OP_R;
  if (ind == "G") return OP_G;
  if (ind == "B") return OP_B;
  if (ind == "RE") return OP_RE;
  if (ind == "NIR") return OP_NIR;
  if (ind == "SWIR") return OP_SWIR;
  if (ind == "NR" || ind == "RCC") return OP_NR;
  if (ind == "NG" || ind == "GCC") return OP_NG;
  if (ind == "NB" || ind == "BCC") return OP_NB;
  if (ind == "GRAY") return OP_GRAY;
  if (ind == "GRAY2") return OP_GRAY2;
  if (ind == "IPCA") return OP_IPCA;
  if (ind == "GLI") return OP_GLI;
  if (ind == "VARI" || ind == "GRVI2" || ind == "NGRDI") return OP_VARI;
  if (ind == "NGBDI") return OP_NGBDI;
  if (ind == "BGI") return OP_BGI;
  if (ind == "BI" || ind == "BI2") return OP_BI;
  if (ind == "CI") return OP_CI;
  if (ind == "CIVE") return OP_CIVE;
  if (ind == "EGVI") return OP_EGVI;
  if (ind == "ERVI") return OP_ERVI;
  if (ind == "GB") return OP_GB;
  if (ind == "GD") return OP_GD;
  if (ind == "GLAI") return OP_GLAI;
  if (ind == "GR") return OP_GR;
  if (ind == "HI") return OP_HI;
  if (ind == "HUE") return OP_HUE;
  if (ind == "HUE2") return OP_HUE2;
  if (ind == "I") return OP_I;
  if (ind == "L") return OP_L;
  if (ind == "MGVRI") return OP_MGVRI;
  if (ind == "RB") return OP_RB;
  if (ind == "RI") return OP_RI;
  if (ind == "S") return OP_S;
  if (ind == "SAVI") return OP_SAVI;
  if (ind == "SCI") return OP_SCI;
  if (ind == "SHP") return OP_SHP;
  if (ind == "SI") return OP_SI;
  if (ind == "BRVI") return OP_BRVI;
  if (ind == "MVARI") return OP_MVARI;
  if (ind == "NDI") return OP_NDI;
  if (ind == "RGBVI") return OP_RGBVI;
  if (ind == "TGI") return OP_TGI;
  if (ind == "VEG") return OP_VEG;
  if (ind == "vNDVI") return OP_VNDVI;
  if (ind == "WI") return OP_WI;
  if (ind == "L*") return OP_LAB_L;
  if (ind == "a") return OP_LAB_A;
  if (ind == "b*") return OP_LAB_B;
  if (ind == "b*-a") return OP_LAB_BA;
  if (ind == "L*-a") return OP_LAB_LA;
  if (ind == "L*-b") return OP_LAB_LB;
  if (ind == "DGCI") return OP_DGCI;
  if (ind == "ARI") return OP_ARI;
  if (ind == "ARVI") return OP_ARVI;
  if (ind == "BAI") return OP_BAI;
  if (ind == "BWDRVI") return OP_BWDRVI;
  if (ind == "CCCI") return OP_CCCI;
  if (ind == "CIG") return OP_CIG;
  if (ind == "CIRE") return OP_CIRE;
  if (ind == "CVI") return OP_CVI;
  if (ind == "EVI") return OP_EVI;
  if (ind == "GARI") return OP_GARI;
  if (ind == "GDVI") return OP_GDVI;
  if (ind == "GEMI") return OP_GEMI;
  if (ind == "GNDVI") return OP_GNDVI;
  if (ind == "GOSAVI") return OP_GOSAVI;
  if (ind == "GRVI") return OP_GRVI;
  if (ind == "GSAVI") return OP_GSAVI;
  if (ind == "IPVI") return OP_IPVI;
  if (ind == "LAI") return OP_LAI;
  if (ind == "MCARI1") return OP_MCARI1;
  if (ind == "MCARI2") return OP_MCARI2;
  if (ind == "MSAVI" || ind == "MSAVI2") return OP_MSAVI;
  if (ind == "MSR") return OP_MSR;
  if (ind == "NDMI") return OP_NDMI;
  if (ind == "NDRE") return OP_NDRE;
  if (ind == "NDVI") return OP_NDVI;
  if (ind == "NDWI") return OP_NDWI;
  if (ind == "NLI") return OP_NLI;
  if (ind == "OSAVI") return OP_OSAVI;
  if (ind == "PNDVI") return OP_PNDVI;
  if (ind == "PSRI") return OP_PSRI;
  if (ind == "RDVI") return OP_RDVI;
  if (ind == "REDVI") return OP_REDVI;
  if (ind == "RESR") return OP_RESR;
  if (ind == "RVI") return OP_RVI;
  if (ind == "SAVI2") return OP_SAVI2;
  if (ind == "TCARI") return OP_TCARI;
  if (ind == "TDVI") return OP_TDVI;
  if (ind == "TSAVI") return OP_TSAVI;
  if (ind == "TVI") return OP_TVI;
  if (ind == "VARIRE") return OP_VARIRE;
  if (ind == "VIG") return OP_VIG;
  if (ind == "VIN") return OP_VIN;
  if (ind == "VIRE") return OP_VIRE;
  if (ind == "WDRVI") return OP_WDRVI;
  if (ind == "WSI") return OP_WSI;
  if (rpn_valid) return OP_RPN;
  return OP_UNKNOWN;
}

inline double exec_pixel_op(IndexOp op, double r, double g, double b, double re, double nir, double swir, const std::vector<RPNToken>& rpn) {
  if (std::isnan(r) || std::isnan(g) || std::isnan(b)) {
    return NA_REAL;
  }
  switch (op) {
    case OP_R: return r;
    case OP_G: return g;
    case OP_B: return b;
    case OP_RE: return re;
    case OP_NIR: return nir;
    case OP_SWIR: return swir;
    case OP_NR: { double s = r + g + b; return (s > 0.0) ? (r / s) : 0.0; }
    case OP_NG: { double s = r + g + b; return (s > 0.0) ? (g / s) : 0.0; }
    case OP_NB: { double s = r + g + b; return (s > 0.0) ? (b / s) : 0.0; }
    case OP_GRAY: return 0.299 * r + 0.587 * g + 0.114 * b;
    case OP_GRAY2: {
      double num = std::pow(r, 2.2) + std::pow(1.5 * g, 2.2) + std::pow(0.6 * b, 2.2);
      double den = 1.0 + std::pow(1.5, 2.2) + std::pow(0.6, 2.2);
      return std::pow(num / den, 1.0 / 2.2);
    }
    case OP_IPCA: return 0.994 * std::abs(r - b) + 0.961 * std::abs(g - b) + 0.914 * std::abs(g - r);
    case OP_GLI: {
      double den = 2.0 * g + r + b;
      return (den != 0.0) ? ((2.0 * g - r - b) / den) : 0.0;
    }
    case OP_VARI: {
      double den = g + r;
      return (den != 0.0) ? ((g - r) / den) : 0.0;
    }
    case OP_NGBDI: {
      double den = g + b;
      return (den != 0.0) ? ((g - b) / den) : 0.0;
    }
    case OP_BGI: return (g != 0.0) ? (b / g) : 0.0;
    case OP_BI: return std::sqrt((r * r + g * g + b * b) / 3.0);
    case OP_CI: return (r != 0.0) ? ((r - b) / r) : 0.0;
    case OP_CIVE: return 0.441 * r - 0.881 * g + 0.385 * b + 18.78745;
    case OP_EGVI: return 2.0 * g - r - b;
    case OP_ERVI: return 1.4 * r - g;
    case OP_GB: return (b != 0.0) ? (g / b) : 0.0;
    case OP_GD: return g - (r + b) / 2.0;
    case OP_GLAI: {
      double den = g + r - b;
      return (den != 0.0) ? (25.0 * (g - r) / den + 1.25) : 0.0;
    }
    case OP_GR: return (r != 0.0) ? (g / r) : 0.0;
    case OP_HI: {
      double den = g - b;
      return (den != 0.0) ? ((2.0 * r - g - b) / den) : 0.0;
    }
    case OP_HUE: return std::atan2(2.0 * (b - g - r), 30.5 * (g - r));
    case OP_HUE2: return std::atan2(2.0 * (r - g - r), 30.5 * (g - b));
    case OP_I: return r + g + b;
    case OP_L: return (r + g + b) / 3.0;
    case OP_MGVRI: {
      double den = g * g + r * r;
      return (den != 0.0) ? ((g * g - r * r) / den) : 0.0;
    }
    case OP_RB: return (b != 0.0) ? (r / b) : 0.0;
    case OP_RI: {
      double den = b * g * g * g;
      return (den != 0.0) ? ((r * r) / den) : 0.0;
    }
    case OP_S: {
      double den = r + g + b;
      return (den != 0.0) ? ((den - 3.0 * b) / den) : 0.0;
    }
    case OP_SAVI: {
      double den = g + r + 0.5;
      return (den != 0.0) ? (1.5 * (g - r) / den) : 0.0;
    }
    case OP_SCI: {
      double den = r + g;
      return (den != 0.0) ? ((r - g) / den) : 0.0;
    }
    case OP_SHP: {
      double den = g - b;
      return (den != 0.0) ? (2.0 * (r - g - b) / den) : 0.0;
    }
    case OP_SI: {
      double den = r + b;
      return (den != 0.0) ? ((r - b) / den) : 0.0;
    }
    case OP_BRVI: {
      double den = b + r;
      return (den != 0.0) ? ((b - r) / den) : 0.0;
    }
    case OP_MVARI: {
      double den = g + r - b;
      return (den != 0.0) ? ((g - b) / den) : 0.0;
    }
    case OP_NDI: {
      double den = g + r;
      return (den != 0.0) ? (128.0 * ((g - r) / den + 1.0)) : 0.0;
    }
    case OP_RGBVI: {
      double den = g * g + b * r;
      return (den != 0.0) ? ((g * g - b * r) / den) : 0.0;
    }
    case OP_TGI: return g - 0.39 * r - 0.61 * b;
    case OP_VEG: {
      double den = std::pow(r, 0.667) * std::pow(b, 0.334);
      return (den != 0.0) ? (g / den) : 0.0;
    }
    case OP_VNDVI: return 0.5268 * (std::pow(r, -0.1294) * std::pow(g, 0.3389) * std::pow(b, -0.3118));
    case OP_WI: {
      double den = r - g;
      return (den != 0.0) ? ((g - b) / den) : 0.0;
    }
    case OP_LAB_L: return 0.2126 * r + 0.7152 * g + 0.0722 * b;
    case OP_LAB_A: {
      double l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      return 0.55 * ((r - l) / (1.0 - 0.2126));
    }
    case OP_LAB_B: {
      double l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      return 0.55 * ((b - l) / (1.0 - 0.0722));
    }
    case OP_LAB_BA: {
      double l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      double a_val = 0.55 * ((r - l) / (1.0 - 0.2126));
      double b_val = 0.55 * ((b - l) / (1.0 - 0.0722));
      return b_val - a_val;
    }
    case OP_LAB_LA: {
      double l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      double a_val = 0.55 * ((r - l) / (1.0 - 0.2126));
      return l - a_val;
    }
    case OP_LAB_LB: {
      double l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      double b_val = 0.55 * ((b - l) / (1.0 - 0.0722));
      return l - b_val;
    }
    case OP_DGCI: {
      double max_v = std::max(r, std::max(g, b));
      double min_v = std::min(r, std::min(g, b));
      double delta = max_v - min_v;
      double val_b = max_v * 100.0;
      double val_s = (max_v > 0.0) ? (delta / max_v * 100.0) : 0.0;
      double val_h = 0.0;
      if (delta > 0.0) {
        if (max_v == r) val_h = 60.0 * std::fmod(((g - b) / delta), 6.0);
        else if (max_v == g) val_h = 60.0 * (((b - r) / delta) + 2.0);
        else val_h = 60.0 * (((r - g) / delta) + 4.0);
      }
      if (val_h < 0.0) val_h += 360.0;
      return ((val_h - 60.0) / 60.0 + (1.0 - val_s / 100.0) + (1.0 - val_b / 100.0)) / 3.0;
    }
    case OP_ARI: return (g != 0.0 && re != 0.0) ? ((1.0 / g) - (1.0 / re)) : 0.0;
    case OP_ARVI: {
      double rb = r - 0.1 * (r - b);
      double den = nir + rb;
      return (den != 0.0) ? ((nir - rb) / den) : 0.0;
    }
    case OP_BAI: {
      double den = (0.1 - r) * (0.1 - r) + (0.06 - nir) * (0.06 - nir);
      return (den != 0.0) ? (1.0 / den) : 0.0;
    }
    case OP_BWDRVI: {
      double den = 0.1 * nir + b;
      return (den != 0.0) ? ((0.1 * nir - b) / den) : 0.0;
    }
    case OP_CCCI: return 1.0;
    case OP_CIG: return (g != 0.0) ? ((nir / g) - 1.0) : 0.0;
    case OP_CIRE: return (re != 0.0) ? ((nir / re) - 1.0) : 0.0;
    case OP_CVI: return (g != 0.0) ? (nir * (r / (g * g))) : 0.0;
    case OP_EVI: {
      double den = nir + 6.0 * r - 7.5 * b + 1.0;
      return (den != 0.0) ? (2.5 * (nir - r) / den) : 0.0;
    }
    case OP_GARI: {
      double br = 1.7 * (b - r);
      double den = nir + br;
      return (den != 0.0) ? ((nir - br) / den) : 0.0;
    }
    case OP_GDVI: return nir - g;
    case OP_GEMI: {
      double den1 = nir + r + 0.5;
      double n1 = (den1 != 0.0) ? ((2.0 * (nir * nir - r * r) + 1.5 * nir + 0.5 * r) / den1) : 0.0;
      double n2 = 1.0 - 0.25 * n1;
      double n3 = (1.0 - r != 0.0) ? ((r - 0.125) / (1.0 - r)) : 0.0;
      return n1 * n2 - n3;
    }
    case OP_GNDVI: {
      double den = nir + g;
      return (den != 0.0) ? ((nir - g) / den) : 0.0;
    }
    case OP_GOSAVI: {
      double den = nir + g + 0.16;
      return (den != 0.0) ? ((nir - g) / den) : 0.0;
    }
    case OP_GRVI: return (g != 0.0) ? (nir / g) : 0.0;
    case OP_GSAVI: {
      double den = nir + g + 0.5;
      return (den != 0.0) ? (((nir - g) / den) * 1.5) : 0.0;
    }
    case OP_IPVI: {
      double den = nir + r;
      return (den != 0.0) ? (nir / den) : 0.0;
    }
    case OP_LAI: {
      double den = nir + 6.0 * r - 7.5 * b + 1.0;
      double evi = (den != 0.0) ? (2.5 * (nir - r) / den) : 0.0;
      return 3.368 * evi - 0.118;
    }
    case OP_MCARI1: return 1.2 * (2.5 * (nir - r) - 1.3 * (nir - g));
    case OP_MCARI2: {
      double num = 1.2 * (2.5 * (nir - r) - 1.3 * (nir - g));
      double sq = std::sqrt(r);
      double den = std::sqrt((2.0 * nir + 1.0) * (2.0 * nir + 1.0) - (6.0 * nir - 5.0 * sq) - 0.5);
      return (den != 0.0) ? (num / den) : 0.0;
    }
    case OP_MSAVI: {
      double term = (2.0 * nir + 1.0) * (2.0 * nir + 1.0) - 8.0 * (nir - r);
      return (term >= 0.0) ? ((2.0 * nir + 1.0 - std::sqrt(term)) / 2.0) : 0.0;
    }
    case OP_MSR: {
      if (r == 0.0) return 0.0;
      double ratio = nir / r;
      double den = std::sqrt(ratio) + 1.0;
      return (den != 0.0) ? ((ratio - 1.0) / den) : 0.0;
    }
    case OP_NDMI: {
      double den = nir + swir;
      return (den != 0.0) ? ((nir - swir) / den) : 0.0;
    }
    case OP_NDRE: {
      double den = nir + re;
      return (den != 0.0) ? ((nir - re) / den) : 0.0;
    }
    case OP_NDVI: {
      double den = nir + r;
      return (den != 0.0) ? ((nir - r) / den) : 0.0;
    }
    case OP_NDWI: {
      double den = g + nir;
      return (den != 0.0) ? ((g - nir) / den) : 0.0;
    }
    case OP_NLI: {
      double den = nir * nir + r;
      return (den != 0.0) ? ((nir * nir - r) / den) : 0.0;
    }
    case OP_OSAVI: {
      double den = nir + r + 0.16;
      return (den != 0.0) ? ((nir - r) / den) : 0.0;
    }
    case OP_PNDVI: {
      double den = nir + g + r + b;
      return (den != 0.0) ? ((nir - (g + r + b)) / den) : 0.0;
    }
    case OP_PSRI: return (re != 0.0) ? ((r - g) / re) : 0.0;
    case OP_RDVI: {
      double den = std::sqrt(nir + r);
      return (den != 0.0) ? ((nir - r) / den) : 0.0;
    }
    case OP_REDVI: return nir - re;
    case OP_RESR: return (re != 0.0) ? (nir / re) : 0.0;
    case OP_RVI: return (nir != 0.0) ? (r / nir) : 0.0;
    case OP_SAVI2: {
      double den = nir + r + 0.5;
      return (den != 0.0) ? (((nir - r) / den) + 0.5) : 0.0;
    }
    case OP_TCARI: return (r != 0.0) ? (3.0 * ((re - r) - 0.2 * (re - g) * (re / r))) : 0.0;
    case OP_TDVI: {
      double den = std::sqrt(nir * nir + r + 0.5);
      return (den != 0.0) ? (1.5 * (nir - r) / den) : 0.0;
    }
    case OP_TSAVI: {
      double den = r + 2.0 * (nir - 1.0) + 2.5;
      return (den != 0.0) ? ((2.0 * (nir - 2.0) * (r - 1.0)) / den) : 0.0;
    }
    case OP_TVI: {
      double den = nir + r;
      double val = (den != 0.0) ? ((nir - r) / den + 0.5) : 0.5;
      return (val >= 0.0) ? std::sqrt(val) : 0.0;
    }
    case OP_VARIRE: {
      double num = re - 1.7 * r + 0.7 * b;
      double den = re + 2.3 * r - 1.3 * b;
      return (den != 0.0) ? (num / den) : 0.0;
    }
    case OP_VIG: {
      double den = g + r;
      return (den != 0.0) ? ((g - r) / den) : 0.0;
    }
    case OP_VIN: return (r != 0.0) ? (nir / r) : 0.0;
    case OP_VIRE: {
      double den = re + r;
      return (den != 0.0) ? ((re - r) / den) : 0.0;
    }
    case OP_WDRVI: {
      double den = 0.2 * nir + r;
      return (den != 0.0) ? ((0.2 * nir - r) / den) : 0.0;
    }
    case OP_WSI: {
      double den1 = nir + r;
      double den2 = re + nir;
      double t1 = (den1 != 0.0) ? (5.0 * (nir - r) / den1) : 0.0;
      double t2 = (den2 != 0.0) ? ((re - nir) / den2) : 0.0;
      return t1 - t2;
    }
    case OP_RPN: return eval_rpn(rpn, r, g, b, re, nir, swir);
    default: return 1e9;
  }
}

// [[Rcpp::export]]
SEXP compute_single_index_cpp(SEXP img_sexp, std::string ind, int r = 1, int g = 2, int b = 3, int re = 4, int nir = 5, int swir = 6, std::string storage = "auto") {
  SEXP dims_sexp = Rf_getAttrib(img_sexp, R_DimSymbol);
  int npix = 0;
  int w = 0, h = 0, nch = 1;

  if (!Rf_isNull(dims_sexp)) {
    int ndim = Rf_length(dims_sexp);
    w = INTEGER(dims_sexp)[0];
    h = INTEGER(dims_sexp)[1];
    nch = (ndim >= 3) ? INTEGER(dims_sexp)[2] : 1;
    npix = w * h;
  } else {
    npix = Rf_length(img_sexp);
    w = npix;
    h = 1;
  }

  if (npix == 0) return R_NilValue;

  int r_idx = std::max(1, std::min(nch, r)) - 1;
  int g_idx = std::max(1, std::min(nch, g)) - 1;
  int b_idx = std::max(1, std::min(nch, b)) - 1;
  int re_idx = std::max(1, std::min(nch, re)) - 1;
  int nir_idx = std::max(1, std::min(nch, nir)) - 1;
  int swir_idx = std::max(1, std::min(nch, swir)) - 1;

  std::vector<RPNToken> rpn_tokens;
  bool rpn_valid = false;
  try {
    rpn_tokens = parse_to_rpn(ind);
    if (!rpn_tokens.empty()) rpn_valid = true;
  } catch (...) {
    rpn_valid = false;
  }

  IndexOp op = get_index_op(ind, rpn_valid);
  if (op == OP_UNKNOWN) return R_NilValue;

  bool is_band = (ind == "R" || ind == "G" || ind == "B" || ind == "RE" || ind == "NIR" || ind == "SWIR");
  bool is_gray  = (ind == "GRAY");
  bool return_raw = false;
  if (storage == "raw") {
    return_raw = true;
  } else if (storage == "auto" || storage == "" || storage == "NULL") {
    if ((is_band || is_gray) && TYPEOF(img_sexp) == RAWSXP) {
      return_raw = true;
    }
  }

  if (return_raw) {
    SEXP res = PROTECT(Rf_allocMatrix(RAWSXP, w, h));
    if (!Rf_isNull(dims_sexp) && Rf_length(dims_sexp) >= 2) {
      SEXP out_dims = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(out_dims)[0] = w;
      INTEGER(out_dims)[1] = h;
      Rf_setAttrib(res, R_DimSymbol, out_dims);
      UNPROTECT(1);
    }
    uint8_t* out = RAW(res);

    if (TYPEOF(img_sexp) == RAWSXP) {
      const uint8_t* ptr = RAW(img_sexp);
      const uint8_t* pChan = (ind == "R") ? ptr + r_idx * npix :
                             (ind == "G") ? ptr + g_idx * npix :
                             (ind == "B") ? ptr + b_idx * npix :
                             (ind == "RE") ? ptr + re_idx * npix :
                             (ind == "NIR") ? ptr + nir_idx * npix :
                             (ind == "SWIR") ? ptr + swir_idx * npix : NULL;
      if (pChan != NULL) {
        std::memcpy(out, pChan, npix);
      } else if (is_gray) {
        const uint8_t* pR = ptr + r_idx * npix;
        const uint8_t* pG = ptr + g_idx * npix;
        const uint8_t* pB = ptr + b_idx * npix;
        #pragma omp parallel for if(npix > 250000)
        for (int i = 0; i < npix; i++) {
          int v = (int)std::round(0.299 * pR[i] + 0.587 * pG[i] + 0.114 * pB[i]);
          out[i] = (uint8_t)std::max(0, std::min(255, v));
        }
      }
    } else if (TYPEOF(img_sexp) == REALSXP) {
      const double* ptr = REAL(img_sexp);
      const double* pChan = (ind == "R") ? ptr + r_idx * npix :
                            (ind == "G") ? ptr + g_idx * npix :
                            (ind == "B") ? ptr + b_idx * npix :
                            (ind == "RE") ? ptr + re_idx * npix :
                            (ind == "NIR") ? ptr + nir_idx * npix :
                            (ind == "SWIR") ? ptr + swir_idx * npix : NULL;
      if (pChan != NULL) {
        #pragma omp parallel for if(npix > 250000)
        for (int i = 0; i < npix; i++) {
          int v = (int)std::round(pChan[i] * 255.0);
          out[i] = (uint8_t)std::max(0, std::min(255, v));
        }
      } else if (is_gray) {
        const double* pR = ptr + r_idx * npix;
        const double* pG = ptr + g_idx * npix;
        const double* pB = ptr + b_idx * npix;
        #pragma omp parallel for if(npix > 250000)
        for (int i = 0; i < npix; i++) {
          int v = (int)std::round((0.299 * pR[i] + 0.587 * pG[i] + 0.114 * pB[i]) * 255.0);
          out[i] = (uint8_t)std::max(0, std::min(255, v));
        }
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

  Rcpp::NumericMatrix res(w, h);
  double* out = REAL(res);

  if (TYPEOF(img_sexp) == RAWSXP) {
    const uint8_t* ptr = RAW(img_sexp);
    const uint8_t* pR = ptr + r_idx * npix;
    const uint8_t* pG = ptr + g_idx * npix;
    const uint8_t* pB = ptr + b_idx * npix;
    const uint8_t* pRE = ptr + re_idx * npix;
    const uint8_t* pNIR = ptr + nir_idx * npix;
    const uint8_t* pSWIR = ptr + swir_idx * npix;

    if (op == OP_NB) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
        out[i] = (s > 0.0) ? ((double)pB[i] / s) : 0.0;
      }
    } else if (op == OP_NR) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
        out[i] = (s > 0.0) ? ((double)pR[i] / s) : 0.0;
      }
    } else if (op == OP_NG) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
        out[i] = (s > 0.0) ? ((double)pG[i] / s) : 0.0;
      }
    } else if (op == OP_GRAY) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        out[i] = (0.299 * pR[i] + 0.587 * pG[i] + 0.114 * pB[i]) / 255.0;
      }
    } else if (op == OP_R) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pR[i] / 255.0;
    } else if (op == OP_G) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pG[i] / 255.0;
    } else if (op == OP_B) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pB[i] / 255.0;
    } else {
      #pragma omp parallel for if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        out[i] = exec_pixel_op(op, pR[i] / 255.0, pG[i] / 255.0, pB[i] / 255.0, pRE[i] / 255.0, pNIR[i] / 255.0, pSWIR[i] / 255.0, rpn_tokens);
      }
    }
  } else if (TYPEOF(img_sexp) == REALSXP) {
    const double* ptr = REAL(img_sexp);
    const double* pR = ptr + r_idx * npix;
    const double* pG = ptr + g_idx * npix;
    const double* pB = ptr + b_idx * npix;
    const double* pRE = ptr + re_idx * npix;
    const double* pNIR = ptr + nir_idx * npix;
    const double* pSWIR = ptr + swir_idx * npix;

    if (op == OP_NB) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (std::isnan(pR[i]) || std::isnan(pG[i]) || std::isnan(pB[i])) {
          out[i] = NA_REAL;
        } else {
          double s = pR[i] + pG[i] + pB[i];
          out[i] = (s > 0.0) ? (pB[i] / s) : 0.0;
        }
      }
    } else if (op == OP_NR) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (std::isnan(pR[i]) || std::isnan(pG[i]) || std::isnan(pB[i])) {
          out[i] = NA_REAL;
        } else {
          double s = pR[i] + pG[i] + pB[i];
          out[i] = (s > 0.0) ? (pR[i] / s) : 0.0;
        }
      }
    } else if (op == OP_NG) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (std::isnan(pR[i]) || std::isnan(pG[i]) || std::isnan(pB[i])) {
          out[i] = NA_REAL;
        } else {
          double s = pR[i] + pG[i] + pB[i];
          out[i] = (s > 0.0) ? (pG[i] / s) : 0.0;
        }
      }
    } else if (op == OP_GRAY) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (std::isnan(pR[i]) || std::isnan(pG[i]) || std::isnan(pB[i])) {
          out[i] = NA_REAL;
        } else {
          out[i] = 0.299 * pR[i] + 0.587 * pG[i] + 0.114 * pB[i];
        }
      }
    } else if (op == OP_R) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pR[i];
    } else if (op == OP_G) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pG[i];
    } else if (op == OP_B) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = pB[i];
    } else {
      #pragma omp parallel for if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        out[i] = exec_pixel_op(op, pR[i], pG[i], pB[i], pRE[i], pNIR[i], pSWIR[i], rpn_tokens);
      }
    }
  } else if (TYPEOF(img_sexp) == INTSXP) {
    const int* ptr = INTEGER(img_sexp);
    const int* pR = ptr + r_idx * npix;
    const int* pG = ptr + g_idx * npix;
    const int* pB = ptr + b_idx * npix;
    const int* pRE = ptr + re_idx * npix;
    const int* pNIR = ptr + nir_idx * npix;
    const int* pSWIR = ptr + swir_idx * npix;

    if (op == OP_NB) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (pR[i] == NA_INTEGER || pG[i] == NA_INTEGER || pB[i] == NA_INTEGER) {
          out[i] = NA_REAL;
        } else {
          double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
          out[i] = (s > 0.0) ? ((double)pB[i] / s) : 0.0;
        }
      }
    } else if (op == OP_NR) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (pR[i] == NA_INTEGER || pG[i] == NA_INTEGER || pB[i] == NA_INTEGER) {
          out[i] = NA_REAL;
        } else {
          double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
          out[i] = (s > 0.0) ? ((double)pR[i] / s) : 0.0;
        }
      }
    } else if (op == OP_NG) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (pR[i] == NA_INTEGER || pG[i] == NA_INTEGER || pB[i] == NA_INTEGER) {
          out[i] = NA_REAL;
        } else {
          double s = (double)pR[i] + (double)pG[i] + (double)pB[i];
          out[i] = (s > 0.0) ? ((double)pG[i] / s) : 0.0;
        }
      }
    } else if (op == OP_GRAY) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (pR[i] == NA_INTEGER || pG[i] == NA_INTEGER || pB[i] == NA_INTEGER) {
          out[i] = NA_REAL;
        } else {
          out[i] = (0.299 * (double)pR[i] + 0.587 * (double)pG[i] + 0.114 * (double)pB[i]) / 255.0;
        }
      }
    } else if (op == OP_R) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = (pR[i] == NA_INTEGER) ? NA_REAL : (double)pR[i] / 255.0;
    } else if (op == OP_G) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = (pG[i] == NA_INTEGER) ? NA_REAL : (double)pG[i] / 255.0;
    } else if (op == OP_B) {
      #pragma omp parallel for schedule(static) if(npix > 20000)
      for (int i = 0; i < npix; i++) out[i] = (pB[i] == NA_INTEGER) ? NA_REAL : (double)pB[i] / 255.0;
    } else {
      #pragma omp parallel for if(npix > 20000)
      for (int i = 0; i < npix; i++) {
        if (pR[i] == NA_INTEGER || pG[i] == NA_INTEGER || pB[i] == NA_INTEGER) {
          out[i] = NA_REAL;
        } else {
          out[i] = exec_pixel_op(op, (double)pR[i] / 255.0, (double)pG[i] / 255.0, (double)pB[i] / 255.0, (double)pRE[i] / 255.0, (double)pNIR[i] / 255.0, (double)pSWIR[i] / 255.0, rpn_tokens);
        }
      }
    }
  }

  res.attr("class") = Rcpp::CharacterVector::create("image", "array");
  res.attr("colormode") = "Grayscale";
  return res;
}
