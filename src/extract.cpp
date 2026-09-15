#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#if defined(__APPLE__) || defined(__MACH__)
#include <mach/mach.h>
#include <mach/host_info.h>
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//' Fast C++ retrieval of available system physical RAM in bytes (Windows, Linux, macOS)
//'
//' @export
// [[Rcpp::export]]
double get_free_ram_cpp() {
#ifdef _WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return (double)status.ullAvailPhys;
    }
#elif defined(__APPLE__) || defined(__MACH__)
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
    vm_statistics_data_t vm_stat;
    if (host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vm_stat, &count) == KERN_SUCCESS) {
        long page_size = sysconf(_SC_PAGE_SIZE);
        return (double)(vm_stat.free_count + vm_stat.inactive_count) * (double)page_size;
    }
#elif defined(__linux__)
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0) {
        return (double)pages * (double)page_size;
    }
#endif
    return -1.0;
}

// Geometry structures for C++ scanline rasterization
struct Ring {
    std::vector<double> x;
    std::vector<double> y;
    int n;
};

struct PolyPart {
    Ring outer;
    std::vector<Ring> holes;
    double bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax;
};

struct FeatureGeom {
    std::vector<PolyPart> parts;
    double bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax;
};

// Thread-safe C++ container for feature raw pixel results
struct RawFeatureResult {
    std::vector<int> ext_cell;
    std::vector<int> ext_r;
    std::vector<int> ext_c;
    std::vector<double> ext_area;
    std::vector<std::vector<double>> ext_vals;
};

// Zero-allocation stack buffer for scanline intersections
struct ScanlineBuffer {
    double x_ints[128];
    std::pair<double, double> intervals[64];
    int n_ints;
    int n_inter;
};

inline void get_scanline_intervals_stack(const Ring& ring, double y, ScanlineBuffer& buf) {
    buf.n_ints = 0;
    buf.n_inter = 0;
    const double* rx = ring.x.data();
    const double* ry = ring.y.data();
    int n = ring.n;

    for (int i = 0, j = n - 1; i < n; j = i++) {
        double y1 = ry[i], y2 = ry[j];
        if ((y1 > y) != (y2 > y)) {
            double x1 = rx[i], x2 = rx[j];
            double x_int = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
            if (buf.n_ints < 128) {
                buf.x_ints[buf.n_ints++] = x_int;
            }
        }
    }
    if (buf.n_ints < 2) return;

    std::sort(buf.x_ints, buf.x_ints + buf.n_ints);

    for (int k = 0; k + 1 < buf.n_ints; k += 2) {
        if (buf.n_inter < 64) {
            buf.intervals[buf.n_inter++] = {buf.x_ints[k], buf.x_ints[k + 1]};
        }
    }
}

// Fast exact analytical 1D area fraction calculation across scanline cell
inline double compute_exact_boundary_fraction(double cell_xmin, double cell_xmax, double x_left, double x_right, double res_x) {
    double overlap_min = cell_xmin > x_left ? cell_xmin : x_left;
    double overlap_max = cell_xmax < x_right ? cell_xmax : x_right;
    double len = overlap_max - overlap_min;
    if (len <= 0.0) return 0.0;
    double frac = len / res_x;
    return frac > 1.0 ? 1.0 : frac;
}

// High-precision 2D boundary cell coverage fraction integration
inline double compute_exact_cell_2d_fraction(double cell_xmin, double cell_xmax, double cell_ymin, double cell_ymax, double res_x, double res_y, const PolyPart& part, int subdiv = 10) {
    double dy = res_y / subdiv;
    double frac_sum = 0.0;
    ScanlineBuffer sub_buf;

    for (int s = 0; s < subdiv; ++s) {
        double sub_y = cell_ymax - (s + 0.5) * dy;
        if (sub_y < part.bbox_ymin || sub_y > part.bbox_ymax) continue;

        get_scanline_intervals_stack(part.outer, sub_y, sub_buf);
        if (sub_buf.n_inter == 0) continue;

        for (int k = 0; k < sub_buf.n_inter; ++k) {
            const auto& inter = sub_buf.intervals[k];
            double fx = compute_exact_boundary_fraction(cell_xmin, cell_xmax, inter.first, inter.second, res_x);
            frac_sum += fx;
        }
    }
    double total_frac = frac_sum / subdiv;
    return total_frac > 1.0 ? 1.0 : total_frac;
}

// Helper to get full row outer X-extent and inner pure-interior X-extent for row [row_ymin, row_ymax]
inline bool get_row_scanline_extents(const PolyPart& part, double row_ymin, double row_ymax, double cell_yc, double& outer_min_x, double& outer_max_x, double& inner_min_x, double& inner_max_x) {
    ScanlineBuffer b_mid, b_top, b_bot;
    double eps_y = (row_ymax - row_ymin) * 1e-4;

    get_scanline_intervals_stack(part.outer, cell_yc, b_mid);
    get_scanline_intervals_stack(part.outer, row_ymax - eps_y, b_top);
    get_scanline_intervals_stack(part.outer, row_ymin + eps_y, b_bot);

    if (b_mid.n_inter == 0 && b_top.n_inter == 0 && b_bot.n_inter == 0) {
        return false;
    }

    outer_min_x = 1e30;
    outer_max_x = -1e30;

    auto expand_outer = [&](const ScanlineBuffer& buf) {
        for (int k = 0; k < buf.n_inter; ++k) {
            if (buf.intervals[k].first < outer_min_x) outer_min_x = buf.intervals[k].first;
            if (buf.intervals[k].second > outer_max_x) outer_max_x = buf.intervals[k].second;
        }
    };
    expand_outer(b_mid);
    expand_outer(b_top);
    expand_outer(b_bot);

    if (b_top.n_inter > 0 && b_bot.n_inter > 0 && b_mid.n_inter > 0) {
        inner_min_x = std::max({b_top.intervals[0].first, b_mid.intervals[0].first, b_bot.intervals[0].first});
        inner_max_x = std::min({b_top.intervals[0].second, b_mid.intervals[0].second, b_bot.intervals[0].second});
    } else {
        inner_min_x = 1e30;
        inner_max_x = -1e30;
    }
    return true;
}

// Helper to compute continuous weighted quantile via linear interpolation
inline double compute_weighted_quantile(const std::vector<std::pair<double, double>>& sorted_pixels, double W, double q_prob) {
    int n_p = sorted_pixels.size();
    if (n_p == 0) return NA_REAL;
    if (n_p == 1 || W <= 0.0) return sorted_pixels[0].first;

    std::vector<double> p_mid(n_p);
    double cum_w = 0.0;
    for (int k = 0; k < n_p; ++k) {
        p_mid[k] = (cum_w + 0.5 * sorted_pixels[k].second) / W;
        cum_w += sorted_pixels[k].second;
    }

    if (q_prob <= p_mid[0]) return sorted_pixels[0].first;
    if (q_prob >= p_mid[n_p - 1]) return sorted_pixels[n_p - 1].first;

    auto it = std::upper_bound(p_mid.begin(), p_mid.end(), q_prob);
    int k2 = std::distance(p_mid.begin(), it);
    int k1 = k2 - 1;
    double t = (q_prob - p_mid[k1]) / (p_mid[k2] - p_mid[k1]);
    return (1.0 - t) * sorted_pixels[k1].first + t * sorted_pixels[k2].first;
}

//' Fast C++ raster extraction engine without GDAL
//'
//' @param values Matrix of raster values (ncell x nlyr).
//' @param n_row Number of raster rows.
//' @param n_col Number of raster columns.
//' @param n_lyr Number of raster layers.
//' @param bbox_raster Vector c(xmin, xmax, ymin, ymax).
//' @param geoms_r Rcpp::List of polygon geometries for each feature.
//' @param fun CharacterVector of summary functions ("none", "mean", "median", "sum", "min", "max", "sd", "count", "quantile", "quantiles").
//' @param exact Logical, whether to compute exact pixel area fraction coverage.
//' @param return_coverage_area Logical, whether to include coverage_area column in raw pixel output.
//' @param subdiv Number of subdivisions per dimension when exact = TRUE.
//' @param summarize_quantiles NumericVector of quantile probabilities (between 0 and 1) when "quantiles" is requested.
//' @export
// [[Rcpp::export]]
RObject cpp_extract_raster(NumericMatrix values,
                           int n_row,
                           int n_col,
                           int n_lyr,
                           NumericVector bbox_raster,
                           List geoms_r,
                           CharacterVector fun = CharacterVector::create("none"),
                           bool exact = false,
                           bool return_coverage_area = false,
                           int subdiv = 10,
                           NumericVector summarize_quantiles = NumericVector::create(0.05, 0.975)) {
    int n_feat = geoms_r.size();
    std::vector<FeatureGeom> feats(n_feat);

    for (int i = 0; i < n_feat; ++i) {
        List f_list = geoms_r[i];
        List parts_r = f_list["parts"];
        int n_parts = parts_r.size();

        double f_xmin = 1e30, f_xmax = -1e30, f_ymin = 1e30, f_ymax = -1e30;

        for (int p = 0; p < n_parts; ++p) {
            List p_list = parts_r[p];
            NumericMatrix outer_m = p_list["outer"];
            List holes_l = p_list["holes"];

            PolyPart part;
            part.outer.n = outer_m.nrow();
            part.outer.x.resize(part.outer.n);
            part.outer.y.resize(part.outer.n);

            double p_xmin = 1e30, p_xmax = -1e30, p_ymin = 1e30, p_ymax = -1e30;
            for (int pt = 0; pt < part.outer.n; ++pt) {
                double vx = outer_m(pt, 0);
                double vy = outer_m(pt, 1);
                part.outer.x[pt] = vx;
                part.outer.y[pt] = vy;
                if (vx < p_xmin) p_xmin = vx;
                if (vx > p_xmax) p_xmax = vx;
                if (vy < p_ymin) p_ymin = vy;
                if (vy > p_ymax) p_ymax = vy;
            }
            part.bbox_xmin = p_xmin;
            part.bbox_xmax = p_xmax;
            part.bbox_ymin = p_ymin;
            part.bbox_ymax = p_ymax;

            if (p_xmin < f_xmin) f_xmin = p_xmin;
            if (p_xmax > f_xmax) f_xmax = p_xmax;
            if (p_ymin < f_ymin) f_ymin = p_ymin;
            if (p_ymax > f_ymax) f_ymax = p_ymax;

            int n_holes = holes_l.size();
            for (int h = 0; h < n_holes; ++h) {
                NumericMatrix hole_m = holes_l[h];
                Ring hole_ring;
                hole_ring.n = hole_m.nrow();
                hole_ring.x.resize(hole_ring.n);
                hole_ring.y.resize(hole_ring.n);
                for (int pt = 0; pt < hole_ring.n; ++pt) {
                    hole_ring.x[pt] = hole_m(pt, 0);
                    hole_ring.y[pt] = hole_m(pt, 1);
                }
                part.holes.push_back(hole_ring);
            }
            feats[i].parts.push_back(part);
        }
        feats[i].bbox_xmin = f_xmin;
        feats[i].bbox_xmax = f_xmax;
        feats[i].bbox_ymin = f_ymin;
        feats[i].bbox_ymax = f_ymax;
    }

    double r_xmin = bbox_raster[0];
    double r_xmax = bbox_raster[1];
    double r_ymin = bbox_raster[2];
    double r_ymax = bbox_raster[3];

    double res_x = (r_xmax - r_xmin) / n_col;
    double res_y = (r_ymax - r_ymin) / n_row;
    double cell_area = std::abs(res_x * res_y);

    int n_funs = fun.size();
    bool is_summary = true;
    if (n_funs == 1 && std::string(fun[0]) == "none") {
        is_summary = false;
    }

    const double* __restrict val_ptr = REAL(values);
    int n_cell_total = n_row * n_col;

    if (is_summary) {
        int n_quantiles = summarize_quantiles.size();
        std::vector<double> probs(n_quantiles);
        for (int q = 0; q < n_quantiles; ++q) probs[q] = summarize_quantiles[q];

        int total_output_cols = 0;
        for (int f = 0; f < n_funs; ++f) {
            std::string s = Rcpp::as<std::string>(fun[f]);
            if (s == "quantile" || s == "quantiles") {
                total_output_cols += n_lyr * n_quantiles;
            } else {
                total_output_cols += n_lyr;
            }
        }
        if (return_coverage_area) {
            total_output_cols += 3; // covered_area, plot_area, coverage
        }

        NumericMatrix summary_mat(n_feat, total_output_cols);
        std::fill(summary_mat.begin(), summary_mat.end(), NA_REAL);

        bool need_sorted_pixels = false;
        bool need_sd = false;
        for (int f = 0; f < n_funs; ++f) {
            std::string s = Rcpp::as<std::string>(fun[f]);
            if (s == "median" || s == "quantile" || s == "quantiles") need_sorted_pixels = true;
            if (s == "sd") need_sd = true;
        }

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_feat; ++i) {
            const auto& feat = feats[i];

            int c_start = std::max(0, (int)std::floor((feat.bbox_xmin - r_xmin) / res_x));
            int c_end   = std::min(n_col - 1, (int)std::floor((feat.bbox_xmax - r_xmin) / res_x));
            int r_start = std::max(0, (int)std::floor((r_ymax - feat.bbox_ymax) / res_y));
            int r_end   = std::min(n_row - 1, (int)std::floor((r_ymax - feat.bbox_ymin) / res_y));

            if (c_start > c_end || r_start > r_end) {
                continue;
            }

            std::vector<double> w_sum(n_lyr, 0.0);
            std::vector<double> val_sum(n_lyr, 0.0);
            std::vector<double> min_v(n_lyr, 1e30);
            std::vector<double> max_v(n_lyr, -1e30);
            std::vector<int> valid_cnt(n_lyr, 0);

            double total_w_poly = 0.0;
            double valid_w_poly = 0.0;

            std::vector<std::vector<std::pair<double, double>>> layer_pixels(n_lyr);

            ScanlineBuffer scanline_buf;
            int cell_count = 0;

            auto process_pixel = [&](int r, int c, double frac) {
                cell_count++;
                total_w_poly += frac;
                int cell_idx = r * n_col + c;

                bool is_valid_any = false;
                for (int l = 0; l < n_lyr; ++l) {
                    double v = val_ptr[l * n_cell_total + cell_idx];
                    if (R_IsNA(v)) continue;
                    is_valid_any = true;
                    valid_cnt[l]++;
                    val_sum[l] += v * frac;
                    w_sum[l] += frac;
                    if (v < min_v[l]) min_v[l] = v;
                    if (v > max_v[l]) max_v[l] = v;
                    if (need_sorted_pixels || need_sd) {
                        layer_pixels[l].push_back({v, frac});
                    }
                }
                if (is_valid_any) {
                    valid_w_poly += frac;
                }
            };

            for (int r = r_start; r <= r_end; ++r) {
                double cell_yc = r_ymax - (r + 0.5) * res_y;
                double row_ymin = r_ymax - (r + 1) * res_y;
                double row_ymax = r_ymax - r * res_y;

                for (const auto& part : feat.parts) {
                    if (row_ymax < part.bbox_ymin || row_ymin > part.bbox_ymax) continue;

                    if (!exact) {
                        get_scanline_intervals_stack(part.outer, cell_yc, scanline_buf);
                        if (scanline_buf.n_inter == 0) continue;
                        for (int k = 0; k < scanline_buf.n_inter; ++k) {
                            const auto& inter = scanline_buf.intervals[k];
                            int c_in_start = std::max(c_start, (int)std::floor((inter.first - r_xmin) / res_x));
                            int c_in_end   = std::min(c_end, (int)std::floor((inter.second - r_xmin) / res_x));
                            for (int c = c_in_start; c <= c_in_end; ++c) {
                                process_pixel(r, c, 1.0);
                            }
                        }
                    } else {
                        double outer_min_x, outer_max_x, inner_min_x, inner_max_x;
                        if (!get_row_scanline_extents(part, row_ymin, row_ymax, cell_yc, outer_min_x, outer_max_x, inner_min_x, inner_max_x)) continue;

                        int c_in_start = std::max(c_start, (int)std::floor((outer_min_x - r_xmin) / res_x));
                        int c_in_end   = std::min(c_end, (int)std::floor((outer_max_x - r_xmin) / res_x));

                        int c_int_start = std::max(c_in_start, (int)std::ceil((inner_min_x - r_xmin) / res_x));
                        int c_int_end   = std::min(c_in_end, (int)std::floor((inner_max_x - r_xmin) / res_x) - 1);

                        bool is_pure_interior_row = (row_ymax <= part.bbox_ymax && row_ymin >= part.bbox_ymin);

                        if (!part.holes.empty() || c_int_start > c_int_end || !is_pure_interior_row) {
                            for (int c = c_in_start; c <= c_in_end; ++c) {
                                double cell_xmin = r_xmin + c * res_x;
                                double cell_xmax = r_xmin + (c + 1) * res_x;
                                double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                if (frac <= 0.0) continue;
                                process_pixel(r, c, frac);
                            }
                        } else {
                            for (int c = c_in_start; c < c_int_start && c <= c_in_end; ++c) {
                                double cell_xmin = r_xmin + c * res_x;
                                double cell_xmax = r_xmin + (c + 1) * res_x;
                                double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                if (frac <= 0.0) continue;
                                process_pixel(r, c, frac);
                            }
                            for (int c = c_int_start; c <= c_int_end; ++c) {
                                process_pixel(r, c, 1.0);
                            }
                            for (int c = std::max(c_in_start, c_int_end + 1); c <= c_in_end; ++c) {
                                double cell_xmin = r_xmin + c * res_x;
                                double cell_xmax = r_xmin + (c + 1) * res_x;
                                double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                if (frac <= 0.0) continue;
                                process_pixel(r, c, frac);
                            }
                        }
                    }
                }
            }

            if (cell_count == 0) continue;

            // Sort layer pixels once if needed
            if (need_sorted_pixels) {
                for (int l = 0; l < n_lyr; ++l) {
                    if (!layer_pixels[l].empty()) {
                        std::sort(layer_pixels[l].begin(), layer_pixels[l].end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                            return a.first < b.first;
                        });
                    }
                }
            }

            int col_offset = 0;
            for (int f_idx = 0; f_idx < n_funs; ++f_idx) {
                std::string f_str = Rcpp::as<std::string>(fun[f_idx]);

                if (f_str == "quantile" || f_str == "quantiles") {
                    for (int q_idx = 0; q_idx < n_quantiles; ++q_idx) {
                        double q_prob = probs[q_idx];
                        for (int l = 0; l < n_lyr; ++l) {
                            int cur_col = col_offset + q_idx * n_lyr + l;
                            if (valid_cnt[l] == 0 || layer_pixels[l].empty()) continue;

                            if (!exact) {
                                int n_p = layer_pixels[l].size();
                                double pos = q_prob * (n_p - 1);
                                int idx = (int)std::floor(pos);
                                double frac = pos - idx;
                                if (idx + 1 < n_p) {
                                    summary_mat(i, cur_col) = (1.0 - frac) * layer_pixels[l][idx].first + frac * layer_pixels[l][idx + 1].first;
                                } else {
                                    summary_mat(i, cur_col) = layer_pixels[l][idx].first;
                                }
                            } else {
                                summary_mat(i, cur_col) = compute_weighted_quantile(layer_pixels[l], w_sum[l], q_prob);
                            }
                        }
                    }
                    col_offset += n_lyr * n_quantiles;
                } else {
                    for (int l = 0; l < n_lyr; ++l) {
                        int cur_col = col_offset + l;
                        if (valid_cnt[l] == 0) continue;

                        if (f_str == "mean") {
                            if (w_sum[l] > 0.0) summary_mat(i, cur_col) = val_sum[l] / w_sum[l];
                        } else if (f_str == "sum") {
                            summary_mat(i, cur_col) = val_sum[l];
                        } else if (f_str == "min") {
                            summary_mat(i, cur_col) = min_v[l];
                        } else if (f_str == "max") {
                            summary_mat(i, cur_col) = max_v[l];
                        } else if (f_str == "count") {
                            summary_mat(i, cur_col) = (double)valid_cnt[l];
                        } else if (f_str == "sd") {
                            if (layer_pixels[l].empty()) continue;
                            double mean_val = val_sum[l] / w_sum[l];
                            double sq_diff_sum = 0.0;
                            for (const auto& vw : layer_pixels[l]) {
                                sq_diff_sum += vw.second * (vw.first - mean_val) * (vw.first - mean_val);
                            }
                            double denom = (exact ? w_sum[l] : (double)(valid_cnt[l] - 1));
                            if (denom > 0.0) {
                                summary_mat(i, cur_col) = std::sqrt(sq_diff_sum / denom);
                            }
                        } else if (f_str == "median") {
                            if (layer_pixels[l].empty()) continue;
                            if (!exact) {
                                int n_p = layer_pixels[l].size();
                                if (n_p % 2 == 1) {
                                    summary_mat(i, cur_col) = layer_pixels[l][n_p / 2].first;
                                } else {
                                    summary_mat(i, cur_col) = 0.5 * (layer_pixels[l][n_p / 2 - 1].first + layer_pixels[l][n_p / 2].first);
                                }
                            } else {
                                summary_mat(i, cur_col) = compute_weighted_quantile(layer_pixels[l], w_sum[l], 0.5);
                            }
                        }
                    }
                    col_offset += n_lyr;
                }
            }

            if (return_coverage_area) {
                double plot_area_val = total_w_poly * cell_area;
                double covered_area_val = valid_w_poly * cell_area;
                double coverage_val = (plot_area_val > 0.0) ? (covered_area_val / plot_area_val) : NA_REAL;

                summary_mat(i, total_output_cols - 3) = covered_area_val;
                summary_mat(i, total_output_cols - 2) = plot_area_val;
                summary_mat(i, total_output_cols - 1) = coverage_val;
            }
        }
        return summary_mat;
    } else {
        std::vector<RawFeatureResult> feat_results(n_feat);

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < n_feat; ++i) {
            const auto& feat = feats[i];
            auto& fr = feat_results[i];
            fr.ext_vals.resize(n_lyr);

            int c_start = std::max(0, (int)std::floor((feat.bbox_xmin - r_xmin) / res_x));
            int c_end   = std::min(n_col - 1, (int)std::floor((feat.bbox_xmax - r_xmin) / res_x));
            int r_start = std::max(0, (int)std::floor((r_ymax - feat.bbox_ymax) / res_y));
            int r_end   = std::min(n_row - 1, (int)std::floor((r_ymax - feat.bbox_ymin) / res_y));

            if (c_start <= c_end && r_start <= r_end) {
                ScanlineBuffer scanline_buf;

                for (int r = r_start; r <= r_end; ++r) {
                    double cell_yc   = r_ymax - (r + 0.5) * res_y;
                    double row_ymin = r_ymax - (r + 1) * res_y;
                    double row_ymax = r_ymax - r * res_y;

                    for (const auto& part : feat.parts) {
                        if (row_ymax < part.bbox_ymin || row_ymin > part.bbox_ymax) continue;

                        if (!exact) {
                            get_scanline_intervals_stack(part.outer, cell_yc, scanline_buf);
                            if (scanline_buf.n_inter == 0) continue;
                            for (int k = 0; k < scanline_buf.n_inter; ++k) {
                                const auto& inter = scanline_buf.intervals[k];
                                int c_in_start = std::max(c_start, (int)std::floor((inter.first - r_xmin) / res_x));
                                int c_in_end   = std::min(c_end, (int)std::floor((inter.second - r_xmin) / res_x));
                                for (int c = c_in_start; c <= c_in_end; ++c) {
                                    int cell_idx = r * n_col + c;
                                    fr.ext_cell.push_back(cell_idx + 1);
                                    if (return_coverage_area) fr.ext_area.push_back(cell_area);
                                    fr.ext_r.push_back(r + 1);
                                    fr.ext_c.push_back(c + 1);
                                    for (int l = 0; l < n_lyr; ++l) {
                                        fr.ext_vals[l].push_back(val_ptr[l * n_cell_total + cell_idx]);
                                    }
                                }
                            }
                        } else {
                            double outer_min_x, outer_max_x, inner_min_x, inner_max_x;
                            if (!get_row_scanline_extents(part, row_ymin, row_ymax, cell_yc, outer_min_x, outer_max_x, inner_min_x, inner_max_x)) continue;

                            int c_in_start = std::max(c_start, (int)std::floor((outer_min_x - r_xmin) / res_x));
                            int c_in_end   = std::min(c_end, (int)std::floor((outer_max_x - r_xmin) / res_x));

                            int c_int_start = std::max(c_in_start, (int)std::ceil((inner_min_x - r_xmin) / res_x));
                            int c_int_end   = std::min(c_in_end, (int)std::floor((inner_max_x - r_xmin) / res_x) - 1);

                            bool is_pure_interior_row = (row_ymax <= part.bbox_ymax && row_ymin >= part.bbox_ymin);

                            if (!part.holes.empty() || c_int_start > c_int_end || !is_pure_interior_row) {
                                for (int c = c_in_start; c <= c_in_end; ++c) {
                                    double cell_xmin = r_xmin + c * res_x;
                                    double cell_xmax = r_xmin + (c + 1) * res_x;
                                    double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                    if (frac <= 0.0) continue;

                                    int cell_idx = r * n_col + c;
                                    fr.ext_cell.push_back(cell_idx + 1);
                                    if (return_coverage_area) fr.ext_area.push_back(frac * cell_area);
                                    fr.ext_r.push_back(r + 1);
                                    fr.ext_c.push_back(c + 1);
                                    for (int l = 0; l < n_lyr; ++l) {
                                        fr.ext_vals[l].push_back(val_ptr[l * n_cell_total + cell_idx]);
                                    }
                                }
                            } else {
                                // Left boundary
                                for (int c = c_in_start; c < c_int_start && c <= c_in_end; ++c) {
                                    double cell_xmin = r_xmin + c * res_x;
                                    double cell_xmax = r_xmin + (c + 1) * res_x;
                                    double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                    if (frac <= 0.0) continue;

                                    int cell_idx = r * n_col + c;
                                    fr.ext_cell.push_back(cell_idx + 1);
                                    if (return_coverage_area) fr.ext_area.push_back(frac * cell_area);
                                    fr.ext_r.push_back(r + 1);
                                    fr.ext_c.push_back(c + 1);
                                    for (int l = 0; l < n_lyr; ++l) {
                                        fr.ext_vals[l].push_back(val_ptr[l * n_cell_total + cell_idx]);
                                    }
                                }

                                // Pure interior
                                for (int c = c_int_start; c <= c_int_end; ++c) {
                                    int cell_idx = r * n_col + c;
                                    fr.ext_cell.push_back(cell_idx + 1);
                                    if (return_coverage_area) fr.ext_area.push_back(cell_area);
                                    fr.ext_r.push_back(r + 1);
                                    fr.ext_c.push_back(c + 1);
                                    for (int l = 0; l < n_lyr; ++l) {
                                        fr.ext_vals[l].push_back(val_ptr[l * n_cell_total + cell_idx]);
                                    }
                                }

                                // Right boundary
                                for (int c = std::max(c_in_start, c_int_end + 1); c <= c_in_end; ++c) {
                                    double cell_xmin = r_xmin + c * res_x;
                                    double cell_xmax = r_xmin + (c + 1) * res_x;
                                    double frac = compute_exact_cell_2d_fraction(cell_xmin, cell_xmax, row_ymin, row_ymax, res_x, res_y, part, subdiv);
                                    if (frac <= 0.0) continue;

                                    int cell_idx = r * n_col + c;
                                    fr.ext_cell.push_back(cell_idx + 1);
                                    if (return_coverage_area) fr.ext_area.push_back(frac * cell_area);
                                    fr.ext_r.push_back(r + 1);
                                    fr.ext_c.push_back(c + 1);
                                    for (int l = 0; l < n_lyr; ++l) {
                                        fr.ext_vals[l].push_back(val_ptr[l * n_cell_total + cell_idx]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Master thread creates Rcpp List out_list(n_feat) of data.frames directly:
        List out_list(n_feat);
        CharacterVector col_names;
        if (return_coverage_area) {
            col_names = CharacterVector::create("cell", "coverage_area", "row", "col");
        } else {
            col_names = CharacterVector::create("cell", "row", "col");
        }
        for (int l = 0; l < n_lyr; ++l) {
            std::string l_name = "lyr." + std::to_string(l + 1);
            col_names.push_back(l_name);
        }

        for (int i = 0; i < n_feat; ++i) {
            const auto& fr = feat_results[i];
            int n_extracted = fr.ext_cell.size();

            if (n_extracted == 0) {
                List df;
                df.attr("class") = "data.frame";
                df.attr("names") = col_names;
                df.attr("row.names") = IntegerVector::create(0);
                out_list[i] = df;
                continue;
            }

            int n_cols = col_names.size();
            List df(n_cols);

            NumericVector v_cell(fr.ext_cell.begin(), fr.ext_cell.end());
            df[0] = v_cell;

            int next_col = 1;
            if (return_coverage_area) {
                NumericVector v_area(fr.ext_area.begin(), fr.ext_area.end());
                df[next_col++] = v_area;
            }

            NumericVector v_r(fr.ext_r.begin(), fr.ext_r.end());
            df[next_col++] = v_r;

            NumericVector v_c(fr.ext_c.begin(), fr.ext_c.end());
            df[next_col++] = v_c;

            for (int l = 0; l < n_lyr; ++l) {
                NumericVector v_l(fr.ext_vals[l].begin(), fr.ext_vals[l].end());
                df[next_col++] = v_l;
            }

            df.attr("class") = "data.frame";
            df.attr("names") = col_names;
            df.attr("row.names") = IntegerVector::create(NA_INTEGER, -n_extracted);
            out_list[i] = df;
        }
        return out_list;
    }
}
