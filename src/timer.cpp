// High-resolution monotonic clock for benchmarking.
//
// proc.time() on Windows ticks in ~10-16 ms steps, which is too coarse to time
// single inferences that take a few milliseconds. std::chrono::steady_clock is
// monotonic (never jumps with system clock changes) and has sub-microsecond
// resolution on all supported platforms (QueryPerformanceCounter on Windows,
// clock_gettime(CLOCK_MONOTONIC) on Linux/macOS).

#include <Rcpp.h>
#include <chrono>
#include <cstdio>
#include <cstring>

#if defined(_WIN32)
  // PSAPI_VERSION 2 maps GetProcessMemoryInfo to K32GetProcessMemoryInfo in
  // kernel32.dll, so no extra library (-lpsapi) is needed.
  #ifndef PSAPI_VERSION
    #define PSAPI_VERSION 2
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #include <psapi.h>
#elif defined(__APPLE__)
  #include <mach/mach.h>
  #include <sys/resource.h>
#elif defined(__linux__)
  #include <sys/resource.h>
#endif

// [[Rcpp::export]]
double pliman_clock_ms_cpp() {
  using namespace std::chrono;
  return duration<double, std::milli>(steady_clock::now().time_since_epoch()).count();
}

// [[Rcpp::export]]
Rcpp::List pliman_clock_info_cpp() {
  using namespace std::chrono;
  // Smallest observable step between two consecutive reads (empirical resolution).
  double min_step = 1e300;
  for (int i = 0; i < 1000; i++) {
    double a = duration<double, std::milli>(steady_clock::now().time_since_epoch()).count();
    double b = a;
    while (b == a) {
      b = duration<double, std::milli>(steady_clock::now().time_since_epoch()).count();
    }
    if (b - a < min_step) min_step = b - a;
  }
  double tick_ns = 1e9 * static_cast<double>(steady_clock::period::num) /
    static_cast<double>(steady_clock::period::den);
  return Rcpp::List::create(
    Rcpp::Named("clock") = "std::chrono::steady_clock",
    Rcpp::Named("is_steady") = steady_clock::is_steady,
    Rcpp::Named("tick_ns") = tick_ns,
    Rcpp::Named("resolution_ns") = min_step * 1e6
  );
}

// Memory used by the whole R process (resident set size), in bytes.
// Unlike R's own allocation counters, this includes memory allocated by C/C++
// code and libraries such as ONNX Runtime. `peak` is the highest RSS reached by
// the process so far (it never decreases). NA when the platform is not supported.
// [[Rcpp::export]]
Rcpp::NumericVector pliman_process_memory_cpp() {
  double rss = NA_REAL;
  double peak = NA_REAL;
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS pmc;
  std::memset(&pmc, 0, sizeof(pmc));
  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
    rss = static_cast<double>(pmc.WorkingSetSize);
    peak = static_cast<double>(pmc.PeakWorkingSetSize);
  }
#elif defined(__APPLE__)
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS) {
    rss = static_cast<double>(info.resident_size);
    peak = static_cast<double>(info.resident_size_max);
  }
#elif defined(__linux__)
  // VmRSS (current) and VmHWM (peak) are reported in kB in /proc/self/status.
  FILE* f = std::fopen("/proc/self/status", "r");
  if (f != NULL) {
    char line[256];
    while (std::fgets(line, sizeof(line), f) != NULL) {
      double kb = 0;
      if (std::strncmp(line, "VmRSS:", 6) == 0 && std::sscanf(line + 6, "%lf", &kb) == 1) {
        rss = kb * 1024.0;
      } else if (std::strncmp(line, "VmHWM:", 6) == 0 && std::sscanf(line + 6, "%lf", &kb) == 1) {
        peak = kb * 1024.0;
      }
    }
    std::fclose(f);
  }
#endif
  return Rcpp::NumericVector::create(Rcpp::Named("rss") = rss, Rcpp::Named("peak") = peak);
}
