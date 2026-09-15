// onnx_inference.cpp — Direct C++ Dynamic Loader and Inference Engine for ONNX Runtime in pliman
// Zero external R packages, zero Python, zero compile-time library linkage.

#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// Ensure SAL macros are cleanly defined and do not emit unused-result warnings under GCC
#ifdef _Check_return_
#undef _Check_return_
#endif
#define _Check_return_

#ifndef _Frees_ptr_opt_
#define _Frees_ptr_opt_
#endif

#include "onnxruntime_c_api.h"

// [[Rcpp::depends(RcppArmadillo)]]

typedef const OrtApiBase* (ORT_API_CALL *OrtGetApiBaseFn)(void);

static void* g_ort_lib_handle = NULL;
static const OrtApi* g_ort_api = NULL;
static std::string g_current_lib_path = "";

#ifdef _WIN32
#include <dxgi.h>
typedef HRESULT (WINAPI *CreateDXGIFactory1Fn)(REFIID, void**);
typedef OrtStatus* (*OrtAppendDMLFn)(OrtSessionOptions*, int);

static int find_best_gpu_id(std::string* out_name = NULL, double* out_vram = NULL) {
  HMODULE hDxgi = LoadLibraryA("dxgi.dll");
  if (!hDxgi) return 0;

  CreateDXGIFactory1Fn create_dxgi1 = (CreateDXGIFactory1Fn)GetProcAddress(hDxgi, "CreateDXGIFactory1");
  if (!create_dxgi1) {
    FreeLibrary(hDxgi);
    return 0;
  }

  static const GUID IID_IDXGIFactory1_Val = {
    0x770aae78, 0xf26f, 0x4dba, { 0xa8, 0x29, 0x25, 0x3c, 0x83, 0xd1, 0xb3, 0x87 }
  };

  IDXGIFactory1* pFactory = NULL;
  HRESULT hr = create_dxgi1(IID_IDXGIFactory1_Val, (void**)&pFactory);
  if (FAILED(hr) || !pFactory) {
    FreeLibrary(hDxgi);
    return 0;
  }

  int best_id = 0;
  SIZE_T max_vram = 0;
  std::string best_name = "Default GPU";

  IDXGIAdapter1* pAdapter = NULL;
  UINT i = 0;
  while (pFactory->EnumAdapters1(i, &pAdapter) != DXGI_ERROR_NOT_FOUND) {
    DXGI_ADAPTER_DESC1 desc;
    pAdapter->GetDesc1(&desc);

    if (!(desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)) {
      if (desc.DedicatedVideoMemory > max_vram) {
        max_vram = desc.DedicatedVideoMemory;
        best_id = (int)i;

        int sz = WideCharToMultiByte(CP_UTF8, 0, desc.Description, -1, NULL, 0, NULL, NULL);
        std::string n(sz - 1, 0);
        WideCharToMultiByte(CP_UTF8, 0, desc.Description, -1, &n[0], sz, NULL, NULL);
        best_name = n;
      }
    }
    pAdapter->Release();
    ++i;
  }
  pFactory->Release();
  FreeLibrary(hDxgi);

  if (out_name) *out_name = best_name;
  if (out_vram) *out_vram = (double)max_vram / (1024.0 * 1024.0);
  return best_id;
}
#endif

// [[Rcpp::export]]
Rcpp::List pliman_gpu_info_cpp() {
#ifdef _WIN32
  HMODULE hDxgi = LoadLibraryA("dxgi.dll");
  if (!hDxgi) return Rcpp::List::create(Rcpp::Named("available") = false);

  CreateDXGIFactory1Fn create_dxgi1 = (CreateDXGIFactory1Fn)GetProcAddress(hDxgi, "CreateDXGIFactory1");
  if (!create_dxgi1) {
    FreeLibrary(hDxgi);
    return Rcpp::List::create(Rcpp::Named("available") = false);
  }

  static const GUID IID_IDXGIFactory1_Val = {
    0x770aae78, 0xf26f, 0x4dba, { 0xa8, 0x29, 0x25, 0x3c, 0x83, 0xd1, 0xb3, 0x87 }
  };

  IDXGIFactory1* pFactory = NULL;
  HRESULT hr = create_dxgi1(IID_IDXGIFactory1_Val, (void**)&pFactory);
  if (FAILED(hr) || !pFactory) {
    FreeLibrary(hDxgi);
    return Rcpp::List::create(Rcpp::Named("available") = false);
  }

  std::vector<int> ids;
  std::vector<std::string> names;
  std::vector<double> vram_mb;
  std::vector<bool> is_dedicated;
  int best_id = 0;
  SIZE_T max_vram = 0;

  IDXGIAdapter1* pAdapter = NULL;
  UINT i = 0;
  while (pFactory->EnumAdapters1(i, &pAdapter) != DXGI_ERROR_NOT_FOUND) {
    DXGI_ADAPTER_DESC1 desc;
    pAdapter->GetDesc1(&desc);

    if (!(desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)) {
      int sz = WideCharToMultiByte(CP_UTF8, 0, desc.Description, -1, NULL, 0, NULL, NULL);
      std::string n(sz - 1, 0);
      WideCharToMultiByte(CP_UTF8, 0, desc.Description, -1, &n[0], sz, NULL, NULL);

      double vram = (double)desc.DedicatedVideoMemory / (1024.0 * 1024.0);
      ids.push_back((int)i);
      names.push_back(n);
      vram_mb.push_back(vram);
      is_dedicated.push_back(desc.DedicatedVideoMemory > 512 * 1024 * 1024);

      if (desc.DedicatedVideoMemory > max_vram) {
        max_vram = desc.DedicatedVideoMemory;
        best_id = (int)i;
      }
    }
    pAdapter->Release();
    ++i;
  }
  pFactory->Release();
  FreeLibrary(hDxgi);

  return Rcpp::List::create(
    Rcpp::Named("available") = true,
    Rcpp::Named("device_id") = ids,
    Rcpp::Named("name") = names,
    Rcpp::Named("vram_mb") = vram_mb,
    Rcpp::Named("is_dedicated") = is_dedicated,
    Rcpp::Named("default_device_id") = best_id
  );
#else
  return Rcpp::List::create(Rcpp::Named("available") = false);
#endif
}

static void clear_onnx_sessions_internal();

static const OrtApi* get_ort_api(const std::string& lib_path) {
  if (g_ort_api != NULL && g_current_lib_path == lib_path) {
    return g_ort_api;
  }
  if (g_ort_api != NULL && g_current_lib_path != lib_path) {
    clear_onnx_sessions_internal();
    g_ort_api = NULL;
    g_ort_lib_handle = NULL;
  }
  g_current_lib_path = lib_path;

#ifdef _WIN32
  int size_needed = MultiByteToWideChar(CP_UTF8, 0, lib_path.c_str(), (int)lib_path.size(), NULL, 0);
  std::wstring wpath(size_needed, 0);
  MultiByteToWideChar(CP_UTF8, 0, lib_path.c_str(), (int)lib_path.size(), &wpath[0], size_needed);
  HMODULE hMod = LoadLibraryW(wpath.c_str());
  if (!hMod) {
    hMod = LoadLibraryA(lib_path.c_str());
  }
  if (!hMod) {
    Rcpp::stop("Failed to load ONNX Runtime C++ library from: " + lib_path + ". Please run `onnx_install()` to download it.");
  }
  OrtGetApiBaseFn get_base = (OrtGetApiBaseFn)GetProcAddress(hMod, "OrtGetApiBase");
  if (!get_base) {
    Rcpp::stop("Failed to locate entry point OrtGetApiBase in: " + lib_path);
  }
  g_ort_lib_handle = (void*)hMod;
#else
  void* hMod = dlopen(lib_path.c_str(), RTLD_NOW | RTLD_GLOBAL);
  if (!hMod) {
    Rcpp::stop("Failed to load ONNX Runtime library from: " + lib_path + ". Please run `onnx_install()` to download it.");
  }
  OrtGetApiBaseFn get_base = (OrtGetApiBaseFn)dlsym(hMod, "OrtGetApiBase");
  if (!get_base) {
    Rcpp::stop("Failed to locate OrtGetApiBase in: " + lib_path);
  }
  g_ort_lib_handle = hMod;
#endif

  const OrtApiBase* api_base = get_base();
  if (!api_base) {
    Rcpp::stop("OrtGetApiBase returned NULL.");
  }

  g_ort_api = api_base->GetApi(ORT_API_VERSION);
  if (!g_ort_api) {
    for (uint32_t v = ORT_API_VERSION; v >= 1; --v) {
      g_ort_api = api_base->GetApi(v);
      if (g_ort_api) break;
    }
  }
  if (!g_ort_api) {
    Rcpp::stop("Failed to obtain OrtApi structure from ONNX Runtime.");
  }
  return g_ort_api;
}

// Helper: IEEE 754 half-precision float (16-bit) to double conversion
static inline double half_to_double(uint16_t h) {
  uint32_t sign = (h >> 15) & 0x0001;
  uint32_t exp  = (h >> 10) & 0x001f;
  uint32_t mant = h & 0x03ff;

  if (exp == 0) {
    if (mant == 0) {
      return sign ? -0.0 : 0.0;
    } else {
      // Subnormal
      while ((mant & 0x0400) == 0) {
        mant <<= 1;
        exp--;
      }
      exp++;
      mant &= ~0x0400;
      uint32_t f = (sign << 31) | ((exp + (127 - 15)) << 23) | (mant << 13);
      float res = 0.0f;
      std::memcpy(&res, &f, 4);
      return (double)res;
    }
  } else if (exp == 31) {
    uint32_t f = (sign << 31) | (0xff << 23) | (mant << 13);
    float res = 0.0f;
    std::memcpy(&res, &f, 4);
    return (double)res;
  }

  uint32_t f = (sign << 31) | ((exp + (127 - 15)) << 23) | (mant << 13);
  float res = 0.0f;
  std::memcpy(&res, &f, 4);
  return (double)res;
}

struct CachedSession {
  OrtSession* session = NULL;
  OrtSessionOptions* options = NULL;
  OrtEnv* env = NULL;
  OrtMemoryInfo* mem_info = NULL;
};

static std::map<std::string, CachedSession> g_session_cache;

static void clear_onnx_sessions_internal() {
  if (g_ort_api == NULL) return;
  for (auto& pair : g_session_cache) {
    CachedSession& cs = pair.second;
    if (cs.session) g_ort_api->ReleaseSession(cs.session);
    if (cs.mem_info) g_ort_api->ReleaseMemoryInfo(cs.mem_info);
    if (cs.options) g_ort_api->ReleaseSessionOptions(cs.options);
    if (cs.env) g_ort_api->ReleaseEnv(cs.env);
  }
  g_session_cache.clear();
}

static CachedSession get_or_create_cached_session(
    const OrtApi* ort,
    const std::string& model_path,
    int num_threads,
    bool use_gpu = false,
    int device_id = -1
) {
  int actual_device_id = -1;
#ifdef _WIN32
  if (use_gpu) {
    actual_device_id = (device_id >= 0) ? device_id : find_best_gpu_id();
  }
#endif

  std::string key = model_path + "_th_" + std::to_string(num_threads) +
                    (use_gpu ? ("_gpu_" + std::to_string(actual_device_id)) : "_cpu");
  auto it = g_session_cache.find(key);
  if (it != g_session_cache.end()) {
    return it->second;
  }

  CachedSession cs;
  (void)ort->CreateEnv(ORT_LOGGING_LEVEL_ERROR, "pliman_onnx", &cs.env);
  (void)ort->CreateSessionOptions(&cs.options);
  if (num_threads > 0) {
    (void)ort->SetIntraOpNumThreads(cs.options, num_threads);
  }
  (void)ort->SetSessionGraphOptimizationLevel(cs.options, ORT_ENABLE_ALL);
  (void)ort->SetSessionLogSeverityLevel(cs.options, ORT_LOGGING_LEVEL_ERROR);
  (void)ort->EnableCpuMemArena(cs.options);
  (void)ort->EnableMemPattern(cs.options);

#ifdef _WIN32
  if (use_gpu && g_ort_lib_handle != NULL) {
    OrtAppendDMLFn append_dml = (OrtAppendDMLFn)GetProcAddress((HMODULE)g_ort_lib_handle, "OrtSessionOptionsAppendExecutionProvider_DML");
    if (append_dml) {
      OrtStatus* dml_st = append_dml(cs.options, actual_device_id);
      if (dml_st != NULL) {
        ort->ReleaseStatus(dml_st);
      }
    }
  }
#endif

#ifdef _WIN32
  int size_needed = MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.size(), NULL, 0);
  std::wstring wmodel_path(size_needed, 0);
  MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.size(), &wmodel_path[0], size_needed);
  OrtStatus* status = ort->CreateSession(cs.env, wmodel_path.c_str(), cs.options, &cs.session);
#else
  OrtStatus* status = ort->CreateSession(cs.env, model_path.c_str(), cs.options, &cs.session);
#endif

  if (status != NULL || !cs.session) {
    std::string msg = status ? ort->GetErrorMessage(status) : "Unknown error";
    if (status) ort->ReleaseStatus(status);
    if (cs.options) ort->ReleaseSessionOptions(cs.options);
    if (cs.env) ort->ReleaseEnv(cs.env);
    Rcpp::stop("Failed to create ONNX session for " + model_path + ": " + msg);
  }

  (void)ort->CreateCpuMemoryInfo(OrtArenaAllocator, OrtMemTypeDefault, &cs.mem_info);

  g_session_cache[key] = cs;
  return cs;
}

// [[Rcpp::export]]
void clear_onnx_sessions_cpp() {
  clear_onnx_sessions_internal();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix run_onnx_inference_cpp(Rcpp::NumericVector tensor_vec,
                                           Rcpp::IntegerVector tensor_dims,
                                           std::string model_path,
                                           std::string lib_path,
                                           int num_threads = 0,
                                           bool use_gpu = false,
                                           int device_id = -1) {
  const OrtApi* ort = get_ort_api(lib_path);

  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  OrtAllocator* allocator = NULL;
  (void)ort->GetAllocatorWithDefaultOptions(&allocator);

  char* input_name = NULL;
  (void)ort->SessionGetInputName(session, 0, allocator, &input_name);

  char* output_name = NULL;
  (void)ort->SessionGetOutputName(session, 0, allocator, &output_name);

  int64_t input_shape[4] = {
    (int64_t)tensor_dims[0],
    (int64_t)tensor_dims[1],
    (int64_t)tensor_dims[2],
    (int64_t)tensor_dims[3]
  };
  size_t total_elements = tensor_dims[0] * tensor_dims[1] * tensor_dims[2] * tensor_dims[3];

  std::vector<float> input_tensor_values(total_elements);
  for (size_t i = 0; i < total_elements; ++i) {
    input_tensor_values[i] = (float)tensor_vec[i];
  }

  OrtValue* input_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info,
    input_tensor_values.data(),
    total_elements * sizeof(float),
    input_shape,
    4,
    ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
    &input_tensor
  );

  const char* input_names[] = { input_name };
  const char* output_names[] = { output_name };
  OrtValue* output_tensor = NULL;

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtStatus* status = ort->Run(
    session,
    run_options,
    input_names,
    (const OrtValue* const*)&input_tensor,
    1,
    output_names,
    1,
    &output_tensor
  );
  ort->ReleaseRunOptions(run_options);

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseValue(input_tensor);
    allocator->Free(allocator, input_name);
    allocator->Free(allocator, output_name);
    Rcpp::stop("ONNX Run failed: " + msg);
  }

  OrtTensorTypeAndShapeInfo* shape_info = NULL;
  (void)ort->GetTensorTypeAndShape(output_tensor, &shape_info);
  size_t num_dims = 0;
  (void)ort->GetDimensionsCount(shape_info, &num_dims);
  std::vector<int64_t> out_dims(num_dims);
  (void)ort->GetDimensions(shape_info, out_dims.data(), num_dims);

  int out_h = 320;
  int out_w = 320;
  if (num_dims == 4) {
    out_h = (int)out_dims[2];
    out_w = (int)out_dims[3];
  } else if (num_dims == 3) {
    out_h = (int)out_dims[1];
    out_w = (int)out_dims[2];
  } else if (num_dims == 2) {
    out_h = (int)out_dims[0];
    out_w = (int)out_dims[1];
  }

  ONNXTensorElementDataType elem_type = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;
  (void)ort->GetTensorElementType(shape_info, &elem_type);

  Rcpp::NumericMatrix out_mat(out_w, out_h);
  if (elem_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16) {
    const uint16_t* half_data = NULL;
    (void)ort->GetTensorMutableData(output_tensor, (void**)&half_data);
    for (int r = 0; r < out_h; ++r) {
      for (int c = 0; c < out_w; ++c) {
        out_mat(c, r) = half_to_double(half_data[r * out_w + c]);
      }
    }
  } else {
    const float* float_data = NULL;
    (void)ort->GetTensorMutableData(output_tensor, (void**)&float_data);
    for (int r = 0; r < out_h; ++r) {
      for (int c = 0; c < out_w; ++c) {
        out_mat(c, r) = (double)float_data[r * out_w + c];
      }
    }
  }

  ort->ReleaseTensorTypeAndShapeInfo(shape_info);
  ort->ReleaseValue(input_tensor);
  ort->ReleaseValue(output_tensor);
  allocator->Free(allocator, input_name);
  allocator->Free(allocator, output_name);

  return out_mat;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix run_sam2_inference_cpp(Rcpp::NumericVector tensor_vec,
                                           Rcpp::NumericVector points_x,
                                           Rcpp::NumericVector points_y,
                                           Rcpp::IntegerVector point_labels,
                                           std::string encoder_path,
                                           std::string decoder_path,
                                           std::string lib_path,
                                           int num_threads = 0,
                                           bool use_gpu = false,
                                           int device_id = -1) {
  const OrtApi* ort = get_ort_api(lib_path);

  CachedSession cs_enc = get_or_create_cached_session(ort, encoder_path, num_threads, use_gpu, device_id);
  CachedSession cs_dec = get_or_create_cached_session(ort, decoder_path, num_threads, use_gpu, device_id);
  OrtSession* session_enc = cs_enc.session;
  OrtSession* session_dec = cs_dec.session;
  OrtMemoryInfo* memory_info = cs_enc.mem_info;

  // 3. Prepare Image Tensor [1, 3, 1024, 1024]
  int64_t img_shape[4] = {1, 3, 1024, 1024};
  size_t total_img_floats = 1 * 3 * 1024 * 1024;
  std::vector<float> img_values(total_img_floats);
  for (size_t i = 0; i < total_img_floats; ++i) {
    img_values[i] = (float)tensor_vec[i];
  }

  OrtValue* input_img_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info,
    img_values.data(),
    total_img_floats * sizeof(float),
    img_shape,
    4,
    ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
    &input_img_tensor
  );

  // 4. Run Encoder
  const char* enc_in_names[] = { "image" };
  const char* enc_out_names[] = { "high_res_feats_0", "high_res_feats_1", "image_embed" };
  OrtValue* enc_outputs[3] = { NULL, NULL, NULL };

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtStatus* status = ort->Run(
    session_enc,
    run_options,
    enc_in_names,
    (const OrtValue* const*)&input_img_tensor,
    1,
    enc_out_names,
    3,
    enc_outputs
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseValue(input_img_tensor);
    Rcpp::stop("SAM 2.1 Encoder Run failed: " + msg);
  }

  // 5. Prepare Decoder Prompts & Inputs
  // 5. Prepare Decoder Constant Inputs
  int mask_h = 256;
  int mask_w = 256;

  int64_t mask_input_shape[4] = { 1, 1, 256, 256 };
  std::vector<float> mask_input_vals(256 * 256, 0.0f);
  OrtValue* mask_input_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info,
    mask_input_vals.data(),
    mask_input_vals.size() * sizeof(float),
    mask_input_shape,
    4,
    ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
    &mask_input_tensor
  );

  int64_t has_mask_shape[1] = { 1 };
  float has_mask_val = 0.0f;
  OrtValue* has_mask_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info,
    &has_mask_val,
    sizeof(float),
    has_mask_shape,
    1,
    ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
    &has_mask_tensor
  );

  const char* dec_in_names[] = {
    "image_embed",
    "high_res_feats_0",
    "high_res_feats_1",
    "point_coords",
    "point_labels",
    "mask_input",
    "has_mask_input"
  };

  const char* dec_out_names[] = { "masks", "iou_predictions" };

  // out_mat has dimensions [width, height] = [256, 256] matching pliman
  Rcpp::NumericMatrix out_mat(mask_w, mask_h);
  for (int i = 0; i < mask_w; ++i) {
    for (int j = 0; j < mask_h; ++j) {
      out_mat(i, j) = 0.0;
    }
  }

  int num_pts = (int)points_x.size();
  if (num_pts < 1) num_pts = 1;

  for (int p = 0; p < num_pts; ++p) {
    float cur_pt_coords[2] = {
      (float)(p < (int)points_x.size() ? points_x[p] : 512.0),
      (float)(p < (int)points_y.size() ? points_y[p] : 512.0)
    };
    float cur_pt_lbl = (float)(p < (int)point_labels.size() ? point_labels[p] : 1);

    int64_t single_pt_shape[3] = { 1, 1, 2 };
    int64_t single_lbl_shape[2] = { 1, 1 };

    OrtValue* single_pts_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info,
      cur_pt_coords,
      2 * sizeof(float),
      single_pt_shape,
      3,
      ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
      &single_pts_tensor
    );

    OrtValue* single_lbl_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info,
      &cur_pt_lbl,
      sizeof(float),
      single_lbl_shape,
      2,
      ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
      &single_lbl_tensor
    );

    const OrtValue* cur_dec_in_values[] = {
      enc_outputs[2], // image_embed
      enc_outputs[0], // high_res_feats_0
      enc_outputs[1], // high_res_feats_1
      single_pts_tensor,
      single_lbl_tensor,
      mask_input_tensor,
      has_mask_tensor
    };

    OrtValue* cur_dec_outputs[2] = { NULL, NULL };

    OrtStatus* dec_status = ort->Run(
      session_dec,
      run_options,
      dec_in_names,
      cur_dec_in_values,
      7,
      dec_out_names,
      2,
      cur_dec_outputs
    );

    if (dec_status == NULL) {
      OrtTensorTypeAndShapeInfo* iou_shape_info = NULL;
      (void)ort->GetTensorTypeAndShape(cur_dec_outputs[1], &iou_shape_info);
      size_t iou_dims_count = 0;
      (void)ort->GetDimensionsCount(iou_shape_info, &iou_dims_count);
      std::vector<int64_t> iou_dims(iou_dims_count);
      (void)ort->GetDimensions(iou_shape_info, iou_dims.data(), iou_dims_count);
      int num_masks = (iou_dims_count >= 2) ? (int)iou_dims[1] : 1;

      float* iou_data = NULL;
      (void)ort->GetTensorMutableData(cur_dec_outputs[1], (void**)&iou_data);

      float* mask_data = NULL;
      (void)ort->GetTensorMutableData(cur_dec_outputs[0], (void**)&mask_data);

      // Select best mask that is not the macro background
      int best_idx = 0;
      float best_score = -1e9f;
      for (int m = 0; m < num_masks; ++m) {
        size_t offset = (size_t)m * mask_h * mask_w;
        int fg_count = 0;
        for (int j = 0; j < mask_h * mask_w; ++j) {
          if (mask_data[offset + j] > 0.0f) fg_count++;
        }
        float fg_ratio = (float)fg_count / (float)(mask_h * mask_w);
        float score = iou_data[m];
        // If mask covers > 80% of canvas, it's the whole scene/background ambiguity mask
        if (fg_ratio > 0.80f) score -= 3.0f;
        if (score > best_score) {
          best_score = score;
          best_idx = m;
        }
      }

      size_t slice_offset = (size_t)best_idx * mask_h * mask_w;
      for (int r = 0; r < mask_h; ++r) {
        for (int c = 0; c < mask_w; ++c) {
          float logit = mask_data[slice_offset + r * mask_w + c];
          double prob = 1.0 / (1.0 + std::exp(-(double)logit));
          // Store in out_mat(c, r) which is (X, Y)
          if (prob > out_mat(c, r)) {
            out_mat(c, r) = prob;
          }
        }
      }

      ort->ReleaseTensorTypeAndShapeInfo(iou_shape_info);
      ort->ReleaseValue(cur_dec_outputs[0]);
      ort->ReleaseValue(cur_dec_outputs[1]);
    } else {
      ort->ReleaseStatus(dec_status);
    }

    ort->ReleaseValue(single_pts_tensor);
    ort->ReleaseValue(single_lbl_tensor);
  }

  ort->ReleaseRunOptions(run_options);
  ort->ReleaseValue(input_img_tensor);
  ort->ReleaseValue(mask_input_tensor);
  ort->ReleaseValue(has_mask_tensor);
  for (int i = 0; i < 3; ++i) ort->ReleaseValue(enc_outputs[i]);

  return out_mat;
}

struct DinoBox {
  double x1, y1, x2, y2;
  double score;
  int token_idx;
};

// [[Rcpp::export]]
Rcpp::List run_grounding_dino_cpp(Rcpp::NumericVector pixel_values,
                                 Rcpp::IntegerVector input_ids,
                                 Rcpp::IntegerVector token_type_ids,
                                 Rcpp::IntegerVector attention_mask,
                                 std::string model_path,
                                 std::string lib_path,
                                 double box_threshold = 0.25,
                                 double text_threshold = 0.25,
                                 double iou_threshold = 0.5,
                                 int num_threads = 0,
                                 bool use_gpu = false,
                                 int device_id = -1) {
  const OrtApi* ort = get_ort_api(lib_path);

  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  // 1. Prepare pixel_values [1, 3, 800, 800]
  int64_t img_shape[4] = { 1, 3, 800, 800 };
  size_t total_pix = 1 * 3 * 800 * 800;
  std::vector<float> img_values(total_pix);
  for (size_t i = 0; i < total_pix; ++i) {
    img_values[i] = (float)pixel_values[i];
  }
  OrtValue* pixel_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, img_values.data(), total_pix * sizeof(float),
    img_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &pixel_tensor
  );

  // 2. Prepare input_ids [1, seq_len]
  int seq_len = input_ids.size();
  int64_t seq_shape[2] = { 1, (int64_t)seq_len };
  std::vector<int64_t> ids_vec(seq_len);
  for (int i = 0; i < seq_len; ++i) ids_vec[i] = (int64_t)input_ids[i];
  OrtValue* ids_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, ids_vec.data(), seq_len * sizeof(int64_t),
    seq_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64, &ids_tensor
  );

  // 3. Prepare token_type_ids [1, seq_len]
  std::vector<int64_t> type_vec(seq_len);
  for (int i = 0; i < seq_len; ++i) type_vec[i] = (int64_t)token_type_ids[i];
  OrtValue* type_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, type_vec.data(), seq_len * sizeof(int64_t),
    seq_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64, &type_tensor
  );

  // 4. Prepare attention_mask [1, seq_len]
  std::vector<int64_t> att_vec(seq_len);
  for (int i = 0; i < seq_len; ++i) att_vec[i] = (int64_t)attention_mask[i];
  OrtValue* att_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, att_vec.data(), seq_len * sizeof(int64_t),
    seq_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64, &att_tensor
  );

  // 5. Prepare pixel_mask [1, 800, 800]
  int64_t mask_shape[3] = { 1, 800, 800 };
  size_t total_mask_pix = 800 * 800;
  std::vector<int64_t> mask_vec(total_mask_pix, 1LL);
  OrtValue* mask_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, mask_vec.data(), total_mask_pix * sizeof(int64_t),
    mask_shape, 3, ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64, &mask_tensor
  );

  const char* in_names[] = {
    "pixel_values", "input_ids", "token_type_ids", "attention_mask", "pixel_mask"
  };
  const OrtValue* in_tensors[] = {
    pixel_tensor, ids_tensor, type_tensor, att_tensor, mask_tensor
  };

  const char* out_names[] = { "logits", "pred_boxes" };
  OrtValue* out_tensors[2] = { NULL, NULL };

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtStatus* status = ort->Run(
    session, run_options, in_names, in_tensors, 5, out_names, 2, out_tensors
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(pixel_tensor);
    ort->ReleaseValue(ids_tensor);
    ort->ReleaseValue(type_tensor);
    ort->ReleaseValue(att_tensor);
    ort->ReleaseValue(mask_tensor);
    Rcpp::stop("Grounding DINO inference failed: " + msg);
  }

  float* logits_data = NULL;
  (void)ort->GetTensorMutableData(out_tensors[0], (void**)&logits_data);

  float* boxes_data = NULL;
  (void)ort->GetTensorMutableData(out_tensors[1], (void**)&boxes_data);

  // Post-processing: 900 queries, 256 token logits
  std::vector<DinoBox> candidates;
  for (int q = 0; q < 900; ++q) {
    float cx = boxes_data[q * 4 + 0];
    float cy = boxes_data[q * 4 + 1];
    float w  = boxes_data[q * 4 + 2];
    float h  = boxes_data[q * 4 + 3];

    double max_prob = -1.0;
    int best_t = -1;
    for (int t = 1; t < seq_len - 1; ++t) {
      float logit = logits_data[q * 256 + t];
      double prob = 1.0 / (1.0 + std::exp(-(double)logit));
      if (prob > max_prob) {
        max_prob = prob;
        best_t = t;
      }
    }

    if (max_prob >= box_threshold) {
      DinoBox b;
      b.x1 = std::max(0.0, (double)(cx - w / 2.0));
      b.y1 = std::max(0.0, (double)(cy - h / 2.0));
      b.x2 = std::min(1.0, (double)(cx + w / 2.0));
      b.y2 = std::min(1.0, (double)(cy + h / 2.0));
      b.score = max_prob;
      b.token_idx = best_t;
      candidates.push_back(b);
    }
  }

  std::sort(candidates.begin(), candidates.end(), [](const DinoBox& a, const DinoBox& b) {
    return a.score > b.score;
  });

  // Apply NMS
  std::vector<DinoBox> kept;
  std::vector<bool> suppressed(candidates.size(), false);
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (suppressed[i]) continue;
    kept.push_back(candidates[i]);
    for (size_t j = i + 1; j < candidates.size(); ++j) {
      if (suppressed[j]) continue;
      double ix1 = std::max(candidates[i].x1, candidates[j].x1);
      double iy1 = std::max(candidates[i].y1, candidates[j].y1);
      double ix2 = std::min(candidates[i].x2, candidates[j].x2);
      double iy2 = std::min(candidates[i].y2, candidates[j].y2);
      double iw = std::max(0.0, ix2 - ix1);
      double ih = std::max(0.0, iy2 - iy1);
      double inter = iw * ih;
      double area1 = (candidates[i].x2 - candidates[i].x1) * (candidates[i].y2 - candidates[i].y1);
      double area2 = (candidates[j].x2 - candidates[j].x1) * (candidates[j].y2 - candidates[j].y1);
      double union_a = area1 + area2 - inter;
      double iou = (union_a > 0.0) ? (inter / union_a) : 0.0;
      if (iou > iou_threshold) {
        suppressed[j] = true;
      }
    }
  }

  int num_kept = (int)kept.size();
  Rcpp::NumericMatrix res_boxes(num_kept, 4);
  Rcpp::NumericVector res_scores(num_kept);
  Rcpp::IntegerVector res_tokens(num_kept);

  for (int i = 0; i < num_kept; ++i) {
    res_boxes(i, 0) = kept[i].x1;
    res_boxes(i, 1) = kept[i].y1;
    res_boxes(i, 2) = kept[i].x2;
    res_boxes(i, 3) = kept[i].y2;
    res_scores[i]   = kept[i].score;
    res_tokens[i]   = kept[i].token_idx;
  }

  // Cleanup
  ort->ReleaseRunOptions(run_options);
  ort->ReleaseValue(pixel_tensor);
  ort->ReleaseValue(ids_tensor);
  ort->ReleaseValue(type_tensor);
  ort->ReleaseValue(att_tensor);
  ort->ReleaseValue(mask_tensor);
  ort->ReleaseValue(out_tensors[0]);
  ort->ReleaseValue(out_tensors[1]);

  return Rcpp::List::create(
    Rcpp::Named("boxes") = res_boxes,
    Rcpp::Named("scores") = res_scores,
    Rcpp::Named("token_indices") = res_tokens
  );
}

// [[Rcpp::export]]
Rcpp::List run_sam2_instances_cpp(Rcpp::NumericVector tensor_vec,
                                  Rcpp::NumericMatrix boxes,
                                  double orig_w,
                                  double orig_h,
                                  std::string encoder_path,
                                  std::string decoder_path,
                                  std::string lib_path,
                                  int num_threads = 0,
                                  bool use_gpu = false,
                                  int device_id = -1) {
  const OrtApi* ort = get_ort_api(lib_path);

  CachedSession cs_enc = get_or_create_cached_session(ort, encoder_path, num_threads, use_gpu, device_id);
  CachedSession cs_dec = get_or_create_cached_session(ort, decoder_path, num_threads, use_gpu, device_id);
  OrtSession* session_enc = cs_enc.session;
  OrtSession* session_dec = cs_dec.session;
  OrtMemoryInfo* memory_info = cs_enc.mem_info;

  // 1. Run Encoder ONCE on [1, 3, 1024, 1024]
  int64_t img_shape[4] = {1, 3, 1024, 1024};
  size_t total_img_floats = 1 * 3 * 1024 * 1024;
  std::vector<float> img_values(total_img_floats);
  for (size_t i = 0; i < total_img_floats; ++i) {
    img_values[i] = (float)tensor_vec[i];
  }

  OrtValue* input_img_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, img_values.data(), total_img_floats * sizeof(float),
    img_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &input_img_tensor
  );

  const char* enc_in_names[] = { "image" };
  const char* enc_out_names[] = { "high_res_feats_0", "high_res_feats_1", "image_embed" };
  OrtValue* enc_outputs[3] = { NULL, NULL, NULL };

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtStatus* status = ort->Run(
    session_enc, run_options, enc_in_names, (const OrtValue* const*)&input_img_tensor,
    1, enc_out_names, 3, enc_outputs
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(input_img_tensor);
    Rcpp::stop("SAM 2.1 Encoder Run failed: " + msg);
  }

  // 2. Decoder Common Setup
  int mask_h = 256;
  int mask_w = 256;

  int64_t mask_input_shape[4] = { 1, 1, 256, 256 };
  std::vector<float> mask_input_vals(256 * 256, 0.0f);
  OrtValue* mask_input_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, mask_input_vals.data(), mask_input_vals.size() * sizeof(float),
    mask_input_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &mask_input_tensor
  );

  int64_t has_mask_shape[1] = { 1 };
  float has_mask_val = 0.0f;
  OrtValue* has_mask_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, &has_mask_val, sizeof(float),
    has_mask_shape, 1, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &has_mask_tensor
  );

  const char* dec_in_names[] = {
    "image_embed", "high_res_feats_0", "high_res_feats_1",
    "point_coords", "point_labels", "mask_input", "has_mask_input"
  };
  const char* dec_out_names[] = { "masks", "iou_predictions" };

  int num_boxes = boxes.nrow();
  Rcpp::List res_masks(num_boxes);

  for (int b = 0; b < num_boxes; ++b) {
    float x1 = (float)(boxes(b, 0) / orig_w * 1024.0);
    float y1 = (float)(boxes(b, 1) / orig_h * 1024.0);
    float x2 = (float)(boxes(b, 2) / orig_w * 1024.0);
    float y2 = (float)(boxes(b, 3) / orig_h * 1024.0);

    float box_coords[4] = { x1, y1, x2, y2 };
    float box_labels[2] = { 2.0f, 3.0f };

    int64_t pts_shape[3] = { 1, 2, 2 };
    int64_t lbls_shape[2] = { 1, 2 };

    OrtValue* pts_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, box_coords, 4 * sizeof(float),
      pts_shape, 3, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &pts_tensor
    );

    OrtValue* lbl_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, box_labels, 2 * sizeof(float),
      lbls_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &lbl_tensor
    );

    const OrtValue* cur_dec_in_values[] = {
      enc_outputs[2], enc_outputs[0], enc_outputs[1],
      pts_tensor, lbl_tensor, mask_input_tensor, has_mask_tensor
    };

    OrtValue* cur_dec_outputs[2] = { NULL, NULL };
    OrtStatus* dec_status = ort->Run(
      session_dec, run_options, dec_in_names, cur_dec_in_values, 7, dec_out_names, 2, cur_dec_outputs
    );

    Rcpp::NumericMatrix inst_mat(mask_w, mask_h);
    for (int i = 0; i < mask_w; ++i) {
      for (int j = 0; j < mask_h; ++j) {
        inst_mat(i, j) = 0.0;
      }
    }

    if (dec_status == NULL) {
      OrtTensorTypeAndShapeInfo* iou_shape_info = NULL;
      (void)ort->GetTensorTypeAndShape(cur_dec_outputs[1], &iou_shape_info);
      size_t iou_dims_count = 0;
      (void)ort->GetDimensionsCount(iou_shape_info, &iou_dims_count);
      std::vector<int64_t> iou_dims(iou_dims_count);
      (void)ort->GetDimensions(iou_shape_info, iou_dims.data(), iou_dims_count);
      int num_masks = (iou_dims_count >= 2) ? (int)iou_dims[1] : 1;

      float* iou_data = NULL;
      (void)ort->GetTensorMutableData(cur_dec_outputs[1], (void**)&iou_data);

      float* mask_data = NULL;
      (void)ort->GetTensorMutableData(cur_dec_outputs[0], (void**)&mask_data);

      int best_idx = 0;
      float best_score = -1e9f;
      for (int m = 0; m < num_masks; ++m) {
        if (iou_data[m] > best_score) {
          best_score = iou_data[m];
          best_idx = m;
        }
      }

      size_t slice_offset = (size_t)best_idx * mask_h * mask_w;
      for (int r = 0; r < mask_h; ++r) {
        for (int c = 0; c < mask_w; ++c) {
          float logit = mask_data[slice_offset + r * mask_w + c];
          double prob = 1.0 / (1.0 + std::exp(-(double)logit));
          inst_mat(c, r) = prob;
        }
      }

      ort->ReleaseTensorTypeAndShapeInfo(iou_shape_info);
      ort->ReleaseValue(cur_dec_outputs[0]);
      ort->ReleaseValue(cur_dec_outputs[1]);
    } else {
      ort->ReleaseStatus(dec_status);
    }

    ort->ReleaseValue(pts_tensor);
    ort->ReleaseValue(lbl_tensor);

    res_masks[b] = inst_mat;
  }

  // Cleanup
  ort->ReleaseRunOptions(run_options);
  ort->ReleaseValue(input_img_tensor);
  ort->ReleaseValue(mask_input_tensor);
  ort->ReleaseValue(has_mask_tensor);
  for (int i = 0; i < 3; ++i) ort->ReleaseValue(enc_outputs[i]);

  return res_masks;
}

struct ZnccKernelEntry {
  int dx;
  int dy;
  float w;
  float t_zero[3];
};

struct DenseExemplarKernel {
  std::vector<ZnccKernelEntry> entries;
  float tmpl_std;
  float sum_w;
  float inv_sum_w;
  int radius;
};

// [[Rcpp::export]]
Rcpp::List run_sam2_persam_cpp(Rcpp::NumericVector tensor_vec,
                               Rcpp::NumericVector exemplar_x,
                               Rcpp::NumericVector exemplar_y,
                               double orig_w,
                               double orig_h,
                               std::string encoder_path,
                               std::string decoder_path,
                               std::string lib_path,
                               double sim_threshold = 0.5,
                               double min_dist = 16.0,
                               double iou_threshold = 0.5,
                               int max_objects = 200,
                               int feat_res = 256,
                               int num_threads = 0,
                               bool use_gpu = false,
                               int device_id = -1) {
  const OrtApi* ort = get_ort_api(lib_path);

  CachedSession cs_enc = get_or_create_cached_session(ort, encoder_path, num_threads, use_gpu, device_id);
  CachedSession cs_dec = get_or_create_cached_session(ort, decoder_path, num_threads, use_gpu, device_id);
  OrtSession* session_enc = cs_enc.session;
  OrtSession* session_dec = cs_dec.session;
  OrtMemoryInfo* memory_info = cs_enc.mem_info;

  // 1. Run SAM 2.1 Encoder ONCE on [1, 3, 1024, 1024]
  int64_t img_shape[4] = {1, 3, 1024, 1024};
  size_t total_img_floats = 1 * 3 * 1024 * 1024;
  std::vector<float> img_values(total_img_floats);
  for (size_t i = 0; i < total_img_floats; ++i) {
    img_values[i] = (float)tensor_vec[i];
  }

  OrtValue* input_img_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, img_values.data(), total_img_floats * sizeof(float),
    img_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &input_img_tensor
  );

  const char* enc_in_names[] = { "image" };
  const char* enc_out_names[] = { "high_res_feats_0", "high_res_feats_1", "image_embed" };
  OrtValue* enc_outputs[3] = { NULL, NULL, NULL };

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtStatus* status = ort->Run(
    session_enc, run_options, enc_in_names, (const OrtValue* const*)&input_img_tensor,
    1, enc_out_names, 3, enc_outputs
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(input_img_tensor);
    Rcpp::stop("SAM 2.1 Encoder Run failed: " + msg);
  }

  // 2. Decoder Common Setup
  int mask_h = 256;
  int mask_w = 256;

  int64_t mask_input_shape[4] = { 1, 1, 256, 256 };
  std::vector<float> mask_input_vals(256 * 256, 0.0f);
  OrtValue* mask_input_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, mask_input_vals.data(), mask_input_vals.size() * sizeof(float),
    mask_input_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &mask_input_tensor
  );

  int64_t has_mask_shape[1] = { 1 };
  float has_mask_val = 0.0f;
  OrtValue* has_mask_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, &has_mask_val, sizeof(float),
    has_mask_shape, 1, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &has_mask_tensor
  );

  const char* dec_in_names[] = {
    "image_embed", "high_res_feats_0", "high_res_feats_1",
    "point_coords", "point_labels", "mask_input", "has_mask_input"
  };
  const char* dec_out_names[] = { "masks", "iou_predictions" };

  float* embed_data = NULL;
  (void)ort->GetTensorMutableData(enc_outputs[2], (void**)&embed_data);

  float* hr0_data = NULL;
  (void)ort->GetTensorMutableData(enc_outputs[0], (void**)&hr0_data);

  // 3. Extract Multi-Prototype Feature Vectors (both 256-dim semantic and 32-dim high-res) and Dense Exemplar Kernels
  std::vector<std::vector<float>> exemplar_prototypes;
  std::vector<std::vector<float>> exemplar_prototypes_hr0;
  std::vector<float> exemplar_areas;
  std::vector<DenseExemplarKernel> exemplar_dense_kernels;
  int num_exemplars = (int)exemplar_x.size();
  if (num_exemplars < 1) {
    exemplar_x = Rcpp::NumericVector::create(orig_w / 2.0);
    exemplar_y = Rcpp::NumericVector::create(orig_h / 2.0);
    num_exemplars = 1;
  }

  for (int ex = 0; ex < num_exemplars; ++ex) {
    float px = (float)(exemplar_x[ex] / orig_w * 1024.0);
    float py = (float)(exemplar_y[ex] / orig_h * 1024.0);
    px = std::max(0.0f, std::min(1023.0f, px));
    py = std::max(0.0f, std::min(1023.0f, py));

    float pt_coords[2] = { px, py };
    float pt_label[1] = { 1.0f };

    int64_t pt_shape[3] = { 1, 1, 2 };
    int64_t lbl_shape[2] = { 1, 1 };

    OrtValue* ex_pts_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, pt_coords, 2 * sizeof(float),
      pt_shape, 3, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &ex_pts_tensor
    );

    OrtValue* ex_lbl_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, pt_label, 1 * sizeof(float),
      lbl_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &ex_lbl_tensor
    );

    const OrtValue* ex_dec_in_values[] = {
      enc_outputs[2], enc_outputs[0], enc_outputs[1],
      ex_pts_tensor, ex_lbl_tensor, mask_input_tensor, has_mask_tensor
    };

    OrtValue* ex_dec_outputs[2] = { NULL, NULL };
    OrtStatus* ex_dec_status = ort->Run(
      session_dec, run_options, dec_in_names, ex_dec_in_values, 7, dec_out_names, 2, ex_dec_outputs
    );

    std::vector<float> proto_v(256, 0.0f);
    std::vector<float> proto_hr0(32, 0.0f);

    if (ex_dec_status == NULL) {
      OrtTensorTypeAndShapeInfo* iou_shape_info = NULL;
      (void)ort->GetTensorTypeAndShape(ex_dec_outputs[1], &iou_shape_info);
      size_t iou_dims_count = 0;
      (void)ort->GetDimensionsCount(iou_shape_info, &iou_dims_count);
      std::vector<int64_t> iou_dims(iou_dims_count);
      (void)ort->GetDimensions(iou_shape_info, iou_dims.data(), iou_dims_count);
      int num_masks = (iou_dims_count >= 2) ? (int)iou_dims[1] : 1;

      float* iou_data = NULL;
      (void)ort->GetTensorMutableData(ex_dec_outputs[1], (void**)&iou_data);

      float* mask_data = NULL;
      (void)ort->GetTensorMutableData(ex_dec_outputs[0], (void**)&mask_data);

      int best_idx = 0;
      float best_score = -1e9f;
      for (int m = 0; m < num_masks; ++m) {
        if (iou_data[m] > best_score) {
          best_score = iou_data[m];
          best_idx = m;
        }
      }

      size_t slice_offset = (size_t)best_idx * mask_h * mask_w;
      float total_weight = 0.0f;
      float total_hr0_weight = 0.0f;

      // Extract 32-D high-resolution prototype directly on 256x256 mask
      for (int r = 0; r < 256; ++r) {
        for (int k = 0; k < 256; ++k) {
          if (mask_data[slice_offset + (size_t)r * 256 + k] > 0.0f) {
            total_hr0_weight += 1.0f;
            size_t cell_offset = (size_t)r * 256 + k;
            for (int c = 0; c < 32; ++c) {
              proto_hr0[c] += hr0_data[(size_t)c * 65536 + cell_offset];
            }
          }
        }
      }
      if (total_hr0_weight > 0.0f) {
        for (int c = 0; c < 32; ++c) proto_hr0[c] /= total_hr0_weight;
        exemplar_areas.push_back(total_hr0_weight);
      } else {
        int k256 = std::min(255, std::max(0, (int)(px / 4.0f)));
        int r256 = std::min(255, std::max(0, (int)(py / 4.0f)));
        size_t cell_offset = (size_t)r256 * 256 + k256;
        for (int c = 0; c < 32; ++c) proto_hr0[c] = hr0_data[(size_t)c * 65536 + cell_offset];
        exemplar_areas.push_back(50.0f);
      }

      // Extract 256-D semantic prototype from 64x64 embed
      for (int r = 0; r < 64; ++r) {
        for (int k = 0; k < 64; ++k) {
          float w = 0.0f;
          for (int dr = 0; dr < 4; ++dr) {
            for (int dc = 0; dc < 4; ++dc) {
              if (mask_data[slice_offset + (r * 4 + dr) * 256 + (k * 4 + dc)] > 0.0f) {
                w += 1.0f;
              }
            }
          }
          w /= 16.0f;
          if (w > 0.1f) {
            total_weight += w;
            size_t cell_offset = (size_t)r * 64 + k;
            for (int c = 0; c < 256; ++c) {
              proto_v[c] += w * embed_data[(size_t)c * 4096 + cell_offset];
            }
          }
        }
      }

      if (total_weight > 0.0f) {
        for (int c = 0; c < 256; ++c) {
          proto_v[c] /= total_weight;
        }
      } else {
        int k_pt = std::min(63, std::max(0, (int)(px / 16.0f)));
        int r_pt = std::min(63, std::max(0, (int)(py / 16.0f)));
        size_t cell_offset = (size_t)r_pt * 64 + k_pt;
        for (int c = 0; c < 256; ++c) {
          proto_v[c] = embed_data[(size_t)c * 4096 + cell_offset];
        }
      }

      // Extract Dense Convolutional Exemplar Kernel from img_values
      float sum_x = 0.0f, sum_y = 0.0f;
      float mask_cnt = 0.0f;
      for (int r = 0; r < 256; ++r) {
        for (int k = 0; k < 256; ++k) {
          if (mask_data[slice_offset + (size_t)r * 256 + k] > 0.0f) {
            sum_x += (float)(k * 4 + 2);
            sum_y += (float)(r * 4 + 2);
            mask_cnt += 1.0f;
          }
        }
      }
      float cen_x = (mask_cnt > 0.0f) ? (sum_x / mask_cnt) : px;
      float cen_y = (mask_cnt > 0.0f) ? (sum_y / mask_cnt) : py;
      float area_1024 = (mask_cnt > 0.0f) ? (mask_cnt * 16.0f) : 256.0f;
      float r_eq = std::sqrt(area_1024 / 3.14159265f);
      int rad = std::min(48, std::max(6, (int)std::round(r_eq)));

      float r_sq = (float)(rad * rad);
      float k_sum_w = 0.0f;
      float k_sum_tc[3] = {0.0f, 0.0f, 0.0f};

      struct RawKernelEntry {
        int dx, dy;
        float w;
        float raw_t[3];
      };
      std::vector<RawKernelEntry> raw_entries;
      size_t plane_size = 1024 * 1024;

      for (int dy = -rad; dy <= rad; ++dy) {
        int sy = std::max(0, std::min(1023, (int)std::round(cen_y) + dy));
        for (int dx = -rad; dx <= rad; ++dx) {
          float d_sq = (float)(dx * dx + dy * dy);
          if (d_sq <= r_sq) {
            int sx = std::max(0, std::min(1023, (int)std::round(cen_x) + dx));
            int r256 = sy / 4;
            int k256 = sx / 4;
            bool in_mask = (mask_cnt > 0.0f) && (mask_data[slice_offset + (size_t)r256 * 256 + k256] > 0.0f);
            float w = in_mask ? (1.0f - 0.5f * (d_sq / r_sq)) : 0.0f;
            if (w > 0.0f) {
              RawKernelEntry re;
              re.dx = dx;
              re.dy = dy;
              re.w = w;
              for (int c = 0; c < 3; ++c) {
                re.raw_t[c] = img_values[(size_t)c * plane_size + (size_t)sy * 1024 + sx];
                k_sum_tc[c] += w * re.raw_t[c];
              }
              k_sum_w += w;
              raw_entries.push_back(re);
            }
          }
        }
      }

      // Fallback to circular disk if mask was empty or too small
      if (raw_entries.size() < 10) {
        raw_entries.clear();
        k_sum_w = 0.0f;
        k_sum_tc[0] = k_sum_tc[1] = k_sum_tc[2] = 0.0f;
        for (int dy = -rad; dy <= rad; ++dy) {
          int sy = std::max(0, std::min(1023, (int)std::round(cen_y) + dy));
          for (int dx = -rad; dx <= rad; ++dx) {
            float d_sq = (float)(dx * dx + dy * dy);
            if (d_sq <= r_sq) {
              int sx = std::max(0, std::min(1023, (int)std::round(cen_x) + dx));
              float w = 1.0f - 0.5f * (d_sq / r_sq);
              RawKernelEntry re;
              re.dx = dx;
              re.dy = dy;
              re.w = w;
              for (int c = 0; c < 3; ++c) {
                re.raw_t[c] = img_values[(size_t)c * plane_size + (size_t)sy * 1024 + sx];
                k_sum_tc[c] += w * re.raw_t[c];
              }
              k_sum_w += w;
              raw_entries.push_back(re);
            }
          }
        }
      }

      if (k_sum_w < 1e-5f) k_sum_w = 1.0f;
      float mu_tc[3] = { k_sum_tc[0] / k_sum_w, k_sum_tc[1] / k_sum_w, k_sum_tc[2] / k_sum_w };

      float tmpl_std_sq = 0.0f;
      DenseExemplarKernel d_kernel;
      d_kernel.radius = rad;
      d_kernel.sum_w = k_sum_w;
      d_kernel.inv_sum_w = 1.0f / k_sum_w;

      for (const auto& re : raw_entries) {
        ZnccKernelEntry ke;
        ke.dx = re.dx;
        ke.dy = re.dy;
        ke.w = re.w;
        for (int c = 0; c < 3; ++c) {
          ke.t_zero[c] = re.raw_t[c] - mu_tc[c];
          tmpl_std_sq += re.w * ke.t_zero[c] * ke.t_zero[c];
        }
        d_kernel.entries.push_back(ke);
      }
      d_kernel.tmpl_std = std::sqrt(tmpl_std_sq);
      if (d_kernel.tmpl_std < 1e-6f) d_kernel.tmpl_std = 1.0f;
      exemplar_dense_kernels.push_back(d_kernel);

      ort->ReleaseTensorTypeAndShapeInfo(iou_shape_info);
      ort->ReleaseValue(ex_dec_outputs[0]);
      ort->ReleaseValue(ex_dec_outputs[1]);
    } else {
      ort->ReleaseStatus(ex_dec_status);
      int k_pt = std::min(63, std::max(0, (int)(px / 16.0f)));
      int r_pt = std::min(63, std::max(0, (int)(py / 16.0f)));
      size_t cell_offset = (size_t)r_pt * 64 + k_pt;
      for (int c = 0; c < 256; ++c) {
        proto_v[c] = embed_data[(size_t)c * 4096 + cell_offset];
      }

      int k256 = std::min(255, std::max(0, (int)(px / 4.0f)));
      int r256 = std::min(255, std::max(0, (int)(py / 4.0f)));
      size_t cell_offset_256 = (size_t)r256 * 256 + k256;
      for (int c = 0; c < 32; ++c) {
        proto_hr0[c] = hr0_data[(size_t)c * 65536 + cell_offset_256];
      }

      int rad = 16;
      float r_sq = (float)(rad * rad);
      float k_sum_w = 0.0f;
      float k_sum_tc[3] = {0.0f, 0.0f, 0.0f};
      struct RawKernelEntry { int dx, dy; float w; float raw_t[3]; };
      std::vector<RawKernelEntry> raw_entries;
      size_t plane_size = 1024 * 1024;
      for (int dy = -rad; dy <= rad; ++dy) {
        int sy = std::max(0, std::min(1023, (int)std::round(py) + dy));
        for (int dx = -rad; dx <= rad; ++dx) {
          float d_sq = (float)(dx * dx + dy * dy);
          if (d_sq <= r_sq) {
            int sx = std::max(0, std::min(1023, (int)std::round(px) + dx));
            float w = 1.0f - 0.5f * (d_sq / r_sq);
            RawKernelEntry re;
            re.dx = dx; re.dy = dy; re.w = w;
            for (int c = 0; c < 3; ++c) {
              re.raw_t[c] = img_values[(size_t)c * plane_size + (size_t)sy * 1024 + sx];
              k_sum_tc[c] += w * re.raw_t[c];
            }
            k_sum_w += w;
            raw_entries.push_back(re);
          }
        }
      }
      if (k_sum_w < 1e-5f) k_sum_w = 1.0f;
      float mu_tc[3] = { k_sum_tc[0] / k_sum_w, k_sum_tc[1] / k_sum_w, k_sum_tc[2] / k_sum_w };
      float tmpl_std_sq = 0.0f;
      DenseExemplarKernel d_kernel;
      d_kernel.radius = rad;
      d_kernel.sum_w = k_sum_w;
      d_kernel.inv_sum_w = 1.0f / k_sum_w;
      for (const auto& re : raw_entries) {
        ZnccKernelEntry ke;
        ke.dx = re.dx; ke.dy = re.dy; ke.w = re.w;
        for (int c = 0; c < 3; ++c) {
          ke.t_zero[c] = re.raw_t[c] - mu_tc[c];
          tmpl_std_sq += re.w * ke.t_zero[c] * ke.t_zero[c];
        }
        d_kernel.entries.push_back(ke);
      }
      d_kernel.tmpl_std = std::sqrt(tmpl_std_sq);
      if (d_kernel.tmpl_std < 1e-6f) d_kernel.tmpl_std = 1.0f;
      exemplar_dense_kernels.push_back(d_kernel);
    }

    ort->ReleaseValue(ex_pts_tensor);
    ort->ReleaseValue(ex_lbl_tensor);

    // Normalize semantic prototype to unit norm
    float proto_norm_sq = 0.0f;
    for (int c = 0; c < 256; ++c) {
      proto_norm_sq += proto_v[c] * proto_v[c];
    }
    float proto_norm = std::sqrt(proto_norm_sq);
    if (proto_norm < 1e-8f) proto_norm = 1.0f;
    for (int c = 0; c < 256; ++c) {
      proto_v[c] /= proto_norm;
    }
    exemplar_prototypes.push_back(proto_v);

    // Normalize high-res prototype to unit norm
    float hr0_norm_sq = 0.0f;
    for (int c = 0; c < 32; ++c) {
      hr0_norm_sq += proto_hr0[c] * proto_hr0[c];
    }
    float hr0_norm = std::sqrt(hr0_norm_sq);
    if (hr0_norm < 1e-8f) hr0_norm = 1.0f;
    for (int c = 0; c < 32; ++c) {
      proto_hr0[c] /= hr0_norm;
    }
    exemplar_prototypes_hr0.push_back(proto_hr0);
  }

  float avg_ex_area = 100.0f;
  if (!exemplar_areas.empty()) {
    float sum_area = 0.0f;
    int cnt_area = 0;
    for (float a : exemplar_areas) {
      if (a > 5.0f) {
        sum_area += a;
        cnt_area++;
      }
    }
    if (cnt_area > 0) avg_ex_area = sum_area / cnt_area;
  }

  // 4. Compute Cosine Similarity Heatmap with Multi-Prototype Max-Pooling (64x64 base):
  std::vector<float> sim_grid(64 * 64, 0.0f);
  int num_protos = (int)exemplar_prototypes.size();

  for (int r = 0; r < 64; ++r) {
    for (int k = 0; k < 64; ++k) {
      size_t cell_offset = (size_t)r * 64 + k;

      float cell_norm_sq = 0.0f;
      for (int c = 0; c < 256; ++c) {
        float val = embed_data[(size_t)c * 4096 + cell_offset];
        cell_norm_sq += val * val;
      }
      float cell_norm = std::sqrt(cell_norm_sq);

      float max_sim = -1.0f;
      if (cell_norm > 1e-8f) {
        for (int p = 0; p < num_protos; ++p) {
          float dot = 0.0f;
          const float* p_vec = exemplar_prototypes[p].data();
          for (int c = 0; c < 256; ++c) {
            dot += p_vec[c] * embed_data[(size_t)c * 4096 + cell_offset];
          }
          float sim = dot / cell_norm;
          if (sim > max_sim) {
            max_sim = sim;
          }
        }
      } else {
        max_sim = 0.0f;
      }

      sim_grid[cell_offset] = max_sim;
    }
  }

  // 5. Build similarity grid at target resolution (feat_res: 64, 128, 256, 512, 1024)
  int max_dim = (feat_res > 64) ? feat_res : 64;
  int sim_w, sim_h;
  if (orig_w >= orig_h) {
    sim_w = max_dim;
    sim_h = std::max(16, (int)std::round((double)max_dim * orig_h / orig_w));
  } else {
    sim_h = max_dim;
    sim_w = std::max(16, (int)std::round((double)max_dim * orig_w / orig_h));
  }
  Rcpp::NumericMatrix sim_mat(sim_w, sim_h);

  int grid_res = (feat_res > 64) ? feat_res : 64;
  std::vector<float> peak_grid;

  if (grid_res <= 64) {
    // Mode 1: Fast 64x64 ViT path (blazing fast, ideal for large objects like people, leaves)
    peak_grid = sim_grid;
    for (int oy = 0; oy < sim_h; ++oy) {
      float src_y = (oy + 0.5f) / (float)sim_h * 64.0f - 0.5f;
      int y0 = std::max(0, std::min(63, (int)std::floor(src_y)));
      int y1 = std::max(0, std::min(63, y0 + 1));
      float dy = src_y - (float)y0;
      for (int ox = 0; ox < sim_w; ++ox) {
        float src_x = (ox + 0.5f) / (float)sim_w * 64.0f - 0.5f;
        int x0 = std::max(0, std::min(63, (int)std::floor(src_x)));
        int x1 = std::max(0, std::min(63, x0 + 1));
        float dx = src_x - (float)x0;

        float v00 = sim_grid[(size_t)y0 * 64 + x0];
        float v01 = sim_grid[(size_t)y0 * 64 + x1];
        float v10 = sim_grid[(size_t)y1 * 64 + x0];
        float v11 = sim_grid[(size_t)y1 * 64 + x1];

        sim_mat(ox, oy) = (1.0f - dy) * ((1.0f - dx) * v00 + dx * v01) +
                          dy * ((1.0f - dx) * v10 + dx * v11);
      }
    }
  } else if (grid_res <= 256) {
    // Mode 2: Multi-scale 256x256 feature path (SAM 2.1 high_res_feats_0 stride 4)
    std::vector<float> sim_grid_256(256 * 256, 0.0f);
    for (int r = 0; r < 256; ++r) {
      // Bilinear interpolation of semantic similarity from 64x64 embed
      float src_r = (r + 0.5f) * 0.25f - 0.5f;
      int r0 = std::max(0, std::min(63, (int)std::floor(src_r)));
      int r1 = std::max(0, std::min(63, r0 + 1));
      float dr = src_r - (float)r0;

      for (int k = 0; k < 256; ++k) {
        float src_k = (k + 0.5f) * 0.25f - 0.5f;
        int k0 = std::max(0, std::min(63, (int)std::floor(src_k)));
        int k1 = std::max(0, std::min(63, k0 + 1));
        float dk = src_k - (float)k0;

        float v00 = sim_grid[(size_t)r0 * 64 + k0];
        float v01 = sim_grid[(size_t)r0 * 64 + k1];
        float v10 = sim_grid[(size_t)r1 * 64 + k0];
        float v11 = sim_grid[(size_t)r1 * 64 + k1];
        float sem_val = (1.0f - dr) * ((1.0f - dk) * v00 + dk * v01) +
                        dr * ((1.0f - dk) * v10 + dk * v11);

        // Native 256x256 high-resolution spatial feature cosine similarity
        size_t cell_offset = (size_t)r * 256 + k;
        float hr0_norm_sq = 0.0f;
        for (int c = 0; c < 32; ++c) {
          float val = hr0_data[(size_t)c * 65536 + cell_offset];
          hr0_norm_sq += val * val;
        }
        float hr0_norm = std::sqrt(hr0_norm_sq);

        float max_hr0_sim = -1.0f;
        if (hr0_norm > 1e-8f && !exemplar_prototypes_hr0.empty()) {
          for (int p = 0; p < num_protos; ++p) {
            float dot = 0.0f;
            const float* p_vec = exemplar_prototypes_hr0[p].data();
            for (int c = 0; c < 32; ++c) {
              dot += p_vec[c] * hr0_data[(size_t)c * 65536 + cell_offset];
            }
            float sim = dot / hr0_norm;
            if (sim > max_hr0_sim) max_hr0_sim = sim;
          }
        } else {
          max_hr0_sim = sem_val;
        }

        // Fused similarity: combination of semantic context and pinpoint spatial resolution
        float fused = 0.5f * sem_val + 0.5f * max_hr0_sim;
        sim_grid_256[cell_offset] = fused;
      }
    }

    // Fill aspect-ratio aligned sim_mat directly:
    for (int oy = 0; oy < sim_h; ++oy) {
      float src_y = (oy + 0.5f) / (float)sim_h * 256.0f - 0.5f;
      int y0 = std::max(0, std::min(255, (int)std::floor(src_y)));
      int y1 = std::max(0, std::min(255, y0 + 1));
      float dy = src_y - (float)y0;
      for (int ox = 0; ox < sim_w; ++ox) {
        float src_x = (ox + 0.5f) / (float)sim_w * 256.0f - 0.5f;
        int x0 = std::max(0, std::min(255, (int)std::floor(src_x)));
        int x1 = std::max(0, std::min(255, x0 + 1));
        float dx = src_x - (float)x0;

        float v00 = sim_grid_256[(size_t)y0 * 256 + x0];
        float v01 = sim_grid_256[(size_t)y0 * 256 + x1];
        float v10 = sim_grid_256[(size_t)y1 * 256 + x0];
        float v11 = sim_grid_256[(size_t)y1 * 256 + x1];

        sim_mat(ox, oy) = (1.0f - dy) * ((1.0f - dx) * v00 + dx * v01) +
                          dy * ((1.0f - dx) * v10 + dx * v11);
      }
    }

    peak_grid = std::move(sim_grid_256);
  } else {
    // Mode 3: Dense Stride-1 Convolutional Exemplar Path (1024x1024)
    // Resolves small touching grains, seeds, and micro-objects with pixel-level precision
    grid_res = 1024;
    size_t plane_size = 1024 * 1024;

    std::vector<float> vit_1024(plane_size, 0.0f);
    for (int y = 0; y < 1024; ++y) {
      float src_y = (y + 0.5f) * (64.0f / 1024.0f) - 0.5f;
      int y0 = std::max(0, std::min(63, (int)std::floor(src_y)));
      int y1 = std::max(0, std::min(63, y0 + 1));
      float dy = src_y - (float)y0;

      for (int x = 0; x < 1024; ++x) {
        float src_x = (x + 0.5f) * (64.0f / 1024.0f) - 0.5f;
        int x0 = std::max(0, std::min(63, (int)std::floor(src_x)));
        int x1 = std::max(0, std::min(63, x0 + 1));
        float dx = src_x - (float)x0;

        float v00 = sim_grid[(size_t)y0 * 64 + x0];
        float v01 = sim_grid[(size_t)y0 * 64 + x1];
        float v10 = sim_grid[(size_t)y1 * 64 + x0];
        float v11 = sim_grid[(size_t)y1 * 64 + x1];

        vit_1024[(size_t)y * 1024 + x] = (1.0f - dy) * ((1.0f - dx) * v00 + dx * v01) +
                                         dy * ((1.0f - dx) * v10 + dx * v11);
      }
    }

    float vit_thresh = std::max(0.15f, (float)sim_threshold * 0.50f);
    std::vector<float> dense_out(plane_size, 0.0f);

    const float* img_c0 = img_values.data();
    const float* img_c1 = img_values.data() + plane_size;
    const float* img_c2 = img_values.data() + 2 * plane_size;
    int num_kernels = (int)exemplar_dense_kernels.size();

    #pragma omp parallel for schedule(dynamic, 16)
    for (int y = 0; y < 1024; ++y) {
      for (int x = 0; x < 1024; ++x) {
        float v_sim = vit_1024[(size_t)y * 1024 + x];
        if (v_sim < vit_thresh) continue;

        float max_zncc = 0.0f;
        for (int p = 0; p < num_kernels; ++p) {
          const DenseExemplarKernel& dk = exemplar_dense_kernels[p];
          int rad = dk.radius;
          size_t num_k_entries = dk.entries.size();

          float dot = 0.0f;
          float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
          float sq0 = 0.0f, sq1 = 0.0f, sq2 = 0.0f;

          if (x >= rad && x < 1024 - rad && y >= rad && y < 1024 - rad) {
            size_t center_idx = (size_t)y * 1024 + x;
            for (size_t k = 0; k < num_k_entries; ++k) {
              const ZnccKernelEntry& ke = dk.entries[k];
              size_t idx = center_idx + (size_t)ke.dy * 1024 + ke.dx;
              float i0 = img_c0[idx];
              float i1 = img_c1[idx];
              float i2 = img_c2[idx];
              float w = ke.w;

              dot += w * (i0 * ke.t_zero[0] + i1 * ke.t_zero[1] + i2 * ke.t_zero[2]);
              s0 += w * i0;
              s1 += w * i1;
              s2 += w * i2;
              sq0 += w * i0 * i0;
              sq1 += w * i1 * i1;
              sq2 += w * i2 * i2;
            }
          } else {
            for (size_t k = 0; k < num_k_entries; ++k) {
              const ZnccKernelEntry& ke = dk.entries[k];
              int px = std::max(0, std::min(1023, x + ke.dx));
              int py = std::max(0, std::min(1023, y + ke.dy));
              size_t idx = (size_t)py * 1024 + px;
              float i0 = img_c0[idx];
              float i1 = img_c1[idx];
              float i2 = img_c2[idx];
              float w = ke.w;

              dot += w * (i0 * ke.t_zero[0] + i1 * ke.t_zero[1] + i2 * ke.t_zero[2]);
              s0 += w * i0;
              s1 += w * i1;
              s2 += w * i2;
              sq0 += w * i0 * i0;
              sq1 += w * i1 * i1;
              sq2 += w * i2 * i2;
            }
          }

          float var = (sq0 - s0 * s0 * dk.inv_sum_w) +
                      (sq1 - s1 * s1 * dk.inv_sum_w) +
                      (sq2 - s2 * s2 * dk.inv_sum_w);
          float patch_std = std::sqrt(std::max(1e-6f, var));
          float zncc = dot / (patch_std * dk.tmpl_std);
          if (zncc > max_zncc) max_zncc = zncc;
        }

        if (max_zncc > 0.0f) {
          dense_out[(size_t)y * 1024 + x] = v_sim * max_zncc;
        }
      }
    }

    for (int oy = 0; oy < sim_h; ++oy) {
      int sy = std::min(1023, std::max(0, (int)std::round((oy + 0.5f) / (float)sim_h * 1024.0f - 0.5f)));
      for (int ox = 0; ox < sim_w; ++ox) {
        int sx = std::min(1023, std::max(0, (int)std::round((ox + 0.5f) / (float)sim_w * 1024.0f - 0.5f)));
        sim_mat(ox, oy) = dense_out[(size_t)sy * 1024 + sx];
      }
    }

    peak_grid = std::move(dense_out);
  }

  // 6. Detect Local Maxima (Peaks) on peak_grid
  struct CandPeak {
    int r;
    int k;
    float sim;
  };
  std::vector<CandPeak> cand_peaks;

  float peak_thresh = (float)sim_threshold;
  if (grid_res >= 512) {
    float max_peak_val = 0.0f;
    for (float v : peak_grid) {
      if (v > max_peak_val) max_peak_val = v;
    }
    if (max_peak_val > 0.0f) {
      peak_thresh = std::min((float)sim_threshold, max_peak_val * (float)sim_threshold);
    }
  }

  for (int r = 1; r < grid_res - 1; ++r) {
    for (int k = 1; k < grid_res - 1; ++k) {
      float cur_sim = peak_grid[(size_t)r * grid_res + k];
      if (cur_sim < peak_thresh) continue;

      bool is_peak = true;
      for (int dr = -1; dr <= 1; ++dr) {
        for (int dc = -1; dc <= 1; ++dc) {
          if (dr == 0 && dc == 0) continue;
          if (peak_grid[(size_t)(r + dr) * grid_res + (k + dc)] > cur_sim) {
            is_peak = false;
            break;
          }
        }
        if (!is_peak) break;
      }
      if (is_peak) {
        cand_peaks.push_back({ r, k, cur_sim });
      }
    }
  }

  std::sort(cand_peaks.begin(), cand_peaks.end(), [](const CandPeak& a, const CandPeak& b) {
    return a.sim > b.sim;
  });

  // Spatial NMS on candidate peaks (min_dist proportionally scaled to 1024-space)
  float cell_size = 1024.0f / (float)grid_res;
  float min_dist_1024 = (float)(min_dist * 1024.0 / std::max(orig_w, orig_h));
  float grid_min_dist = min_dist_1024 / cell_size;
  if (grid_min_dist < 1.0f) grid_min_dist = 1.0f;

  std::vector<CandPeak> accepted_peaks;
  for (size_t i = 0; i < cand_peaks.size(); ++i) {
    bool suppressed = false;
    for (size_t j = 0; j < accepted_peaks.size(); ++j) {
      float dr = (float)(cand_peaks[i].r - accepted_peaks[j].r);
      float dc = (float)(cand_peaks[i].k - accepted_peaks[j].k);
      if (std::sqrt(dr * dr + dc * dc) < grid_min_dist) {
        suppressed = true;
        break;
      }
    }
    if (!suppressed) {
      accepted_peaks.push_back(cand_peaks[i]);
      if ((int)accepted_peaks.size() >= max_objects) break;
    }
  }

  // 7. Run SAM 2.1 Decoder for Each Accepted Peak
  struct PersamInst {
    double x1, y1, x2, y2;
    float sim_score;
    float iou_pred;
    Rcpp::NumericMatrix mask;
    std::vector<uint8_t> bin_mask;
  };
  std::vector<PersamInst> instances;

  for (size_t p = 0; p < accepted_peaks.size(); ++p) {
    float cand_px = (accepted_peaks[p].k + 0.5f) * cell_size;
    float cand_py = (accepted_peaks[p].r + 0.5f) * cell_size;

    float pt_coords[2] = { cand_px, cand_py };
    float pt_label[1] = { 1.0f };

    int64_t pt_shape[3] = { 1, 1, 2 };
    int64_t lbl_shape[2] = { 1, 1 };

    OrtValue* cand_pts_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, pt_coords, 2 * sizeof(float),
      pt_shape, 3, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &cand_pts_tensor
    );

    OrtValue* cand_lbl_tensor = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, pt_label, 1 * sizeof(float),
      lbl_shape, 2, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &cand_lbl_tensor
    );

    const OrtValue* cand_dec_in_values[] = {
      enc_outputs[2], enc_outputs[0], enc_outputs[1],
      cand_pts_tensor, cand_lbl_tensor, mask_input_tensor, has_mask_tensor
    };

    OrtValue* cand_dec_outputs[2] = { NULL, NULL };
    OrtStatus* dec_status = ort->Run(
      session_dec, run_options, dec_in_names, cand_dec_in_values, 7, dec_out_names, 2, cand_dec_outputs
    );

    if (dec_status == NULL) {
      OrtTensorTypeAndShapeInfo* iou_shape_info = NULL;
      (void)ort->GetTensorTypeAndShape(cand_dec_outputs[1], &iou_shape_info);
      size_t iou_dims_count = 0;
      (void)ort->GetDimensionsCount(iou_shape_info, &iou_dims_count);
      std::vector<int64_t> iou_dims(iou_dims_count);
      (void)ort->GetDimensions(iou_shape_info, iou_dims.data(), iou_dims_count);
      int num_masks = (iou_dims_count >= 2) ? (int)iou_dims[1] : 1;

      float* iou_data = NULL;
      (void)ort->GetTensorMutableData(cand_dec_outputs[1], (void**)&iou_data);

      float* mask_data = NULL;
      (void)ort->GetTensorMutableData(cand_dec_outputs[0], (void**)&mask_data);

      int best_idx = 0;
      float best_score = -1e9f;
      for (int m = 0; m < num_masks; ++m) {
        size_t offset = (size_t)m * mask_h * mask_w;
        int fg_c = 0;
        for (int j = 0; j < mask_h * mask_w; ++j) {
          if (mask_data[offset + j] > 0.0f) fg_c++;
        }
        float fg_ratio = (float)fg_c / (float)(mask_h * mask_w);
        float score = iou_data[m];
        if (fg_ratio > 0.80f) score -= 5.0f;

        // Area prior from exemplar: penalize masks that merge adjacent objects into clusters
        if (avg_ex_area > 5.0f && fg_c > 0) {
          float area_ratio = (float)fg_c / avg_ex_area;
          if (area_ratio > 1.35f) {
            // Mask is noticeably larger than exemplar: penalize multi-object cluster
            score -= (area_ratio - 1.35f) * 2.0f;
          } else if (area_ratio < 0.25f) {
            // Tiny fragment
            score -= (0.25f - area_ratio) * 2.0f;
          }
        }

        if (score > best_score) {
          best_score = score;
          best_idx = m;
        }
      }

      size_t slice_offset = (size_t)best_idx * mask_h * mask_w;
      Rcpp::NumericMatrix inst_mat(mask_w, mask_h);
      std::vector<uint8_t> bin_vec(mask_w * mask_h, 0);
      int min_r = mask_h, max_r = -1, min_c = mask_w, max_c = -1;
      int fg_total = 0;

      for (int r = 0; r < mask_h; ++r) {
        for (int c = 0; c < mask_w; ++c) {
          float logit = mask_data[slice_offset + r * mask_w + c];
          double prob = 1.0 / (1.0 + std::exp(-(double)logit));
          inst_mat(c, r) = prob;
          if (prob > 0.5) {
            bin_vec[r * mask_w + c] = 1;
            fg_total++;
            if (r < min_r) min_r = r;
            if (r > max_r) max_r = r;
            if (c < min_c) min_c = c;
            if (c > max_c) max_c = c;
          }
        }
      }

      if (fg_total >= 10 && fg_total <= (int)(mask_h * mask_w * 0.85f)) {
        double bx1 = (double)min_c / 256.0 * orig_w;
        double by1 = (double)min_r / 256.0 * orig_h;
        double bx2 = (double)(max_c + 1) / 256.0 * orig_w;
        double by2 = (double)(max_r + 1) / 256.0 * orig_h;

        instances.push_back({
          bx1, by1, bx2, by2,
          accepted_peaks[p].sim,
          best_score,
          inst_mat,
          bin_vec
        });
      }

      ort->ReleaseTensorTypeAndShapeInfo(iou_shape_info);
      ort->ReleaseValue(cand_dec_outputs[0]);
      ort->ReleaseValue(cand_dec_outputs[1]);
    } else {
      ort->ReleaseStatus(dec_status);
    }

    ort->ReleaseValue(cand_pts_tensor);
    ort->ReleaseValue(cand_lbl_tensor);
  }

  // 7. Instance IoU NMS
  std::vector<bool> suppressed(instances.size(), false);
  for (size_t i = 0; i < instances.size(); ++i) {
    if (suppressed[i]) continue;
    for (size_t j = i + 1; j < instances.size(); ++j) {
      if (suppressed[j]) continue;
      int intersection = 0;
      int union_cnt = 0;
      for (size_t k = 0; k < instances[i].bin_mask.size(); ++k) {
        uint8_t a = instances[i].bin_mask[k];
        uint8_t b = instances[j].bin_mask[k];
        if (a && b) intersection++;
        if (a || b) union_cnt++;
      }
      double iou = (union_cnt > 0) ? ((double)intersection / (double)union_cnt) : 0.0;
      if (iou > iou_threshold) {
        if (instances[i].sim_score >= instances[j].sim_score) {
          suppressed[j] = true;
        } else {
          suppressed[i] = true;
          break;
        }
      }
    }
  }

  std::vector<PersamInst> kept_instances;
  for (size_t i = 0; i < instances.size(); ++i) {
    if (!suppressed[i]) kept_instances.push_back(instances[i]);
  }

  // 8. Prepare Return Structures
  int num_final = (int)kept_instances.size();
  Rcpp::NumericMatrix res_boxes(num_final, 4);
  Rcpp::NumericVector res_scores(num_final);
  Rcpp::List res_masks(num_final);

  for (int i = 0; i < num_final; ++i) {
    res_boxes(i, 0) = kept_instances[i].x1;
    res_boxes(i, 1) = kept_instances[i].y1;
    res_boxes(i, 2) = kept_instances[i].x2;
    res_boxes(i, 3) = kept_instances[i].y2;
    res_scores[i]   = kept_instances[i].sim_score;
    res_masks[i]    = kept_instances[i].mask;
  }

  // Cleanup
  ort->ReleaseRunOptions(run_options);
  ort->ReleaseValue(input_img_tensor);
  ort->ReleaseValue(mask_input_tensor);
  ort->ReleaseValue(has_mask_tensor);
  for (int i = 0; i < 3; ++i) ort->ReleaseValue(enc_outputs[i]);

  return Rcpp::List::create(
    Rcpp::Named("boxes") = res_boxes,
    Rcpp::Named("scores") = res_scores,
    Rcpp::Named("masks") = res_masks,
    Rcpp::Named("similarity_map") = sim_mat
  );
}

// [[Rcpp::export]]
Rcpp::List inspect_onnx_model_cpp(std::string model_path, std::string lib_path) {
  const OrtApi* ort = get_ort_api(lib_path);
  OrtEnv* env = NULL;
  (void)ort->CreateEnv(ORT_LOGGING_LEVEL_ERROR, "pliman_inspect", &env);
  OrtSessionOptions* session_options = NULL;
  (void)ort->CreateSessionOptions(&session_options);
  (void)ort->SetSessionLogSeverityLevel(session_options, ORT_LOGGING_LEVEL_ERROR);
  OrtSession* session = NULL;
#ifdef _WIN32
  int size_needed = MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.size(), NULL, 0);
  std::wstring wmodel_path(size_needed, 0);
  MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.size(), &wmodel_path[0], size_needed);
  (void)ort->CreateSession(env, wmodel_path.c_str(), session_options, &session);
#else
  (void)ort->CreateSession(env, model_path.c_str(), session_options, &session);
#endif

  if (!session) {
    ort->ReleaseSessionOptions(session_options);
    ort->ReleaseEnv(env);
    return Rcpp::List::create(Rcpp::Named("error") = "Failed to create session");
  }

  OrtAllocator* allocator = NULL;
  (void)ort->GetAllocatorWithDefaultOptions(&allocator);

  size_t num_inputs = 0;
  (void)ort->SessionGetInputCount(session, &num_inputs);
  Rcpp::CharacterVector in_names(num_inputs);
  Rcpp::IntegerVector in_types(num_inputs);
  Rcpp::List in_shapes(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    char* name = NULL;
    (void)ort->SessionGetInputName(session, i, allocator, &name);
    in_names[i] = std::string(name);
    allocator->Free(allocator, name);

    OrtTypeInfo* type_info = NULL;
    (void)ort->SessionGetInputTypeInfo(session, i, &type_info);
    const OrtTensorTypeAndShapeInfo* shape_info = NULL;
    (void)ort->CastTypeInfoToTensorInfo(type_info, &shape_info);
    ONNXTensorElementDataType elem_type = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;
    if (shape_info) {
      (void)ort->GetTensorElementType(shape_info, &elem_type);
      size_t num_dims = 0;
      (void)ort->GetDimensionsCount(shape_info, &num_dims);
      std::vector<int64_t> dims(num_dims);
      (void)ort->GetDimensions(shape_info, dims.data(), num_dims);
      Rcpp::NumericVector d_vec(num_dims);
      for (size_t d = 0; d < num_dims; ++d) d_vec[d] = (double)dims[d];
      in_shapes[i] = d_vec;
    } else {
      in_shapes[i] = Rcpp::NumericVector::create();
    }
    in_types[i] = (int)elem_type;
    ort->ReleaseTypeInfo(type_info);
  }

  size_t num_outputs = 0;
  (void)ort->SessionGetOutputCount(session, &num_outputs);
  Rcpp::CharacterVector out_names(num_outputs);
  Rcpp::List out_shapes(num_outputs);
  Rcpp::IntegerVector out_types(num_outputs);
  for (size_t i = 0; i < num_outputs; ++i) {
    char* name = NULL;
    (void)ort->SessionGetOutputName(session, i, allocator, &name);
    out_names[i] = std::string(name);
    allocator->Free(allocator, name);

    OrtTypeInfo* type_info = NULL;
    (void)ort->SessionGetOutputTypeInfo(session, i, &type_info);
    const OrtTensorTypeAndShapeInfo* shape_info = NULL;
    (void)ort->CastTypeInfoToTensorInfo(type_info, &shape_info);
    ONNXTensorElementDataType elem_type = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;
    if (shape_info) {
      (void)ort->GetTensorElementType(shape_info, &elem_type);
      size_t num_dims = 0;
      (void)ort->GetDimensionsCount(shape_info, &num_dims);
      std::vector<int64_t> dims(num_dims);
      (void)ort->GetDimensions(shape_info, dims.data(), num_dims);
      Rcpp::NumericVector d_vec(num_dims);
      for (size_t d = 0; d < num_dims; ++d) d_vec[d] = (double)dims[d];
      out_shapes[i] = d_vec;
    } else {
      out_shapes[i] = Rcpp::NumericVector::create();
    }
    out_types[i] = (int)elem_type;
    ort->ReleaseTypeInfo(type_info);
  }

  ort->ReleaseSession(session);
  ort->ReleaseSessionOptions(session_options);
  ort->ReleaseEnv(env);

  return Rcpp::List::create(
    Rcpp::Named("inputs") = in_names,
    Rcpp::Named("types") = in_types,
    Rcpp::Named("shapes") = in_shapes,
    Rcpp::Named("outputs") = out_names,
    Rcpp::Named("output_shapes") = out_shapes,
    Rcpp::Named("output_types") = out_types
  );
}

// --------------------------------------------------------------------------
// Helper: Dynamically retrieve all input and output names from session
// --------------------------------------------------------------------------
static void get_session_io_names(const OrtApi* ort, OrtSession* session,
                                 std::vector<std::string>& in_names_str,
                                 std::vector<const char*>& in_names,
                                 std::vector<std::string>& out_names_str,
                                 std::vector<const char*>& out_names) {
  OrtAllocator* allocator = NULL;
  (void)ort->GetAllocatorWithDefaultOptions(&allocator);
  size_t num_inputs = 0;
  (void)ort->SessionGetInputCount(session, &num_inputs);
  in_names_str.resize(num_inputs);
  in_names.resize(num_inputs);
  for (size_t i = 0; i < num_inputs; ++i) {
    char* name = NULL;
    (void)ort->SessionGetInputName(session, i, allocator, &name);
    in_names_str[i] = std::string(name ? name : "input");
    if (name) allocator->Free(allocator, name);
    in_names[i] = in_names_str[i].c_str();
  }

  size_t num_outputs = 0;
  (void)ort->SessionGetOutputCount(session, &num_outputs);
  out_names_str.resize(num_outputs);
  out_names.resize(num_outputs);
  for (size_t i = 0; i < num_outputs; ++i) {
    char* name = NULL;
    (void)ort->SessionGetOutputName(session, i, allocator, &name);
    out_names_str[i] = std::string(name ? name : "output");
    if (name) allocator->Free(allocator, name);
    out_names[i] = out_names_str[i].c_str();
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix run_depth_anything_cpp(
    Rcpp::NumericVector tensor_vec,
    int in_w,
    int in_h,
    int orig_w,
    int orig_h,
    std::string model_path,
    std::string lib_path,
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  int64_t in_shape[4] = {1, 3, in_h, in_w};
  size_t total_floats = (size_t)(3 * in_h * in_w);
  std::vector<float> input_vals(total_floats);
  for (size_t i = 0; i < total_floats; ++i) {
    input_vals[i] = (float)tensor_vec[i];
  }

  OrtValue* in_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, input_vals.data(), total_floats * sizeof(float),
    in_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_tensor
  );

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtValue* out_tensor = NULL;
  OrtStatus* status = ort->Run(
    session, run_options, in_names.data(), (const OrtValue* const*)&in_tensor,
    1, out_names.data(), 1, &out_tensor
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(in_tensor);
    Rcpp::stop("Depth Anything V2 run failed: " + msg);
  }

  float* out_data = NULL;
  (void)ort->GetTensorMutableData(out_tensor, (void**)&out_data);

  Rcpp::NumericMatrix depth_mat(orig_w, orig_h);
  for (int oy = 0; oy < orig_h; ++oy) {
    float src_y = (oy + 0.5f) / (float)orig_h * (float)in_h - 0.5f;
    int y0 = std::max(0, std::min(in_h - 1, (int)std::floor(src_y)));
    int y1 = std::max(0, std::min(in_h - 1, y0 + 1));
    float dy = src_y - (float)y0;

    for (int ox = 0; ox < orig_w; ++ox) {
      float src_x = (ox + 0.5f) / (float)orig_w * (float)in_w - 0.5f;
      int x0 = std::max(0, std::min(in_w - 1, (int)std::floor(src_x)));
      int x1 = std::max(0, std::min(in_w - 1, x0 + 1));
      float dx = src_x - (float)x0;

      float v00 = out_data[(size_t)y0 * in_w + x0];
      float v01 = out_data[(size_t)y0 * in_w + x1];
      float v10 = out_data[(size_t)y1 * in_w + x0];
      float v11 = out_data[(size_t)y1 * in_w + x1];

      depth_mat(ox, oy) = (1.0f - dy) * ((1.0f - dx) * v00 + dx * v01) +
                          dy * ((1.0f - dx) * v10 + dx * v11);
    }
  }

  ort->ReleaseValue(out_tensor);
  ort->ReleaseValue(in_tensor);
  ort->ReleaseRunOptions(run_options);

  return depth_mat;
}

// [[Rcpp::export]]
Rcpp::List run_dinov2_cpp(
    Rcpp::NumericVector tensor_vec,
    int in_w,
    int in_h,
    int patch_size = 14,
    bool return_pca = true,
    std::string model_path = "",
    std::string lib_path = "",
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  int64_t in_shape[4] = {1, 3, in_h, in_w};
  size_t total_floats = (size_t)(3 * in_h * in_w);
  std::vector<float> input_vals(total_floats);
  for (size_t i = 0; i < total_floats; ++i) input_vals[i] = (float)tensor_vec[i];

  OrtValue* in_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, input_vals.data(), total_floats * sizeof(float),
    in_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_tensor
  );

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  OrtValue* out_tensor = NULL;
  OrtStatus* status = ort->Run(
    session, run_options, in_names.data(), (const OrtValue* const*)&in_tensor,
    1, out_names.data(), 1, &out_tensor
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(in_tensor);
    Rcpp::stop("DINOv2 run failed: " + msg);
  }

  OrtTensorTypeAndShapeInfo* shape_info = NULL;
  (void)ort->GetTensorTypeAndShape(out_tensor, &shape_info);
  size_t num_dims = 0;
  (void)ort->GetDimensionsCount(shape_info, &num_dims);
  std::vector<int64_t> dims(num_dims);
  (void)ort->GetDimensions(shape_info, dims.data(), num_dims);
  ort->ReleaseTensorTypeAndShapeInfo(shape_info);

  int total_tokens = (num_dims >= 2) ? (int)dims[1] : 0;
  int embed_dim = (num_dims >= 3) ? (int)dims[2] : 0;
  int wp = in_w / patch_size;
  int hp = in_h / patch_size;
  int num_patches = wp * hp;

  float* out_data = NULL;
  (void)ort->GetTensorMutableData(out_tensor, (void**)&out_data);

  int cls_idx = 0;
  int patch_start_idx = (total_tokens > num_patches) ? 1 : 0;

  Rcpp::NumericVector cls_vec(embed_dim);
  if (total_tokens > 0 && embed_dim > 0) {
    for (int c = 0; c < embed_dim; ++c) {
      cls_vec[c] = out_data[(size_t)cls_idx * embed_dim + c];
    }
  }

  Rcpp::NumericMatrix r_mat(wp, hp);
  Rcpp::NumericMatrix g_mat(wp, hp);
  Rcpp::NumericMatrix b_mat(wp, hp);

  if (return_pca && num_patches > 3 && embed_dim >= 3 && out_data != NULL) {
    arma::mat X(num_patches, embed_dim);
    for (int p = 0; p < num_patches; ++p) {
      size_t tok_offset = (size_t)(patch_start_idx + p) * embed_dim;
      for (int c = 0; c < embed_dim; ++c) {
        X(p, c) = out_data[tok_offset + c];
      }
    }

    arma::rowvec mean_x = arma::mean(X, 0);
    X.each_row() -= mean_x;

    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd_econ(U, s, V, X);

    arma::mat P = X * V.cols(0, 2);

    for (int comp = 0; comp < 3; ++comp) {
      double min_val = P.col(comp).min();
      double max_val = P.col(comp).max();
      double range = (max_val > min_val) ? (max_val - min_val) : 1.0;
      for (int p = 0; p < num_patches; ++p) {
        double norm_val = (P(p, comp) - min_val) / range;
        int py = p / wp;
        int px = p % wp;
        if (comp == 0) r_mat(px, py) = norm_val;
        else if (comp == 1) g_mat(px, py) = norm_val;
        else if (comp == 2) b_mat(px, py) = norm_val;
      }
    }
  }

  ort->ReleaseValue(out_tensor);
  ort->ReleaseValue(in_tensor);
  ort->ReleaseRunOptions(run_options);

  return Rcpp::List::create(
    Rcpp::Named("cls_token") = cls_vec,
    Rcpp::Named("pca_r") = r_mat,
    Rcpp::Named("pca_g") = g_mat,
    Rcpp::Named("pca_b") = b_mat,
    Rcpp::Named("num_patches") = num_patches,
    Rcpp::Named("embed_dim") = embed_dim,
    Rcpp::Named("wp") = wp,
    Rcpp::Named("hp") = hp
  );
}

// [[Rcpp::export]]
Rcpp::List run_yolo_cpp(
    Rcpp::NumericVector tensor_vec,
    double orig_w,
    double orig_h,
    double conf_threshold = 0.25,
    double iou_threshold = 0.45,
    std::string model_path = "",
    std::string lib_path = "",
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  int64_t in_shape[4] = {1, 3, 640, 640};
  size_t total_floats = 3 * 640 * 640;
  std::vector<float> input_vals(total_floats);
  for (size_t i = 0; i < total_floats; ++i) input_vals[i] = (float)tensor_vec[i];

  OrtValue* in_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, input_vals.data(), total_floats * sizeof(float),
    in_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_tensor
  );

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  size_t num_outputs = out_names.size();
  std::vector<OrtValue*> out_tensors(num_outputs, NULL);
  OrtStatus* status = ort->Run(
    session, run_options, in_names.data(), (const OrtValue* const*)&in_tensor,
    1, out_names.data(), num_outputs, out_tensors.data()
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(in_tensor);
    Rcpp::stop("YOLO run failed: " + msg);
  }

  OrtTensorTypeAndShapeInfo* shape_info0 = NULL;
  (void)ort->GetTensorTypeAndShape(out_tensors[0], &shape_info0);
  size_t dims0_cnt = 0;
  (void)ort->GetDimensionsCount(shape_info0, &dims0_cnt);
  std::vector<int64_t> dims0(dims0_cnt);
  (void)ort->GetDimensions(shape_info0, dims0.data(), dims0_cnt);
  ort->ReleaseTensorTypeAndShapeInfo(shape_info0);

  bool is_end2end = false;
  int e2e_num_dets = 0;
  int e2e_num_cols = 0;
  if (dims0_cnt == 3 && (dims0[2] == 6 || dims0[2] == 38 || dims0[2] == 57)) {
    is_end2end = true;
    e2e_num_dets = (int)dims0[1];
    e2e_num_cols = (int)dims0[2];
  }

  bool is_transposed = false;
  int num_features = 0;
  int num_anchors = 0;
  if (!is_end2end && dims0_cnt == 3) {
    if (dims0[1] < dims0[2]) {
      num_features = (int)dims0[1];
      num_anchors  = (int)dims0[2];
      is_transposed = false;
    } else {
      num_anchors  = (int)dims0[1];
      num_features = (int)dims0[2];
      is_transposed = true;
    }
  }

  bool is_seg = (num_outputs >= 2);
  int num_mask_coeffs = is_seg ? 32 : 0;
  int num_classes = is_end2end ? 80 : (num_features - 4 - num_mask_coeffs);
  if (num_classes < 1) num_classes = 1;

  float* data0 = NULL;
  (void)ort->GetTensorMutableData(out_tensors[0], (void**)&data0);

  struct YoloCand {
    float x1, y1, x2, y2;
    float score;
    int class_id;
    std::vector<float> mask_coeffs;
    std::vector<float> keypoints;
  };
  std::vector<YoloCand> candidates;

  if (is_end2end) {
    for (int i = 0; i < e2e_num_dets; ++i) {
      size_t offset = (size_t)i * e2e_num_cols;
      float x1 = data0[offset + 0];
      float y1 = data0[offset + 1];
      float x2 = data0[offset + 2];
      float y2 = data0[offset + 3];
      float s = data0[offset + 4];
      int cls = (int)std::round(data0[offset + 5]);

      if (s >= (float)conf_threshold) {
        YoloCand cand;
        cand.x1 = x1;
        cand.y1 = y1;
        cand.x2 = x2;
        cand.y2 = y2;
        cand.score = s;
        cand.class_id = cls;

        if (e2e_num_cols >= 38 && is_seg) {
          cand.mask_coeffs.resize(32);
          for (int m = 0; m < 32; ++m) {
            cand.mask_coeffs[m] = data0[offset + 6 + m];
          }
        }
        if (e2e_num_cols >= 57) {
          cand.keypoints.resize(51);
          for (int k = 0; k < 51; ++k) {
            cand.keypoints[k] = data0[offset + 6 + k];
          }
        }
        candidates.push_back(cand);
      }
    }
  } else {
    for (int a = 0; a < num_anchors; ++a) {
      float cx, cy, w, h;
      if (!is_transposed) {
        cx = data0[(size_t)0 * num_anchors + a];
        cy = data0[(size_t)1 * num_anchors + a];
        w  = data0[(size_t)2 * num_anchors + a];
        h  = data0[(size_t)3 * num_anchors + a];
      } else {
        size_t anchor_offset = (size_t)a * num_features;
        cx = data0[anchor_offset + 0];
        cy = data0[anchor_offset + 1];
        w  = data0[anchor_offset + 2];
        h  = data0[anchor_offset + 3];
      }

      float max_s = 0.0f;
      int best_cls = 0;
      for (int c = 0; c < num_classes; ++c) {
        float s = !is_transposed ? data0[(size_t)(4 + c) * num_anchors + a] :
                                   data0[(size_t)a * num_features + 4 + c];
        if (s > max_s) {
          max_s = s;
          best_cls = c;
        }
      }

      if (max_s >= (float)conf_threshold) {
        float x1 = cx - w / 2.0f;
        float y1 = cy - h / 2.0f;
        float x2 = cx + w / 2.0f;
        float y2 = cy + h / 2.0f;

        YoloCand cand;
        cand.x1 = x1; cand.y1 = y1; cand.x2 = x2; cand.y2 = y2;
        cand.score = max_s;
        cand.class_id = best_cls;

        if (is_seg) {
          cand.mask_coeffs.resize(num_mask_coeffs);
          for (int m = 0; m < num_mask_coeffs; ++m) {
            cand.mask_coeffs[m] = !is_transposed ? data0[(size_t)(4 + num_classes + m) * num_anchors + a] :
                                                   data0[(size_t)a * num_features + 4 + num_classes + m];
          }
        }
        candidates.push_back(cand);
      }
    }
  }

  std::sort(candidates.begin(), candidates.end(), [](const YoloCand& a, const YoloCand& b) {
    return a.score > b.score;
  });

  std::vector<bool> suppressed(candidates.size(), false);
  std::vector<YoloCand> kept;

  for (size_t i = 0; i < candidates.size(); ++i) {
    if (suppressed[i]) continue;
    kept.push_back(candidates[i]);

    float area_i = std::max(0.0f, candidates[i].x2 - candidates[i].x1) *
                   std::max(0.0f, candidates[i].y2 - candidates[i].y1);

    for (size_t j = i + 1; j < candidates.size(); ++j) {
      if (suppressed[j]) continue;
      if (candidates[i].class_id != candidates[j].class_id) continue;

      float inter_x1 = std::max(candidates[i].x1, candidates[j].x1);
      float inter_y1 = std::max(candidates[i].y1, candidates[j].y1);
      float inter_x2 = std::min(candidates[i].x2, candidates[j].x2);
      float inter_y2 = std::min(candidates[i].y2, candidates[j].y2);

      float inter_w = std::max(0.0f, inter_x2 - inter_x1);
      float inter_h = std::max(0.0f, inter_y2 - inter_y1);
      float inter_area = inter_w * inter_h;

      float area_j = std::max(0.0f, candidates[j].x2 - candidates[j].x1) *
                     std::max(0.0f, candidates[j].y2 - candidates[j].y1);

      float union_area = area_i + area_j - inter_area;
      float iou = (union_area > 0.0f) ? (inter_area / union_area) : 0.0f;

      if (iou > (float)iou_threshold) {
        suppressed[j] = true;
      }
    }
  }

  double gain = std::min(640.0 / orig_w, 640.0 / orig_h);
  double pad_x = (640.0 - orig_w * gain) / 2.0;
  double pad_y = (640.0 - orig_h * gain) / 2.0;

  int num_kept = (int)kept.size();
  Rcpp::NumericMatrix boxes(num_kept, 4);
  Rcpp::NumericVector scores(num_kept);
  Rcpp::IntegerVector class_ids(num_kept);

  bool has_kpts = false;
  for (int i = 0; i < num_kept; ++i) {
    if (!kept[i].keypoints.empty()) {
      has_kpts = true;
      break;
    }
  }

  Rcpp::NumericMatrix keypoints_mat(num_kept, has_kpts ? 51 : 0);

  for (int i = 0; i < num_kept; ++i) {
    double x1 = (kept[i].x1 - pad_x) / gain;
    double y1 = (kept[i].y1 - pad_y) / gain;
    double x2 = (kept[i].x2 - pad_x) / gain;
    double y2 = (kept[i].y2 - pad_y) / gain;

    boxes(i, 0) = std::max(0.0, std::min(orig_w, x1));
    boxes(i, 1) = std::max(0.0, std::min(orig_h, y1));
    boxes(i, 2) = std::max(0.0, std::min(orig_w, x2));
    boxes(i, 3) = std::max(0.0, std::min(orig_h, y2));
    scores[i] = kept[i].score;
    class_ids[i] = kept[i].class_id;

    if (has_kpts && kept[i].keypoints.size() >= 51) {
      for (int k = 0; k < 17; ++k) {
        float kx = kept[i].keypoints[k * 3 + 0];
        float ky = kept[i].keypoints[k * 3 + 1];
        float kconf = kept[i].keypoints[k * 3 + 2];
        double x_orig = (kx - pad_x) / gain;
        double y_orig = (ky - pad_y) / gain;
        keypoints_mat(i, k * 3 + 0) = std::max(0.0, std::min(orig_w, x_orig));
        keypoints_mat(i, k * 3 + 1) = std::max(0.0, std::min(orig_h, y_orig));
        keypoints_mat(i, k * 3 + 2) = kconf;
      }
    }
  }

  int ow = (int)std::round(orig_w);
  int oh = (int)std::round(orig_h);
  Rcpp::IntegerMatrix labels(ow, oh);
  Rcpp::LogicalMatrix mask(ow, oh);

  if (is_seg && num_kept > 0) {
    float* proto_data = NULL;
    (void)ort->GetTensorMutableData(out_tensors[1], (void**)&proto_data);
    const size_t proto_plane = 160 * 160;

    for (int i = 0; i < num_kept; ++i) {
      if (kept[i].mask_coeffs.empty()) continue;
      int ox1 = std::max(0, (int)std::floor(boxes(i, 0)));
      int oy1 = std::max(0, (int)std::floor(boxes(i, 1)));
      int ox2 = std::min(ow - 1, (int)std::ceil(boxes(i, 2)));
      int oy2 = std::min(oh - 1, (int)std::ceil(boxes(i, 3)));
      if (ox1 > ox2 || oy1 > oy2) continue;

      // 1. Precompute linear combination of 32 prototype channels for instance i (160x160)
      std::vector<float> inst_proto(proto_plane, 0.0f);
      for (int m = 0; m < 32; ++m) {
        float coeff = kept[i].mask_coeffs[m];
        if (std::abs(coeff) < 1e-7f) continue;
        const float* p_plane = proto_data + (size_t)m * proto_plane;
        for (size_t p = 0; p < proto_plane; ++p) {
          inst_proto[p] += coeff * p_plane[p];
        }
      }

      // 2. High-resolution continuous bilinear interpolation on the original image grid
      for (int oy = oy1; oy <= oy2; ++oy) {
        float can_y = ((float)oy + 0.5f) * (float)gain + (float)pad_y;
        float proto_y = (can_y * 0.25f) - 0.5f;
        int y0 = (int)std::floor(proto_y);
        int y1 = y0 + 1;
        float wy = proto_y - (float)y0;
        int cy0 = std::max(0, std::min(159, y0));
        int cy1 = std::max(0, std::min(159, y1));

        for (int ox = ox1; ox <= ox2; ++ox) {
          float can_x = ((float)ox + 0.5f) * (float)gain + (float)pad_x;
          float proto_x = (can_x * 0.25f) - 0.5f;
          int x0 = (int)std::floor(proto_x);
          int x1 = x0 + 1;
          float wx = proto_x - (float)x0;
          int cx0 = std::max(0, std::min(159, x0));
          int cx1 = std::max(0, std::min(159, x1));

          float v00 = inst_proto[(size_t)cy0 * 160 + cx0];
          float v10 = inst_proto[(size_t)cy0 * 160 + cx1];
          float v01 = inst_proto[(size_t)cy1 * 160 + cx0];
          float v11 = inst_proto[(size_t)cy1 * 160 + cx1];

          float val = (1.0f - wx) * (1.0f - wy) * v00 +
                      wx * (1.0f - wy) * v10 +
                      (1.0f - wx) * wy * v01 +
                      wx * wy * v11;

          if (val > 0.0f) {
            labels(ox, oy) = i + 1;
            mask(ox, oy) = true;
          }
        }
      }
    }
  }

  for (size_t i = 0; i < num_outputs; ++i) ort->ReleaseValue(out_tensors[i]);
  ort->ReleaseValue(in_tensor);
  ort->ReleaseRunOptions(run_options);

  return Rcpp::List::create(
    Rcpp::Named("boxes") = boxes,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("class_ids") = class_ids,
    Rcpp::Named("labels") = labels,
    Rcpp::Named("mask") = mask,
    Rcpp::Named("keypoints") = keypoints_mat
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector run_yolo_cls_cpp(
    Rcpp::NumericVector tensor_vec,
    std::string model_path = "",
    std::string lib_path = "",
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  int64_t in_shape[4] = {1, 3, 640, 640};
  size_t total_floats = 3 * 640 * 640;
  std::vector<float> input_vals(total_floats);
  for (size_t i = 0; i < total_floats; ++i) input_vals[i] = (float)tensor_vec[i];

  OrtValue* in_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, input_vals.data(), total_floats * sizeof(float),
    in_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_tensor
  );

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  size_t num_outputs = out_names.size();
  std::vector<OrtValue*> out_tensors(num_outputs, NULL);
  OrtStatus* status = ort->Run(
    session, run_options, in_names.data(), (const OrtValue* const*)&in_tensor,
    1, out_names.data(), num_outputs, out_tensors.data()
  );

  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(in_tensor);
    Rcpp::stop("YOLO classification run failed: " + msg);
  }

  OrtTensorTypeAndShapeInfo* shape_info0 = NULL;
  (void)ort->GetTensorTypeAndShape(out_tensors[0], &shape_info0);
  size_t dims0_cnt = 0;
  (void)ort->GetDimensionsCount(shape_info0, &dims0_cnt);
  std::vector<int64_t> dims0(dims0_cnt);
  (void)ort->GetDimensions(shape_info0, dims0.data(), dims0_cnt);
  ort->ReleaseTensorTypeAndShapeInfo(shape_info0);

  int num_classes = 1000;
  if (dims0_cnt == 2) {
    num_classes = (int)dims0[1];
  } else if (dims0_cnt == 1) {
    num_classes = (int)dims0[0];
  }

  float* data0 = NULL;
  (void)ort->GetTensorMutableData(out_tensors[0], (void**)&data0);

  float max_val = data0[0];
  for (int c = 1; c < num_classes; ++c) {
    if (data0[c] > max_val) max_val = data0[c];
  }

  double sum_exp = 0.0;
  std::vector<double> probs(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    probs[c] = std::exp((double)(data0[c] - max_val));
    sum_exp += probs[c];
  }

  Rcpp::NumericVector res(num_classes);
  for (int c = 0; c < num_classes; ++c) {
    res[c] = (sum_exp > 0.0) ? (probs[c] / sum_exp) : 0.0;
  }

  for (size_t i = 0; i < num_outputs; ++i) ort->ReleaseValue(out_tensors[i]);
  ort->ReleaseValue(in_tensor);
  ort->ReleaseRunOptions(run_options);

  return res;
}

// [[Rcpp::export]]
Rcpp::List run_stardist_cpp(
    Rcpp::NumericVector tensor_vec,
    int in_w,
    int in_h,
    double orig_w,
    double orig_h,
    double prob_threshold = 0.5,
    double nms_threshold = 0.3,
    std::string model_path = "",
    std::string lib_path = "",
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  OrtTypeInfo* in_ti = NULL;
  (void)ort->SessionGetInputTypeInfo(session, 0, &in_ti);
  const OrtTensorTypeAndShapeInfo* in_si = NULL;
  (void)ort->CastTypeInfoToTensorInfo(in_ti, &in_si);
  size_t in_dims_cnt = 0;
  (void)ort->GetDimensionsCount(in_si, &in_dims_cnt);
  std::vector<int64_t> in_dims(in_dims_cnt);
  (void)ort->GetDimensions(in_si, in_dims.data(), in_dims_cnt);
  ort->ReleaseTypeInfo(in_ti);

  size_t total_floats = (size_t)tensor_vec.size();

  std::vector<float> input_vals;
  int64_t in_shape[4];

  if (in_dims_cnt == 4 && in_dims[3] == 1) {
    // NHWC format: [1, H, W, 1]
    in_shape[0] = 1; in_shape[1] = in_h; in_shape[2] = in_w; in_shape[3] = 1;
    size_t plane = (size_t)in_h * in_w;
    input_vals.resize(plane);
    if (total_floats >= 3 * plane) {
      for (size_t i = 0; i < plane; ++i) {
        input_vals[i] = 0.2989f * (float)tensor_vec[i] +
                        0.5870f * (float)tensor_vec[plane + i] +
                        0.1140f * (float)tensor_vec[2 * plane + i];
      }
    } else {
      for (size_t i = 0; i < plane && i < total_floats; ++i) {
        input_vals[i] = (float)tensor_vec[i];
      }
    }
  } else {
    // NCHW format: [1, C, H, W]
    int in_c = (total_floats == (size_t)in_h * in_w) ? 1 : 3;
    in_shape[0] = 1; in_shape[1] = in_c; in_shape[2] = in_h; in_shape[3] = in_w;
    input_vals.resize(total_floats);
    for (size_t i = 0; i < total_floats; ++i) input_vals[i] = (float)tensor_vec[i];
  }

  OrtValue* in_tensor = NULL;
  (void)ort->CreateTensorWithDataAsOrtValue(
    memory_info, input_vals.data(), input_vals.size() * sizeof(float),
    in_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_tensor
  );

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  size_t num_outputs = out_names.size();
  std::vector<OrtValue*> out_tensors(num_outputs, NULL);
  OrtStatus* status = ort->Run(
    session, run_options, in_names.data(), (const OrtValue* const*)&in_tensor,
    1, out_names.data(), num_outputs, out_tensors.data()
  );



  if (status != NULL) {
    std::string msg = ort->GetErrorMessage(status);
    ort->ReleaseStatus(status);
    ort->ReleaseRunOptions(run_options);
    ort->ReleaseValue(in_tensor);
    Rcpp::stop("StarDist run failed: " + msg);
  }

  int prob_idx = 0;
  int dist_idx = 1;
  int num_rays = 32;
  int out_h = in_h;
  int out_w = in_w;
  bool dist_is_nhwc = false;

  for (size_t i = 0; i < num_outputs; ++i) {
    OrtTensorTypeAndShapeInfo* si = NULL;
    (void)ort->GetTensorTypeAndShape(out_tensors[i], &si);
    size_t d_cnt = 0;
    (void)ort->GetDimensionsCount(si, &d_cnt);
    std::vector<int64_t> d(d_cnt);
    (void)ort->GetDimensions(si, d.data(), d_cnt);
    ort->ReleaseTensorTypeAndShapeInfo(si);

    std::string o_name = (i < out_names_str.size()) ? out_names_str[i] : "";
    std::transform(o_name.begin(), o_name.end(), o_name.begin(), ::tolower);

    if (o_name.find("prob") != std::string::npos || (d_cnt == 4 && d[3] == 1)) {
      prob_idx = (int)i;
      if (d_cnt == 4) {
        out_h = (int)d[1];
        out_w = (int)d[2];
      }
    } else if (o_name.find("dist") != std::string::npos || (d_cnt == 4 && d[3] > 1)) {
      dist_idx = (int)i;
      if (d_cnt == 4) {
        if (d[3] > 1) {
          num_rays = (int)d[3];
          dist_is_nhwc = true;
          out_h = (int)d[1];
          out_w = (int)d[2];
        } else if (d[1] > 1) {
          num_rays = (int)d[1];
          dist_is_nhwc = false;
          out_h = (int)d[2];
          out_w = (int)d[3];
        }
      }
    }
  }

  float* prob_data = NULL;
  (void)ort->GetTensorMutableData(out_tensors[prob_idx], (void**)&prob_data);

  float* dist_data = NULL;
  (void)ort->GetTensorMutableData(out_tensors[dist_idx], (void**)&dist_data);

  struct StarCand {
    float cx, cy;
    float prob;
    float x1, y1, x2, y2;
    std::vector<float> vx;
    std::vector<float> vy;
  };
  std::vector<StarCand> candidates;
  size_t out_plane = (size_t)out_h * out_w;
  float grid_x = (float)in_w / (float)out_w;
  float grid_y = (float)in_h / (float)out_h;

  std::vector<float> cos_phi(num_rays);
  std::vector<float> sin_phi(num_rays);
  for (int k = 0; k < num_rays; ++k) {
    float phi = 2.0f * 3.14159265f * (float)k / (float)num_rays;
    cos_phi[k] = std::cos(phi);
    sin_phi[k] = std::sin(phi);
  }

  for (int y = 0; y < out_h; ++y) {
    for (int x = 0; x < out_w; ++x) {
      float p = prob_data[(size_t)y * out_w + x];
      if (p >= (float)prob_threshold) {
        StarCand cand;
        cand.cx = ((float)x + 0.5f) * grid_x;
        cand.cy = ((float)y + 0.5f) * grid_y;
        cand.prob = p;
        cand.vx.resize(num_rays);
        cand.vy.resize(num_rays);

        float min_x = 1e9f, min_y = 1e9f, max_x = -1e9f, max_y = -1e9f;
        for (int k = 0; k < num_rays; ++k) {
          float r = dist_is_nhwc ?
            dist_data[((size_t)y * out_w + x) * (size_t)num_rays + k] :
            dist_data[(size_t)k * out_plane + (size_t)y * out_w + x];
          float px = cand.cx + r * cos_phi[k];
          float py = cand.cy + r * sin_phi[k];
          cand.vx[k] = px;
          cand.vy[k] = py;
          if (px < min_x) min_x = px;
          if (py < min_y) min_y = py;
          if (px > max_x) max_x = px;
          if (py > max_y) max_y = py;
        }
        cand.x1 = min_x; cand.y1 = min_y; cand.x2 = max_x; cand.y2 = max_y;
        candidates.push_back(cand);
      }
    }
  }



  std::sort(candidates.begin(), candidates.end(), [](const StarCand& a, const StarCand& b) {
    return a.prob > b.prob;
  });

  std::vector<bool> suppressed(candidates.size(), false);
  std::vector<StarCand> kept;

  for (size_t i = 0; i < candidates.size(); ++i) {
    if (suppressed[i]) continue;
    kept.push_back(candidates[i]);

    float area_i = std::max(0.0f, candidates[i].x2 - candidates[i].x1) *
                   std::max(0.0f, candidates[i].y2 - candidates[i].y1);

    for (size_t j = i + 1; j < candidates.size(); ++j) {
      if (suppressed[j]) continue;

      float inter_x1 = std::max(candidates[i].x1, candidates[j].x1);
      float inter_y1 = std::max(candidates[i].y1, candidates[j].y1);
      float inter_x2 = std::min(candidates[i].x2, candidates[j].x2);
      float inter_y2 = std::min(candidates[i].y2, candidates[j].y2);

      float inter_w = std::max(0.0f, inter_x2 - inter_x1);
      float inter_h = std::max(0.0f, inter_y2 - inter_y1);
      float inter_area = inter_w * inter_h;

      float area_j = std::max(0.0f, candidates[j].x2 - candidates[j].x1) *
                     std::max(0.0f, candidates[j].y2 - candidates[j].y1);

      float union_area = area_i + area_j - inter_area;
      float iou = (union_area > 0.0f) ? (inter_area / union_area) : 0.0f;

      if (iou > (float)nms_threshold) {
        suppressed[j] = true;
      }
    }
  }

  double sx = orig_w / (double)in_w;
  double sy = orig_h / (double)in_h;
  int num_kept = (int)kept.size();

  Rcpp::NumericMatrix boxes(num_kept, 4);
  Rcpp::NumericVector scores(num_kept);
  Rcpp::NumericVector centers_x(num_kept);
  Rcpp::NumericVector centers_y(num_kept);
  Rcpp::List polys_x(num_kept);
  Rcpp::List polys_y(num_kept);

  for (int i = 0; i < num_kept; ++i) {
    boxes(i, 0) = std::max(0.0, std::min(orig_w, kept[i].x1 * sx));
    boxes(i, 1) = std::max(0.0, std::min(orig_h, kept[i].y1 * sy));
    boxes(i, 2) = std::max(0.0, std::min(orig_w, kept[i].x2 * sx));
    boxes(i, 3) = std::max(0.0, std::min(orig_h, kept[i].y2 * sy));
    scores[i] = kept[i].prob;
    centers_x[i] = kept[i].cx * sx;
    centers_y[i] = kept[i].cy * sy;

    Rcpp::NumericVector px_vec(num_rays);
    Rcpp::NumericVector py_vec(num_rays);
    for (int k = 0; k < num_rays; ++k) {
      px_vec[k] = kept[i].vx[k] * sx;
      py_vec[k] = kept[i].vy[k] * sy;
    }
    polys_x[i] = px_vec;
    polys_y[i] = py_vec;
  }

  int ow = (int)std::round(orig_w);
  int oh = (int)std::round(orig_h);
  Rcpp::IntegerMatrix labels(ow, oh);
  Rcpp::LogicalMatrix mask(ow, oh);

  for (int i = 0; i < num_kept; ++i) {
    int bx1 = std::max(0, (int)std::floor(boxes(i, 0)));
    int by1 = std::max(0, (int)std::floor(boxes(i, 1)));
    int bx2 = std::min(ow - 1, (int)std::ceil(boxes(i, 2)));
    int by2 = std::min(oh - 1, (int)std::ceil(boxes(i, 3)));

    Rcpp::NumericVector px_vec = polys_x[i];
    Rcpp::NumericVector py_vec = polys_y[i];

    for (int y = by1; y <= by2; ++y) {
      double py = (double)y;
      for (int x = bx1; x <= bx2; ++x) {
        double px = (double)x;
        bool inside = false;
        for (int j = 0, k = num_rays - 1; j < num_rays; k = j++) {
          if (((py_vec[j] > py) != (py_vec[k] > py)) &&
              (px < (px_vec[k] - px_vec[j]) * (py - py_vec[j]) / (py_vec[k] - py_vec[j]) + px_vec[j])) {
            inside = !inside;
          }
        }
        if (inside) {
          labels(x, y) = i + 1;
          mask(x, y) = true;
        }
      }
    }
  }

  for (size_t i = 0; i < num_outputs; ++i) ort->ReleaseValue(out_tensors[i]);
  ort->ReleaseValue(in_tensor);
  ort->ReleaseRunOptions(run_options);

  return Rcpp::List::create(
    Rcpp::Named("boxes") = boxes,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("centers_x") = centers_x,
    Rcpp::Named("centers_y") = centers_y,
    Rcpp::Named("polygons_x") = polys_x,
    Rcpp::Named("polygons_y") = polys_y,
    Rcpp::Named("labels") = labels,
    Rcpp::Named("mask") = mask
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector run_super_resolution_cpp(
    Rcpp::NumericVector tensor_vec,
    int in_w,
    int in_h,
    int scale = 4,
    int tile_size = 256,
    int tile_pad = 16,
    std::string model_path = "",
    std::string lib_path = "",
    int num_threads = 0,
    bool use_gpu = false,
    int device_id = -1
) {
  const OrtApi* ort = get_ort_api(lib_path);
  CachedSession cs = get_or_create_cached_session(ort, model_path, num_threads, use_gpu, device_id);
  OrtSession* session = cs.session;
  OrtMemoryInfo* memory_info = cs.mem_info;

  std::vector<std::string> in_names_str, out_names_str;
  std::vector<const char*> in_names, out_names;
  get_session_io_names(ort, session, in_names_str, in_names, out_names_str, out_names);

  int out_w = in_w * scale;
  int out_h = in_h * scale;
  size_t out_plane = (size_t)out_h * out_w;

  std::vector<float> full_out(3 * out_plane, 0.0f);
  std::vector<float> weight_map(out_plane, 0.0f);

  OrtRunOptions* run_options = NULL;
  (void)ort->CreateRunOptions(&run_options);
  (void)ort->RunOptionsSetRunLogSeverityLevel(run_options, ORT_LOGGING_LEVEL_ERROR);

  size_t in_plane = (size_t)in_h * in_w;
  const double* tensor_data = tensor_vec.begin();

  if (in_w <= tile_size && in_h <= tile_size) {
    int64_t shape[4] = {1, 3, in_h, in_w};
    size_t floats = 3 * in_plane;
    std::vector<float> in_vals(floats);
    for (size_t i = 0; i < floats; ++i) in_vals[i] = (float)tensor_data[i];

    OrtValue* in_val = NULL;
    (void)ort->CreateTensorWithDataAsOrtValue(
      memory_info, in_vals.data(), floats * sizeof(float),
      shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_val
    );

    OrtValue* out_val = NULL;
    OrtStatus* st = ort->Run(session, run_options, in_names.data(), (const OrtValue* const*)&in_val, 1, out_names.data(), 1, &out_val);
    if (st == NULL) {
      float* out_f = NULL;
      (void)ort->GetTensorMutableData(out_val, (void**)&out_f);
      for (size_t i = 0; i < 3 * out_plane; ++i) {
        full_out[i] = std::max(0.0f, std::min(1.0f, out_f[i]));
      }
      ort->ReleaseValue(out_val);
    } else {
      std::string msg = ort->GetErrorMessage(st);
      ort->ReleaseStatus(st);
      ort->ReleaseValue(in_val);
      ort->ReleaseRunOptions(run_options);
      Rcpp::stop("Super resolution run failed: " + msg);
    }
    ort->ReleaseValue(in_val);
  } else {
    int step = tile_size - 2 * tile_pad;
    if (step < 16) step = 16;

    for (int ty = 0; ty < in_h; ty += step) {
      for (int tx = 0; tx < in_w; tx += step) {
        int x1 = std::max(0, tx - tile_pad);
        int y1 = std::max(0, ty - tile_pad);
        int x2 = std::min(in_w, tx + step + tile_pad);
        int y2 = std::min(in_h, ty + step + tile_pad);

        int cur_w = x2 - x1;
        int cur_h = y2 - y1;
        size_t cur_plane = (size_t)cur_w * cur_h;

        int64_t cur_shape[4] = {1, 3, cur_h, cur_w};
        size_t cur_floats = 3 * cur_plane;
        std::vector<float> cur_in(cur_floats);

        for (int c = 0; c < 3; ++c) {
          for (int py = 0; py < cur_h; ++py) {
            for (int px = 0; px < cur_w; ++px) {
              cur_in[(size_t)c * cur_plane + (size_t)py * cur_w + px] =
                (float)tensor_data[(size_t)c * in_plane + (size_t)(y1 + py) * in_w + (x1 + px)];
            }
          }
        }

        OrtValue* in_val = NULL;
        (void)ort->CreateTensorWithDataAsOrtValue(
          memory_info, cur_in.data(), cur_floats * sizeof(float),
          cur_shape, 4, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in_val
        );

        OrtValue* out_val = NULL;
        OrtStatus* st = ort->Run(session, run_options, in_names.data(), (const OrtValue* const*)&in_val, 1, out_names.data(), 1, &out_val);
        if (st == NULL) {
          float* out_f = NULL;
          (void)ort->GetTensorMutableData(out_val, (void**)&out_f);

          int cur_out_w = cur_w * scale;
          int cur_out_h = cur_h * scale;
          size_t cur_out_plane = (size_t)cur_out_w * cur_out_h;

          int dst_x1 = x1 * scale;
          int dst_y1 = y1 * scale;

          for (int py = 0; py < cur_out_h; ++py) {
            int gy = dst_y1 + py;
            if (gy >= out_h) continue;

            float wy = 1.0f;
            if (py < tile_pad * scale && dst_y1 > 0) wy = (float)py / (float)(tile_pad * scale);
            else if (py > cur_out_h - tile_pad * scale && dst_y1 + cur_out_h < out_h)
              wy = (float)(cur_out_h - py) / (float)(tile_pad * scale);

            for (int px = 0; px < cur_out_w; ++px) {
              int gx = dst_x1 + px;
              if (gx >= out_w) continue;

              float wx = 1.0f;
              if (px < tile_pad * scale && dst_x1 > 0) wx = (float)px / (float)(tile_pad * scale);
              else if (px > cur_out_w - tile_pad * scale && dst_x1 + cur_out_w < out_w)
                wx = (float)(cur_out_w - px) / (float)(tile_pad * scale);

              float w = wx * wy;
              if (w <= 0.0f) w = 1e-4f;

              size_t dst_idx = (size_t)gy * out_w + gx;
              size_t src_idx = (size_t)py * cur_out_w + px;

              weight_map[dst_idx] += w;
              for (int c = 0; c < 3; ++c) {
                full_out[(size_t)c * out_plane + dst_idx] += w * out_f[(size_t)c * cur_out_plane + src_idx];
              }
            }
          }
          ort->ReleaseValue(out_val);
        } else {
          ort->ReleaseStatus(st);
        }
        ort->ReleaseValue(in_val);
      }
    }

    for (size_t i = 0; i < out_plane; ++i) {
      float total_w = weight_map[i];
      if (total_w > 0.0f) {
        for (int c = 0; c < 3; ++c) {
          full_out[(size_t)c * out_plane + i] /= total_w;
          full_out[(size_t)c * out_plane + i] = std::max(0.0f, std::min(1.0f, full_out[(size_t)c * out_plane + i]));
        }
      }
    }
  }

  ort->ReleaseRunOptions(run_options);

  Rcpp::NumericVector res(3 * out_plane);
  for (size_t i = 0; i < 3 * out_plane; ++i) res[i] = (double)full_out[i];
  return res;
}

