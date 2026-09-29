#include <Rcpp.h>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include "stb_image_write.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <mfapi.h>
#include <mfidl.h>
#include <mfreadwrite.h>
#include <mferror.h>

static bool g_mf_started = false;
static const wchar_t* PLIMAN_WINDOW_CLASS = L"PlimanNativeVideoWindow";
static bool g_window_class_registered = false;
static const GUID PLIMAN_GUID_NULL = { 0, 0, 0, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static void ensure_mf_started() {
  if (!g_mf_started) {
    CoInitializeEx(NULL, COINIT_MULTITHREADED);
    MFStartup(MF_VERSION);
    g_mf_started = true;
  }
}

struct NativeCamera {
  std::thread worker_thread;
  std::atomic<bool> is_running{false};
  std::atomic<bool> is_open{false};
  std::atomic<bool> has_new_frame{false};
  std::mutex frame_mutex;
  std::vector<unsigned char> latest_frame_bgr32;
  UINT32 width = 0;
  UINT32 height = 0;
  INT32 stride = 0;

  ~NativeCamera() {
    close();
  }

  void close() {
    is_running = false;
    if (worker_thread.joinable()) {
      worker_thread.join();
    }
    is_open = false;
  }
};

static void camera_worker_loop(NativeCamera* cam, int cam_id, int target_w, int target_h) {
  HRESULT hr_co = CoInitializeEx(NULL, COINIT_MULTITHREADED);
  HRESULT hr_mf = MFStartup(MF_VERSION);

  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 1);
  if (FAILED(hr)) {
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    cam->is_running = false;
    return;
  }
  pAttributes->SetGUID(MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE, MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE_VIDCAP_GUID);
  IMFActivate** ppDevices = nullptr;
  UINT32 count = 0;
  hr = MFEnumDeviceSources(pAttributes, &ppDevices, &count);
  pAttributes->Release();
  if (FAILED(hr) || count == 0 || (UINT32)cam_id >= count) {
    if (ppDevices) {
      for (UINT32 i = 0; i < count; ++i) ppDevices[i]->Release();
      CoTaskMemFree(ppDevices);
    }
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    cam->is_running = false;
    return;
  }

  IMFActivate* pActivate = ppDevices[cam_id];
  for (UINT32 i = 0; i < count; ++i) {
    if (i != (UINT32)cam_id) ppDevices[i]->Release();
  }
  CoTaskMemFree(ppDevices);

  IMFMediaSource* pSource = nullptr;
  hr = pActivate->ActivateObject(IID_PPV_ARGS(&pSource));
  if (FAILED(hr)) {
    pActivate->Release();
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    cam->is_running = false;
    return;
  }

  IMFAttributes* pReaderAttributes = nullptr;
  MFCreateAttributes(&pReaderAttributes, 2);
  pReaderAttributes->SetUINT32(MF_SOURCE_READER_ENABLE_VIDEO_PROCESSING, TRUE);
  pReaderAttributes->SetUINT32(MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS, FALSE);

  IMFSourceReader* pReader = nullptr;
  hr = MFCreateSourceReaderFromMediaSource(pSource, pReaderAttributes, &pReader);
  pReaderAttributes->Release();
  if (FAILED(hr)) {
    pSource->Shutdown();
    pSource->Release();
    pActivate->Release();
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    cam->is_running = false;
    return;
  }

  // If target resolution requested, find the closest matching native resolution first
  if (target_w > 0 && target_h > 0) {
    DWORD typeIdx = 0;
    while (true) {
      IMFMediaType* pNative = nullptr;
      HRESULT hr_nt = pReader->GetNativeMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, typeIdx, &pNative);
      if (FAILED(hr_nt)) break;
      UINT32 nw = 0, nh = 0;
      MFGetAttributeSize(pNative, MF_MT_FRAME_SIZE, &nw, &nh);
      if ((int)nw == target_w && (int)nh == target_h) {
        pReader->SetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, NULL, pNative);
        pNative->Release();
        break;
      }
      pNative->Release();
      typeIdx++;
    }
  }

  IMFMediaType* pMediaType = nullptr;
  MFCreateMediaType(&pMediaType);
  pMediaType->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
  pMediaType->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_RGB32);
  if (target_w > 0 && target_h > 0) {
    MFSetAttributeSize(pMediaType, MF_MT_FRAME_SIZE, target_w, target_h);
  }
  hr = pReader->SetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, NULL, pMediaType);
  pMediaType->Release();
  if (FAILED(hr)) {
    pReader->Release();
    pSource->Shutdown();
    pSource->Release();
    pActivate->Release();
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    cam->is_running = false;
    return;
  }

  IMFMediaType* pActualType = nullptr;
  hr = pReader->GetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, &pActualType);
  UINT32 width = 0, height = 0;
  INT32 stride = 0;
  if (SUCCEEDED(hr) && pActualType) {
    MFGetAttributeSize(pActualType, MF_MT_FRAME_SIZE, &width, &height);
    pActualType->GetUINT32(MF_MT_DEFAULT_STRIDE, (UINT32*)&stride);
    pActualType->Release();
  }
  if (stride == 0) stride = (INT32)(width * 4);

  {
    std::lock_guard<std::mutex> lock(cam->frame_mutex);
    cam->width = width;
    cam->height = height;
    cam->stride = stride;
    cam->latest_frame_bgr32.resize((size_t)width * height * 4);
    cam->is_open = true;
  }

  // Camera streaming loop in dedicated thread
  while (cam->is_running) {
    DWORD streamIndex = 0, flags = 0;
    LONGLONG timestamp = 0;
    IMFSample* pSample = nullptr;
    HRESULT hr_rs = pReader->ReadSample(
        (DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM,
        0,
        &streamIndex,
        &flags,
        &timestamp,
        &pSample
    );
    if (FAILED(hr_rs) || (flags & MF_SOURCE_READERF_ENDOFSTREAM)) {
      break;
    }
    if (pSample) {
      IMFMediaBuffer* pBuffer = nullptr;
      if (SUCCEEDED(pSample->ConvertToContiguousBuffer(&pBuffer)) && pBuffer) {
        BYTE* pData = nullptr;
        DWORD maxLen = 0, curLen = 0;
        if (SUCCEEDED(pBuffer->Lock(&pData, &maxLen, &curLen)) && pData) {
          {
            std::lock_guard<std::mutex> lock(cam->frame_mutex);
            size_t copy_size = std::min((size_t)curLen, cam->latest_frame_bgr32.size());
            memcpy(cam->latest_frame_bgr32.data(), pData, copy_size);
            cam->has_new_frame = true;
          }
          pBuffer->Unlock();
        }
        pBuffer->Release();
      }
      pSample->Release();
    } else {
      Sleep(2);
    }
  }

  pReader->Release();
  pSource->Shutdown();
  pSource->Release();
  pActivate->Release();
  if (SUCCEEDED(hr_mf)) MFShutdown();
  if (SUCCEEDED(hr_co)) CoUninitialize();
  cam->is_open = false;
}

struct NativeWindow {
  HWND hwnd = NULL;
  HDC hdc = NULL;
  bool should_close = false;
  bool is_paused = false;
  bool is_fullscreen = false;
  RECT normal_rect = {0, 0, 0, 0};
  int last_key = 0;
  int frame_w = 0;
  int frame_h = 0;
  std::vector<unsigned char> bgr_buffer;

  ~NativeWindow() {
    close();
  }

  void close() {
    if (hdc && hwnd) {
      ReleaseDC(hwnd, hdc);
      hdc = NULL;
    }
    if (hwnd) {
      DestroyWindow(hwnd);
      hwnd = NULL;
    }
  }
};

static LRESULT CALLBACK PlimanWindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
  NativeWindow* ctx = (NativeWindow*)GetWindowLongPtr(hwnd, GWLP_USERDATA);
  switch (uMsg) {
    case WM_CLOSE:
      if (ctx) ctx->should_close = true;
      return 0;
    case WM_DESTROY:
      if (ctx) ctx->should_close = true;
      return 0;
    case WM_KEYDOWN:
      if (ctx) {
        ctx->last_key = (int)wParam;
        if (wParam == VK_ESCAPE || wParam == 'Q' || wParam == 'q') {
          ctx->should_close = true;
        } else if (wParam == VK_SPACE) {
          ctx->is_paused = !ctx->is_paused;
        } else if (wParam == 'F' || wParam == 'f') {
          if (!ctx->is_fullscreen) {
            GetWindowRect(hwnd, &ctx->normal_rect);
            SetWindowLongPtr(hwnd, GWL_STYLE, WS_POPUP | WS_VISIBLE);
            ShowWindow(hwnd, SW_MAXIMIZE);
            ctx->is_fullscreen = true;
          } else {
            SetWindowLongPtr(hwnd, GWL_STYLE, WS_OVERLAPPEDWINDOW | WS_VISIBLE);
            ShowWindow(hwnd, SW_RESTORE);
            SetWindowPos(hwnd, NULL, ctx->normal_rect.left, ctx->normal_rect.top,
                         ctx->normal_rect.right - ctx->normal_rect.left,
                         ctx->normal_rect.bottom - ctx->normal_rect.top,
                         SWP_NOZORDER | SWP_FRAMECHANGED);
            ctx->is_fullscreen = false;
          }
        }
      }
      return 0;
    case WM_ERASEBKGND:
      return 1;
  }
  return DefWindowProcW(hwnd, uMsg, wParam, lParam);
}

static void camera_finalizer(SEXP ptr) {
  NativeCamera* cam = (NativeCamera*)R_ExternalPtrAddr(ptr);
  if (cam) {
    delete cam;
    R_ClearExternalPtr(ptr);
  }
}

static void window_finalizer(SEXP ptr) {
  NativeWindow* win = (NativeWindow*)R_ExternalPtrAddr(ptr);
  if (win) {
    delete win;
    R_ClearExternalPtr(ptr);
  }
}
#endif

// [[Rcpp::export]]
bool has_native_camera_cpp() {
#ifdef _WIN32
  ensure_mf_started();
  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 1);
  if (FAILED(hr)) return false;
  pAttributes->SetGUID(MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE, MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE_VIDCAP_GUID);
  IMFActivate** ppDevices = nullptr;
  UINT32 count = 0;
  hr = MFEnumDeviceSources(pAttributes, &ppDevices, &count);
  pAttributes->Release();
  if (FAILED(hr)) return false;
  for (UINT32 i = 0; i < count; ++i) ppDevices[i]->Release();
  CoTaskMemFree(ppDevices);
  return count > 0;
#else
  return false;
#endif
}

// [[Rcpp::export]]
Rcpp::List list_native_cameras_cpp() {
#ifdef _WIN32
  ensure_mf_started();
  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 1);
  if (FAILED(hr)) return Rcpp::List::create();
  pAttributes->SetGUID(MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE, MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE_VIDCAP_GUID);
  IMFActivate** ppDevices = nullptr;
  UINT32 count = 0;
  hr = MFEnumDeviceSources(pAttributes, &ppDevices, &count);
  pAttributes->Release();
  if (FAILED(hr) || count == 0) return Rcpp::List::create();

  Rcpp::CharacterVector names(count);
  Rcpp::IntegerVector ids(count);

  for (UINT32 i = 0; i < count; ++i) {
    ids[i] = (int)i;
    WCHAR* friendlyName = nullptr;
    UINT32 nameLen = 0;
    ppDevices[i]->GetAllocatedString(MF_DEVSOURCE_ATTRIBUTE_FRIENDLY_NAME, &friendlyName, &nameLen);
    if (friendlyName) {
      char mbName[512] = {0};
      WideCharToMultiByte(CP_UTF8, 0, friendlyName, -1, mbName, sizeof(mbName) - 1, NULL, NULL);
      names[i] = mbName;
      CoTaskMemFree(friendlyName);
    } else {
      names[i] = "Camera " + std::to_string(i);
    }
    ppDevices[i]->Release();
  }
  CoTaskMemFree(ppDevices);
  return Rcpp::List::create(Rcpp::Named("id") = ids, Rcpp::Named("name") = names);
#else
  return Rcpp::List::create();
#endif
}

// [[Rcpp::export]]
Rcpp::DataFrame list_camera_formats_cpp(int cam_id = 0) {
#ifdef _WIN32
  HRESULT hr_co = CoInitializeEx(NULL, COINIT_MULTITHREADED);
  HRESULT hr_mf = MFStartup(MF_VERSION);

  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 1);
  if (FAILED(hr)) {
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    return Rcpp::DataFrame::create();
  }
  pAttributes->SetGUID(MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE, MF_DEVSOURCE_ATTRIBUTE_SOURCE_TYPE_VIDCAP_GUID);
  IMFActivate** ppDevices = nullptr;
  UINT32 count = 0;
  hr = MFEnumDeviceSources(pAttributes, &ppDevices, &count);
  pAttributes->Release();
  if (FAILED(hr) || count == 0 || (UINT32)cam_id >= count) {
    if (ppDevices) {
      for (UINT32 i = 0; i < count; ++i) ppDevices[i]->Release();
      CoTaskMemFree(ppDevices);
    }
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    return Rcpp::DataFrame::create();
  }

  IMFActivate* pActivate = ppDevices[cam_id];
  for (UINT32 i = 0; i < count; ++i) {
    if (i != (UINT32)cam_id) ppDevices[i]->Release();
  }
  CoTaskMemFree(ppDevices);

  IMFMediaSource* pSource = nullptr;
  hr = pActivate->ActivateObject(IID_PPV_ARGS(&pSource));
  pActivate->Release();
  if (FAILED(hr)) {
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    return Rcpp::DataFrame::create();
  }

  IMFSourceReader* pReader = nullptr;
  hr = MFCreateSourceReaderFromMediaSource(pSource, NULL, &pReader);
  if (FAILED(hr)) {
    pSource->Shutdown();
    pSource->Release();
    if (SUCCEEDED(hr_mf)) MFShutdown();
    if (SUCCEEDED(hr_co)) CoUninitialize();
    return Rcpp::DataFrame::create();
  }

  std::vector<int> widths, heights;
  std::vector<double> fps_list;
  std::vector<std::string> formats;

  DWORD typeIndex = 0;
  while (true) {
    IMFMediaType* pNativeType = nullptr;
    hr = pReader->GetNativeMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, typeIndex, &pNativeType);
    if (FAILED(hr)) break;

    UINT32 w = 0, h = 0;
    MFGetAttributeSize(pNativeType, MF_MT_FRAME_SIZE, &w, &h);

    UINT32 num = 0, den = 0;
    MFGetAttributeRatio(pNativeType, MF_MT_FRAME_RATE, &num, &den);
    double fps = (den > 0) ? ((double)num / den) : 0.0;

    GUID subtype = {0};
    pNativeType->GetGUID(MF_MT_SUBTYPE, &subtype);
    std::string fmt = "Other";
    if (IsEqualGUID(subtype, MFVideoFormat_MJPG)) fmt = "MJPG";
    else if (IsEqualGUID(subtype, MFVideoFormat_YUY2)) fmt = "YUY2";
    else if (IsEqualGUID(subtype, MFVideoFormat_NV12)) fmt = "NV12";
    else if (IsEqualGUID(subtype, MFVideoFormat_RGB24)) fmt = "RGB24";
    else if (IsEqualGUID(subtype, MFVideoFormat_RGB32)) fmt = "RGB32";

    widths.push_back((int)w);
    heights.push_back((int)h);
    fps_list.push_back(fps);
    formats.push_back(fmt);

    pNativeType->Release();
    typeIndex++;
  }

  pReader->Release();
  pSource->Shutdown();
  pSource->Release();
  if (SUCCEEDED(hr_mf)) MFShutdown();
  if (SUCCEEDED(hr_co)) CoUninitialize();

  return Rcpp::DataFrame::create(
    Rcpp::Named("width") = widths,
    Rcpp::Named("height") = heights,
    Rcpp::Named("fps") = fps_list,
    Rcpp::Named("format") = formats,
    Rcpp::Named("stringsAsFactors") = false
  );
#else
  return Rcpp::DataFrame::create();
#endif
}

// [[Rcpp::export]]
SEXP open_native_camera_cpp(int cam_id = 0, int target_w = 0, int target_h = 0) {
#ifdef _WIN32
  NativeCamera* cam = new NativeCamera();
  cam->is_running = true;
  cam->is_open = false;
  cam->has_new_frame = false;
  cam->worker_thread = std::thread(camera_worker_loop, cam, cam_id, target_w, target_h);

  // Wait for worker thread to initialize camera and receive first frame (up to 3 seconds)
  int wait_ms = 0;
  while ((!cam->is_open || !cam->has_new_frame) && cam->is_running && wait_ms < 3000) {
    Sleep(15);
    wait_ms += 15;
  }

  if (!cam->is_open) {
    cam->close();
    delete cam;
    Rcpp::stop("Failed to open camera %d or camera did not start streaming.", cam_id);
  }

  SEXP ptr = PROTECT(R_MakeExternalPtr(cam, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, camera_finalizer, (Rboolean)TRUE);
  UNPROTECT(1);
  return ptr;
#else
  Rcpp::stop("Native camera driver is only supported on Windows.");
  return R_NilValue;
#endif
}

// [[Rcpp::export]]
Rcpp::List get_camera_dims_cpp(SEXP cam_ptr) {
#ifdef _WIN32
  if (TYPEOF(cam_ptr) != EXTPTRSXP) Rcpp::stop("Invalid camera pointer");
  NativeCamera* cam = (NativeCamera*)R_ExternalPtrAddr(cam_ptr);
  if (!cam || !cam->is_open) Rcpp::stop("Camera is not open");
  return Rcpp::List::create(Rcpp::Named("width") = (int)cam->width, Rcpp::Named("height") = (int)cam->height);
#else
  return Rcpp::List::create(Rcpp::Named("width") = 0, Rcpp::Named("height") = 0);
#endif
}

// [[Rcpp::export]]
bool grab_native_frame_cpp(SEXP cam_ptr, Rcpp::RawVector out_bm, int crop_x1 = 0, int crop_y1 = 0, int crop_w = 0, int crop_h = 0) {
#ifdef _WIN32
  if (TYPEOF(cam_ptr) != EXTPTRSXP) return false;
  NativeCamera* cam = (NativeCamera*)R_ExternalPtrAddr(cam_ptr);
  if (!cam || !cam->is_open || !cam->is_running) return false;

  // Wait for a fresh frame (up to 1000ms)
  int wait_ms = 0;
  while (!cam->has_new_frame && cam->is_running && wait_ms < 1000) {
    Sleep(2);
    wait_ms += 2;
  }
  if (!cam->is_running) return false;

  int w = (int)cam->width;
  int h = (int)cam->height;
  int stride = std::abs((int)cam->stride);
  bool flip_vertical = (cam->stride < 0);
  unsigned char* dst = RAW(out_bm);

  int out_w = (crop_w > 0) ? crop_w : w;
  int out_h = (crop_h > 0) ? crop_h : h;
  int start_x = (crop_w > 0) ? std::max(0, std::min(w - 1, crop_x1)) : 0;
  int start_y = (crop_h > 0) ? std::max(0, std::min(h - 1, crop_y1)) : 0;
  out_w = std::min(out_w, w - start_x);
  out_h = std::min(out_h, h - start_y);

  {
    std::lock_guard<std::mutex> lock(cam->frame_mutex);
    const unsigned char* pData = cam->latest_frame_bgr32.data();

    #pragma omp parallel for schedule(static) if(out_h > 100)
    for (int y = 0; y < out_h; ++y) {
      int full_y = start_y + y;
      int src_y = flip_vertical ? (h - 1 - full_y) : full_y;
      const unsigned char* src_row = pData + (size_t)src_y * stride;
      unsigned char* dst_row = dst + (size_t)y * ((size_t)out_w * 3);
      for (int x = 0; x < out_w; ++x) {
        int full_x = start_x + x;
        size_t src_idx = (size_t)full_x * 4;
        size_t dst_idx = (size_t)x * 3;
        dst_row[dst_idx]     = src_row[src_idx + 2]; // R
        dst_row[dst_idx + 1] = src_row[src_idx + 1]; // G
        dst_row[dst_idx + 2] = src_row[src_idx];     // B
      }
    }
    cam->has_new_frame = false;
  }

  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
void close_native_camera_cpp(SEXP cam_ptr) {
#ifdef _WIN32
  if (TYPEOF(cam_ptr) == EXTPTRSXP) {
    NativeCamera* cam = (NativeCamera*)R_ExternalPtrAddr(cam_ptr);
    if (cam) {
      cam->close();
    }
  }
#endif
}

// [[Rcpp::export]]
SEXP create_native_window_cpp(std::string title = "pliman - YOLO Detection", int width = 640, int height = 480, bool fullscreen = false) {
#ifdef _WIN32
  HINSTANCE hInstance = GetModuleHandle(NULL);
  if (!g_window_class_registered) {
    WNDCLASSEXW wc = {0};
    wc.cbSize = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc = PlimanWindowProc;
    wc.hInstance = hInstance;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
    wc.lpszClassName = PLIMAN_WINDOW_CLASS;
    RegisterClassExW(&wc);
    g_window_class_registered = true;
  }

  std::wstring wtitle(title.begin(), title.end());
  DWORD style = fullscreen ? (WS_POPUP | WS_VISIBLE) : (WS_OVERLAPPEDWINDOW | WS_VISIBLE);

  // Ensure default window fits comfortably within desktop work area (never overflowing screen)
  RECT work_area = {0, 0, 1920, 1080};
  SystemParametersInfoW(SPI_GETWORKAREA, 0, &work_area, 0);
  int max_client_w = (int)((work_area.right - work_area.left) * 0.85);
  int max_client_h = (int)((work_area.bottom - work_area.top) * 0.85);

  int client_w = width;
  int client_h = height;
  if (!fullscreen && (client_w > max_client_w || client_h > max_client_h)) {
    float scale = std::min((float)max_client_w / (float)client_w, (float)max_client_h / (float)client_h);
    client_w = std::max(320, (int)std::round(client_w * scale));
    client_h = std::max(240, (int)std::round(client_h * scale));
  }

  RECT wr = {0, 0, client_w, client_h};
  if (!fullscreen) {
    AdjustWindowRect(&wr, WS_OVERLAPPEDWINDOW, FALSE);
  }

  int win_w = wr.right - wr.left;
  int win_h = wr.bottom - wr.top;
  int win_x = fullscreen ? 0 : (work_area.left + ((work_area.right - work_area.left) - win_w) / 2);
  int win_y = fullscreen ? 0 : (work_area.top + ((work_area.bottom - work_area.top) - win_h) / 2);

  HWND hwnd = CreateWindowExW(
      0, PLIMAN_WINDOW_CLASS, wtitle.c_str(),
      style,
      win_x, win_y, win_w, win_h,
      NULL, NULL, hInstance, NULL
  );
  if (!hwnd) Rcpp::stop("Failed to create native Win32 window");

  if (fullscreen) {
    ShowWindow(hwnd, SW_MAXIMIZE);
  } else {
    ShowWindow(hwnd, SW_SHOW);
  }
  UpdateWindow(hwnd);

  HDC hdc = GetDC(hwnd);
  NativeWindow* win = new NativeWindow();
  win->hwnd = hwnd;
  win->hdc = hdc;
  win->should_close = false;
  win->is_paused = false;
  win->is_fullscreen = fullscreen;
  win->frame_w = width;
  win->frame_h = height;
  win->bgr_buffer.resize((size_t)width * height * 3);

  SetWindowLongPtr(hwnd, GWLP_USERDATA, (LONG_PTR)win);

  SEXP ptr = PROTECT(R_MakeExternalPtr(win, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, window_finalizer, (Rboolean)TRUE);
  UNPROTECT(1);
  return ptr;
#else
  Rcpp::stop("Native window display is only supported on Windows.");
  return R_NilValue;
#endif
}

// [[Rcpp::export]]
void show_native_frame_cpp(SEXP win_ptr, Rcpp::RawVector bm, bool preserve_aspect = true) {
#ifdef _WIN32
  if (TYPEOF(win_ptr) != EXTPTRSXP) return;
  NativeWindow* win = (NativeWindow*)R_ExternalPtrAddr(win_ptr);
  if (!win || !win->hwnd || !win->hdc) return;

  SEXP dim_attr = Rf_getAttrib(bm, R_DimSymbol);
  if (dim_attr == R_NilValue || Rf_length(dim_attr) < 3) return;
  int* dims = INTEGER(dim_attr);
  int w = dims[1];
  int h = dims[2];

  size_t buf_size = (size_t)w * h * 3;
  if (win->bgr_buffer.size() < buf_size) {
    win->bgr_buffer.resize(buf_size);
  }

  const unsigned char* p_rgb = RAW(bm);
  unsigned char* p_bgr = win->bgr_buffer.data();

  // Convert RGB to BGR for StretchDIBits (BI_RGB expects BGR)
  #pragma omp parallel for schedule(static) if(h > 100)
  for (int y = 0; y < h; ++y) {
    const unsigned char* src_row = p_rgb + (size_t)y * ((size_t)w * 3);
    unsigned char* dst_row = p_bgr + (size_t)y * ((size_t)w * 3);
    for (int x = 0; x < w; ++x) {
      size_t idx = (size_t)x * 3;
      dst_row[idx]     = src_row[idx + 2]; // B
      dst_row[idx + 1] = src_row[idx + 1]; // G
      dst_row[idx + 2] = src_row[idx];     // R
    }
  }

  BITMAPINFO bmi = {0};
  bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
  bmi.bmiHeader.biWidth = w;
  bmi.bmiHeader.biHeight = -h; // Top-down
  bmi.bmiHeader.biPlanes = 1;
  bmi.bmiHeader.biBitCount = 24;
  bmi.bmiHeader.biCompression = BI_RGB;

  RECT rc;
  GetClientRect(win->hwnd, &rc);
  int win_w = rc.right - rc.left;
  int win_h = rc.bottom - rc.top;

  SetStretchBltMode(win->hdc, COLORONCOLOR);

  if (preserve_aspect && win_w > 0 && win_h > 0) {
    float aspect_img = (float)w / (float)h;
    float aspect_win = (float)win_w / (float)win_h;
    int dst_x = 0, dst_y = 0, dst_w = win_w, dst_h = win_h;

    if (aspect_win > aspect_img) {
      dst_w = (int)std::round(win_h * aspect_img);
      dst_x = (win_w - dst_w) / 2;
      // Clear left and right pillarbox margins if any
      if (dst_x > 0) {
        RECT r_left = {0, 0, dst_x, win_h};
        FillRect(win->hdc, &r_left, (HBRUSH)GetStockObject(BLACK_BRUSH));
        RECT r_right = {dst_x + dst_w, 0, win_w, win_h};
        FillRect(win->hdc, &r_right, (HBRUSH)GetStockObject(BLACK_BRUSH));
      }
    } else {
      dst_h = (int)std::round(win_w / aspect_img);
      dst_y = (win_h - dst_h) / 2;
      // Clear top and bottom letterbox margins if any
      if (dst_y > 0) {
        RECT r_top = {0, 0, win_w, dst_y};
        FillRect(win->hdc, &r_top, (HBRUSH)GetStockObject(BLACK_BRUSH));
        RECT r_bottom = {0, dst_y + dst_h, win_w, win_h};
        FillRect(win->hdc, &r_bottom, (HBRUSH)GetStockObject(BLACK_BRUSH));
      }
    }

    StretchDIBits(win->hdc, dst_x, dst_y, dst_w, dst_h,
                  0, 0, w, h,
                  p_bgr, &bmi, DIB_RGB_COLORS, SRCCOPY);
  } else {
    StretchDIBits(win->hdc, 0, 0, win_w, win_h,
                  0, 0, w, h,
                  p_bgr, &bmi, DIB_RGB_COLORS, SRCCOPY);
  }
#endif
}

// [[Rcpp::export]]
int poll_native_window_cpp(SEXP win_ptr) {
#ifdef _WIN32
  if (TYPEOF(win_ptr) != EXTPTRSXP) return 1;
  NativeWindow* win = (NativeWindow*)R_ExternalPtrAddr(win_ptr);
  if (!win || !win->hwnd) return 1;

  MSG msg;
  while (PeekMessageW(&msg, win->hwnd, 0, 0, PM_REMOVE)) {
    TranslateMessage(&msg);
    DispatchMessageW(&msg);
  }

  if (win->should_close) return 1;
  if (win->is_paused) return 2;
  return 0;
#else
  return 1;
#endif
}

// [[Rcpp::export]]
void close_native_window_cpp(SEXP win_ptr) {
#ifdef _WIN32
  if (TYPEOF(win_ptr) == EXTPTRSXP) {
    NativeWindow* win = (NativeWindow*)R_ExternalPtrAddr(win_ptr);
    if (win) {
      win->close();
    }
  }
#endif
}

// [[Rcpp::export]]
void set_native_window_title_cpp(SEXP win_ptr, std::string title) {
#ifdef _WIN32
  if (TYPEOF(win_ptr) == EXTPTRSXP) {
    NativeWindow* win = (NativeWindow*)R_ExternalPtrAddr(win_ptr);
    if (win && win->hwnd) {
      std::wstring wtitle(title.begin(), title.end());
      SetWindowTextW(win->hwnd, wtitle.c_str());
    }
  }
#endif
}

// [[Rcpp::export]]
bool has_native_video_cpp() {
#ifdef _WIN32
  ensure_mf_started();
  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
Rcpp::List get_native_video_info_cpp(std::string video_path) {
#ifdef _WIN32
  ensure_mf_started();

  int wlen = MultiByteToWideChar(CP_UTF8, 0, video_path.c_str(), -1, NULL, 0);
  if (wlen <= 0) {
    Rcpp::stop("Invalid video path encoding: %s", video_path.c_str());
  }
  std::wstring wpath(wlen, 0);
  MultiByteToWideChar(CP_UTF8, 0, video_path.c_str(), -1, &wpath[0], wlen);

  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 1);
  if (SUCCEEDED(hr) && pAttributes) {
    pAttributes->SetUINT32(MF_SOURCE_READER_ENABLE_VIDEO_PROCESSING, TRUE);
  }

  IMFSourceReader* pReader = nullptr;
  hr = MFCreateSourceReaderFromURL(wpath.c_str(), pAttributes, &pReader);
  if (pAttributes) pAttributes->Release();

  if (FAILED(hr) || !pReader) {
    Rcpp::stop("Cannot open video file: %s (HRESULT 0x%08X)", video_path.c_str(), hr);
  }

  PROPVARIANT var;
  PropVariantInit(&var);
  hr = pReader->GetPresentationAttribute((DWORD)MF_SOURCE_READER_MEDIASOURCE, MF_PD_DURATION, &var);
  double duration_sec = 0.0;
  if (SUCCEEDED(hr) && var.vt == VT_UI8) {
    duration_sec = (double)var.uhVal.QuadPart / 10000000.0;
  }
  PropVariantClear(&var);

  IMFMediaType* pNativeType = nullptr;
  hr = pReader->GetNativeMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, 0, &pNativeType);
  UINT32 width = 0, height = 0;
  UINT32 fps_num = 0, fps_den = 1;
  if (SUCCEEDED(hr) && pNativeType) {
    MFGetAttributeSize(pNativeType, MF_MT_FRAME_SIZE, &width, &height);
    MFGetAttributeRatio(pNativeType, MF_MT_FRAME_RATE, &fps_num, &fps_den);
    pNativeType->Release();
  }
  pReader->Release();

  double fps = (fps_den > 0 && fps_num > 0) ? ((double)fps_num / (double)fps_den) : 30.0;
  int est_frames = (int)std::round(duration_sec * fps);

  return Rcpp::List::create(
    Rcpp::Named("width") = (int)width,
    Rcpp::Named("height") = (int)height,
    Rcpp::Named("fps") = fps,
    Rcpp::Named("duration") = duration_sec,
    Rcpp::Named("estimated_frames") = est_frames
  );
#else
  Rcpp::stop("Native video reader is only supported on Windows.");
  return Rcpp::List::create();
#endif
}

static void render_c_progress_bar(int current, int total, double elapsed_sec) {
  const int bar_width = 30;
  int display_total = (total > current) ? total : current;
  float ratio = (display_total > 0) ? (float)current / (float)display_total : 0.0f;
  if (ratio > 1.0f) ratio = 1.0f;
  if (ratio < 0.0f) ratio = 0.0f;
  int filled = (int)(bar_width * ratio);

  char bar_str[64];
  for (int i = 0; i < bar_width; ++i) {
    if (i < filled) bar_str[i] = '=';
    else if (i == filled && filled < bar_width) bar_str[i] = '>';
    else bar_str[i] = ' ';
  }
  bar_str[bar_width] = '\0';

  int pct = (int)(ratio * 100.0f);
  int el_min = (int)elapsed_sec / 60;
  int el_sec = (int)elapsed_sec % 60;

  if (display_total > 0) {
    double eta = (current > 0 && ratio < 1.0f) ? (elapsed_sec / (double)ratio) * (1.0 - (double)ratio) : 0.0;
    int eta_min = (int)eta / 60;
    int eta_sec = (int)eta % 60;
    Rprintf("\rExtracting frames [%s] %3d%% (%d/%d) [%02d:%02d<--%02d:%02d]",
            bar_str, pct, current, display_total, el_min, el_sec, eta_min, eta_sec);
  } else {
    Rprintf("\rExtracting frames [%s] %d frames [%02d:%02d]",
            bar_str, current, el_min, el_sec);
  }
  R_FlushConsole();
}

// [[Rcpp::export]]
Rcpp::DataFrame split_video_native_cpp(std::string video_path,
                                       std::string output_dir,
                                       std::string prefix = "frame_",
                                       int digits = 6,
                                       std::string format = "jpg",
                                       int quality = 95,
                                       double target_fps = 0.0,
                                       int target_n = 0,
                                       bool verbose = true,
                                       Rcpp::Nullable<Rcpp::Function> callback = R_NilValue) {
#ifdef _WIN32
  ensure_mf_started();

  int wlen = MultiByteToWideChar(CP_UTF8, 0, video_path.c_str(), -1, NULL, 0);
  if (wlen <= 0) {
    Rcpp::stop("Invalid video path encoding: %s", video_path.c_str());
  }
  std::wstring wpath(wlen, 0);
  MultiByteToWideChar(CP_UTF8, 0, video_path.c_str(), -1, &wpath[0], wlen);

  IMFAttributes* pAttributes = nullptr;
  HRESULT hr = MFCreateAttributes(&pAttributes, 2);
  if (SUCCEEDED(hr) && pAttributes) {
    pAttributes->SetUINT32(MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS, TRUE);
    pAttributes->SetUINT32(MF_SOURCE_READER_ENABLE_VIDEO_PROCESSING, TRUE);
  }

  IMFSourceReader* pReader = nullptr;
  hr = MFCreateSourceReaderFromURL(wpath.c_str(), pAttributes, &pReader);
  if (pAttributes) pAttributes->Release();

  if (FAILED(hr) || !pReader) {
    Rcpp::stop("Cannot open video file: %s (HRESULT 0x%08X)", video_path.c_str(), hr);
  }

  struct SafeReader {
    IMFSourceReader* p = nullptr;
    ~SafeReader() { if (p) { p->Release(); p = nullptr; } }
  } reader_guard;
  reader_guard.p = pReader;

  pReader->SetStreamSelection((DWORD)MF_SOURCE_READER_ALL_STREAMS, FALSE);
  pReader->SetStreamSelection((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, TRUE);

  IMFMediaType* pMediaType = nullptr;
  hr = MFCreateMediaType(&pMediaType);
  if (SUCCEEDED(hr) && pMediaType) {
    pMediaType->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
    pMediaType->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_RGB32);
    hr = pReader->SetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, NULL, pMediaType);
    pMediaType->Release();
  }
  if (FAILED(hr)) {
    Rcpp::stop("Failed to configure native video decoder for RGB format (HRESULT 0x%08X)", hr);
  }

  IMFMediaType* pActualType = nullptr;
  hr = pReader->GetCurrentMediaType((DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM, &pActualType);
  UINT32 width = 0, height = 0;
  INT32 stride = 0;
  UINT32 fps_num = 0, fps_den = 1;
  if (SUCCEEDED(hr) && pActualType) {
    MFGetAttributeSize(pActualType, MF_MT_FRAME_SIZE, &width, &height);
    pActualType->GetUINT32(MF_MT_DEFAULT_STRIDE, (UINT32*)&stride);
    MFGetAttributeRatio(pActualType, MF_MT_FRAME_RATE, &fps_num, &fps_den);
    pActualType->Release();
  }
  if (stride == 0) stride = (INT32)(width * 4);
  int abs_stride = std::abs((int)stride);
  bool is_bottom_up = (stride < 0);
  double native_fps = (fps_den > 0 && fps_num > 0) ? ((double)fps_num / (double)fps_den) : 30.0;

  PROPVARIANT var;
  PropVariantInit(&var);
  hr = pReader->GetPresentationAttribute((DWORD)MF_SOURCE_READER_MEDIASOURCE, MF_PD_DURATION, &var);
  LONGLONG duration_100ns = 0;
  if (SUCCEEDED(hr) && var.vt == VT_UI8) {
    duration_100ns = var.uhVal.QuadPart;
  }
  PropVariantClear(&var);

  bool sample_by_n = (target_n > 0);
  bool sample_by_fps = (!sample_by_n && target_fps > 0.0);

  std::vector<LONGLONG> target_times;
  size_t target_idx = 0;
  LONGLONG fps_step_100ns = 0;
  LONGLONG next_fps_target = 0;

  if (sample_by_n) {
    target_times.resize(target_n);
    if (target_n == 1) {
      target_times[0] = 0;
    } else {
      LONGLONG frame_dur = (LONGLONG)(10000000.0 / native_fps);
      LONGLONG effective_duration = (duration_100ns > frame_dur) ? (duration_100ns - frame_dur) : duration_100ns;
      for (int i = 0; i < target_n; ++i) {
        target_times[i] = (LONGLONG)((double)i * (double)effective_duration / (double)(target_n - 1));
      }
    }
  } else if (sample_by_fps) {
    fps_step_100ns = (LONGLONG)(10000000.0 / target_fps);
    next_fps_target = 0;
  }

  std::vector<unsigned char> rgb_buffer((size_t)width * height * 3);
  std::vector<int> saved_frames;
  std::vector<std::string> saved_files;
  std::vector<double> saved_times;
  int saved_count = 0;
  bool is_png = (format == "png");
  std::string ext = is_png ? "png" : "jpg";
  if (digits < 1) digits = 6;
  if (quality < 1) quality = 1;
  if (quality > 100) quality = 100;

  int check_interrupt_counter = 0;
  LONGLONG last_sample_time = 0;
  bool last_frame_valid = false;

  auto t_start = std::chrono::steady_clock::now();
  int total_expected = 0;
  if (sample_by_n) {
    total_expected = target_n;
  } else if (sample_by_fps) {
    double dur_s = (double)duration_100ns / 10000000.0;
    total_expected = (int)std::max(1.0, std::round(dur_s * target_fps));
  } else {
    double dur_s = (double)duration_100ns / 10000000.0;
    total_expected = (int)std::max(1.0, std::round(dur_s * native_fps));
  }

  if (verbose) {
    render_c_progress_bar(0, total_expected, 0.0);
  }

  while (true) {
    if (++check_interrupt_counter % 20 == 0) {
      Rcpp::checkUserInterrupt();
    }

    if (sample_by_n && target_idx >= (size_t)target_n) {
      break;
    }

    DWORD streamIndex = 0, flags = 0;
    LONGLONG timestamp = 0;
    IMFSample* pSample = nullptr;

    HRESULT hr_rs = pReader->ReadSample(
      (DWORD)MF_SOURCE_READER_FIRST_VIDEO_STREAM,
      0,
      &streamIndex,
      &flags,
      &timestamp,
      &pSample
    );

    if (FAILED(hr_rs) || (flags & MF_SOURCE_READERF_ENDOFSTREAM)) {
      if (pSample) pSample->Release();
      break;
    }

    if (!pSample) {
      continue;
    }

    last_sample_time = timestamp;

    bool should_save = false;
    if (sample_by_n) {
      if (target_idx < (size_t)target_n && timestamp >= target_times[target_idx]) {
        should_save = true;
        target_idx++;
      }
    } else if (sample_by_fps) {
      if (timestamp >= next_fps_target) {
        should_save = true;
        while (next_fps_target <= timestamp) {
          next_fps_target += fps_step_100ns;
        }
      }
    } else {
      should_save = true;
    }

    // CRITICAL OPTIMIZATION: Do not convert or lock large pixel buffers for discarded frames!
    if (!should_save) {
      pSample->Release();
      continue;
    }

    IMFMediaBuffer* pBuffer = nullptr;
    if (SUCCEEDED(pSample->ConvertToContiguousBuffer(&pBuffer)) && pBuffer) {
      BYTE* pData = nullptr;
      DWORD maxLen = 0, curLen = 0;
      if (SUCCEEDED(pBuffer->Lock(&pData, &maxLen, &curLen)) && pData) {
        if (curLen >= (DWORD)(abs_stride * height)) {
          #pragma omp parallel for schedule(static) if(height > 100)
          for (int y = 0; y < (int)height; ++y) {
            int src_y = is_bottom_up ? ((int)height - 1 - y) : y;
            const unsigned char* src_row = pData + (size_t)src_y * abs_stride;
            unsigned char* dst_row = rgb_buffer.data() + (size_t)y * width * 3;
            for (int x = 0; x < (int)width; ++x) {
              dst_row[x * 3 + 0] = src_row[x * 4 + 2]; // R
              dst_row[x * 3 + 1] = src_row[x * 4 + 1]; // G
              dst_row[x * 3 + 2] = src_row[x * 4 + 0]; // B
            }
          }

          char fname[256];
          snprintf(fname, sizeof(fname), "%s%0*d.%s", prefix.c_str(), digits, saved_count + 1, ext.c_str());
          std::string out_filepath = output_dir + "/" + std::string(fname);

          int write_res = 0;
          if (is_png) {
            write_res = stbi_write_png(out_filepath.c_str(), (int)width, (int)height, 3, rgb_buffer.data(), (int)width * 3);
          } else {
            write_res = stbi_write_jpg(out_filepath.c_str(), (int)width, (int)height, 3, rgb_buffer.data(), quality);
          }

          if (write_res) {
            saved_count++;
            saved_frames.push_back(saved_count);
            saved_files.push_back(out_filepath);
            saved_times.push_back((double)timestamp / 10000000.0);
            if (verbose) {
              auto t_now = std::chrono::steady_clock::now();
              double el = std::chrono::duration<double>(t_now - t_start).count();
              render_c_progress_bar(saved_count, total_expected, el);
            }
            if (callback.isNotNull()) {
              Rcpp::Function cb(callback.get());
              cb();
            }
          }
          last_frame_valid = true;
        }
        pBuffer->Unlock();
      }
      pBuffer->Release();
    }
    pSample->Release();
  }

  // If sample_by_n and only target_n - 1 frames were saved due to stream end timing,
  // save the last valid frame so the user gets exactly target_n frames.
  if (sample_by_n && saved_count == (target_n - 1) && last_frame_valid) {
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%0*d.%s", prefix.c_str(), digits, saved_count + 1, ext.c_str());
    std::string out_filepath = output_dir + "/" + std::string(fname);
    int write_res = 0;
    if (is_png) {
      write_res = stbi_write_png(out_filepath.c_str(), (int)width, (int)height, 3, rgb_buffer.data(), (int)width * 3);
    } else {
      write_res = stbi_write_jpg(out_filepath.c_str(), (int)width, (int)height, 3, rgb_buffer.data(), quality);
    }
    if (write_res) {
      saved_count++;
      saved_frames.push_back(saved_count);
      saved_files.push_back(out_filepath);
      saved_times.push_back((double)last_sample_time / 10000000.0);
      if (verbose) {
        auto t_now = std::chrono::steady_clock::now();
        double el = std::chrono::duration<double>(t_now - t_start).count();
        render_c_progress_bar(saved_count, total_expected, el);
      }
      if (callback.isNotNull()) {
        Rcpp::Function cb(callback.get());
        cb();
      }
    }
  }

  if (verbose) {
    auto t_now = std::chrono::steady_clock::now();
    double el = std::chrono::duration<double>(t_now - t_start).count();
    render_c_progress_bar(saved_count, (total_expected > saved_count ? saved_count : total_expected), el);
    Rprintf("\n");
    R_FlushConsole();
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("frame") = saved_frames,
    Rcpp::Named("file") = saved_files,
    Rcpp::Named("time") = saved_times,
    Rcpp::Named("width") = Rcpp::IntegerVector(saved_frames.size(), (int)width),
    Rcpp::Named("height") = Rcpp::IntegerVector(saved_frames.size(), (int)height),
    Rcpp::Named("stringsAsFactors") = false
  );
#else
  Rcpp::stop("Native video reader is only supported on Windows.");
  return Rcpp::DataFrame::create();
#endif
}
