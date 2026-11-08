import time
import torch
import numpy as np
from ultralytics import YOLO
import onnxruntime as ort

pt_path = r"D:/Downloads/pliman4/wood_log_model/weights/best.pt"
onnx_path = r"D:/Downloads/pliman4/wood_log_model/weights/best.onnx"

print("=" * 65)
print("BENCHMARK PYTHON ULTRALYTICS & ONNX RUNTIME (PURAMENTE INFERÊNCIA)")
print("=" * 65)

# 1. Ultralytics PyTorch (.pt) - GPU CUDA
model_pt = YOLO(pt_path)
model_pt.to("cuda")
dummy_np = np.zeros((640, 640, 3), dtype=np.uint8)

# Warmup
for _ in range(15):
    _ = model_pt.predict(dummy_np, device=0, verbose=False)
torch.cuda.synchronize()

n_iter = 100
t0 = time.perf_counter()
for _ in range(n_iter):
    _ = model_pt.predict(dummy_np, device=0, verbose=False)
torch.cuda.synchronize()
t_pt_predict = (time.perf_counter() - t0) / n_iter * 1000.0

# Pure model forward pass (sem preprocess/NMS do ultralytics wrapper)
dummy_tensor = torch.zeros((1, 3, 640, 640), dtype=torch.float32, device="cuda")
for _ in range(15):
    _ = model_pt.model(dummy_tensor)
torch.cuda.synchronize()

t0 = time.perf_counter()
for _ in range(n_iter):
    _ = model_pt.model(dummy_tensor)
torch.cuda.synchronize()
t_pt_forward = (time.perf_counter() - t0) / n_iter * 1000.0

print(f"1) Python Ultralytics (best.pt | GPU CUDA):")
print(f"   - model.predict() completo : {t_pt_predict:.2f} ms/frame  ({1000.0/t_pt_predict:.1f} FPS)")
print(f"   - pure forward pass model() : {t_pt_forward:.2f} ms/frame  ({1000.0/t_pt_forward:.1f} FPS)")

# 2. Python ONNX Runtime (best.onnx | GPU CUDA)
try:
    providers = ['CUDAExecutionProvider', 'CPUExecutionProvider']
    session_cuda = ort.InferenceSession(onnx_path, providers=providers)
    input_name = session_cuda.get_inputs()[0].name
    input_feed = {input_name: np.zeros((1, 3, 640, 640), dtype=np.float32)}

    for _ in range(15):
        _ = session_cuda.run(None, input_feed)

    t0 = time.perf_counter()
    for _ in range(n_iter):
        _ = session_cuda.run(None, input_feed)
    t_onnx_cuda = (time.perf_counter() - t0) / n_iter * 1000.0

    print(f"\n2) Python ONNXRuntime (best.onnx | CUDAExecutionProvider):")
    print(f"   - session.run() puro        : {t_onnx_cuda:.2f} ms/frame  ({1000.0/t_onnx_cuda:.1f} FPS)")
except Exception as e:
    print(f"\n2) Python ONNXRuntime CUDA falhou: {e}")

# 3. Python ONNX Runtime (best.onnx | CPU)
session_cpu = ort.InferenceSession(onnx_path, providers=['CPUExecutionProvider'])
input_name = session_cpu.get_inputs()[0].name
input_feed = {input_name: np.zeros((1, 3, 640, 640), dtype=np.float32)}
for _ in range(5):
    _ = session_cpu.run(None, input_feed)

t0 = time.perf_counter()
for _ in range(30):
    _ = session_cpu.run(None, input_feed)
t_onnx_cpu = (time.perf_counter() - t0) / 30 * 1000.0
print(f"\n3) Python ONNXRuntime (best.onnx | CPUExecutionProvider):")
print(f"   - session.run() CPU         : {t_onnx_cpu:.2f} ms/frame  ({1000.0/t_onnx_cpu:.1f} FPS)")
print("=" * 65)
