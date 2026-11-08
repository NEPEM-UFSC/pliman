
import time
import sys
from ultralytics import YOLO

model_path = r"D:\Downloads\pliman4\wood_log_model\weights\best.pt"
video_path = r"D:\Downloads\pliman4\toras2.mp4"

print("Loading model...", flush=True)
model = YOLO(model_path)

print("Starting inference on 100 frames...", flush=True)
t0 = time.perf_counter()
count = 0
for r in model.predict(source=video_path, conf=0.5, device=0, stream=True, verbose=False):
    count += 1
    if count % 20 == 0:
        print(f"Frame {count}: {count / (time.perf_counter() - t0):.1f} FPS", flush=True)
    if count >= 100:
        break

t1 = time.perf_counter()
elapsed = t1 - t0
fps = count / elapsed
print(f"Total: {count} frames in {elapsed:.3f} s -> {fps:.1f} FPS ({1000/fps:.1f} ms/frame)", flush=True)
