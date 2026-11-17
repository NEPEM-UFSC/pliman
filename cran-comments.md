## Test environments
* local Windows 11 x64, R 4.5.1
* win-builder (devel and release)
* R-hub builder

## R CMD check results
0 errors | 0 warnings | 0 notes

## Release summary
This is a major release (`4.0.0`) establishing a comprehensive deep learning and computer vision framework natively in R with zero mandatory Python dependencies.

Key highlights of this release:
* **Deep Learning Foundation Architectures**:
  - Segment Anything 2.1 (SAM 2.1) and PerSAM (Personalized SAM) for one-shot visual exemplar segmentation and density counting.
  - Grounded-SAM (Grounding DINO + SAM 2.1) for zero-shot text-prompted instance segmentation.
  - Optical Character Recognition (OCR) via DBNet, PaddleOCRv4/v5, and TrOCR with automated batch regex image renaming and rollback.
  - Zero-DCE and NAFNet neural image restoration for low-light enhancement, denoising, and deblurring.
  - Deep semantic feature embeddings via OpenAI CLIP (512-D) and Meta DINOv2 (384-D) with spatial PCA visualization.
  - Real-ESRGAN neural super-resolution for image upscaling and similarity map refinement.
  - YOLO26 and StarDist real-time instance segmentation and star-convex polygon detection.
  - Keypoint pose estimation and joint skeleton modeling.

* **End-to-End YOLO Dataset & Training Ecosystem**:
  - Interactive manual bounding box labeling (`yolo_dataset_export(model = "manual")`) with lightweight 2-click box definition, live canvas updates, and zero DL overhead.
  - Semi-automated active learning annotation workflows (PerSAM exemplars, Grounded-SAM prompts, YOLO pseudo-labeling, few-shot model fitting).
  - Full YOLO dataset management suite: auditing, deduplication, merging, splitting, offline augmentation, copy-paste foreground augmentation, high-resolution tiling, balancing, and anonymization.
  - In-R YOLO model training, checkpoint recovery, evaluation, error analysis, and INT8/FP16/TensorRT edge quantization.

* **3D Reconstruction, Video & High-Throughput Phenotyping**:
  - 3D Gaussian Splatting (`gs_*`) for plant volumetric modeling and rendering.
  - Video multi-object tracking (`video_detect_dl()`) with native C++ state handling.
  - Panoramic image stitching and focus stacking.

* **Dependencies & Architecture**:
  - Standalone C++ algorithms and ONNX Runtime C API integration eliminating heavy external dependencies.
  - Scalable background execution via `mirai`.
  - Comprehensive `cli` console reporting for batch operations.

## Note on Bioconductor / Optional Dependencies
* {pliman} optionally suggests {EBImage} for specific legacy workflows. All required biological and image processing algorithms are implemented natively in C++. If {EBImage} is absent, clear internal instructions guide the user on optional installation.
