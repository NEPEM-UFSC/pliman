# Licenses of the pre-trained models

The deep learning functions of **pliman** (`image_detect_dl()`, `image_segment_dl()`,
`image_remove_bg_dl()`, `image_depth_dl()`, `video_detect_dl()`, ...) use pre-trained
models that are **not part of the package**. They are downloaded on first use, and each
one keeps the license chosen by its authors. The pliman source code is distributed under
the GNU General Public License (>= 3); this does not change the license of any model.

> **Non-commercial use.** Several of the models listed below forbid commercial use or
> are released under the AGPL-3.0, and some licenses could not be verified. For this
> reason, the models available through pliman are provided for **research, educational
> and non-commercial use only**. If you plan to use them in a commercial product or
> service, check the license of every model you use and obtain the necessary permissions.

This table is informational and is **not legal advice**. Licenses change over time:
always read the license distributed with the weights you download.

## Restricted: non-commercial use or copyleft

| Model name(s) in pliman | License | Commercial use |
|:---|:---|:---|
| `rmbg-2.0` (BRIA RMBG 2.0) | CC BY-NC 4.0 | No (commercial agreement with BRIA required) |
| `rmbg-1.4` (BRIA RMBG 1.4) | BRIA custom license (`bria-rmbg-1.4`) | No (non-commercial) |
| `yolo26*` (detection, segmentation, pose, classification, OBB) | AGPL-3.0 (Ultralytics) | Only with an open-source release of the whole project or an Ultralytics Enterprise License |
| `yolov8s/m/l/x-worldv2` (YOLO-World) | AGPL-3.0 (Ultralytics weights); the original repository is GPL-3.0 | Copyleft terms apply |
| `yolo_nas_s`, `yolo_nas_l` (YOLO-NAS) | Code Apache-2.0, **weights under the Deci license** | No (research and non-commercial) |

## Permissive licenses (verified)

| Model name(s) in pliman | License | Commercial use |
|:---|:---|:---|
| `depth-anything-v2` (Depth Anything V2 **Small**) | Apache-2.0 | Yes. The Base, Large and Giant sizes are CC BY-NC 4.0 |
| `ben2` (BEN2 Base) | MIT | Yes. The full BEN2 model of PramaLLC is a separate commercial product |
| `grounded-sam` (Grounding DINO Tiny) | Apache-2.0 | Yes |
| `ppocr_*`, `ppocrv5_*` (PaddleOCR PP-OCRv4 / v5) | Apache-2.0 | Yes |

## Not verified (license expected from the upstream project, please confirm)

| Model name(s) in pliman | Expected license |
|:---|:---|
| `u2net`, `u2netp`, `silueta` (U-2-Net and derivatives) | Apache-2.0 |
| `isnet-general-use` (IS-Net) | Apache-2.0 |
| `birefnet-lite` (BiRefNet) | MIT |
| `sam2.1`, `persam` (SAM 2.1) | Apache-2.0 |
| `dinov2` (DINOv2 ViT-S/14) | Apache-2.0 |
| `realesrgan`, `realesrgan_x2_full`, `realesrgan_x4_full` (Real-ESRGAN) | BSD-3-Clause |
| `clip-vit-b32` (CLIP) | MIT |
| `craft` (CRAFT) | MIT |
| `trocr_*` (TrOCR) | MIT |

## Unknown (license not found)

| Model name(s) in pliman | Note |
|:---|:---|
| `withoutbg` | Combines Depth Anything V2 features with the withoutBG open weights, whose license was not checked |
| `ifcnn` | The upstream repository does not state a license |
| `stardist` (DSB 2018) | Not verified (the StarDist code is usually BSD-3-Clause) |
| `SAFMN2x`, `SAFMN4x`, `SAFMN2x_sharper` (PureScale) | Not verified: community-trained weights on top of SAFMN |

## Using your own models

Functions that take a `model` argument also accept the path to your own `.onnx` file.
In that case, the license of that file is your responsibility.
