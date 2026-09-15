# utils_neural.R — Pre-trained Deep Learning Segmentation and Background Removal for pliman
#
# ==============================================================================
# DEEP LEARNING ARCHITECTURE & ONNX RUNTIME DOCUMENTATION
# ==============================================================================
#' @title Deep Learning Background Removal & Salient Object Segmentation in pliman
#' @name utils_neural
#' @description
#' The `pliman` neural module provides high-accuracy, zero-Python foreground segmentation
#' and background removal powered by pre-trained Deep Learning models in **ONNX**
#' (Open Neural Network Exchange) format.
#'
#' @details
#' # 1. Supported Neural Architectures
#' * **`"u2netp"` (Default for CPU)**:
#'   * **Size:** ~4.4 MB | **Input Resolution:** 320 × 320
#'   * **Architecture:** Ultra-lightweight U2-Net with nested U-structure residual blocks.
#'   * **Best For:** Fast batch processing on standard laptops and CPUs without dedicated GPU.
#'
#' * **`"ben2"` (Boundary-Aware Extraction Network)**:
#'   * **Size:** ~212.6 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** BEN2 Base (Prama LLC) with boundary-focused refinement stream.
#'   * **Best For:** Crisp sub-pixel object boundaries, fine petioles, leaf margins, roots, and hairs.
#'
#' * **`"grounded-sam"` (Zero-Shot Text-Prompted Detection & Instance Segmentation)**:
#'   * **Size:** ~314.8 MB (Grounding DINO ~194.4 MB + SAM 2.1 ~120.2 MB) | **Input Resolution:** 800 × 800 (DINO) & 1024 × 1024 (SAM)
#'   * **Architecture:** Combines Grounding DINO Tiny for open-vocabulary text-prompted bounding box detection with SAM 2.1 Hiera-Tiny for instance mask extraction.
#'   * **Best For:** Detecting and segmenting specific objects using natural language prompts (e.g. `"cow"`, `"leaf"`, `"lesion"`), with support for bounding-box-only detection (`mask = FALSE, bbox = TRUE`).
#'
#' * **`"sam2.1"` (Segment Anything 2.1 by Meta AI, alias `"sam2"`)**:
#'   * **Size:** ~120.2 MB (Encoder + Decoder) | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** Meta AI Hiera-Tiny foundation vision transformer with promptable mask decoder.
#'   * **Best For:** Zero-shot foundation segmentation with point prompts (`pick_object = TRUE`), box prompts, or center point defaults.
#'
#' * **`"sam3.1"` (Segment Anything 3.1 by Meta AI, alias `"sam3"`)**:
#'   * **Size:** ~868.1 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** Meta AI Segment Anything 3.1 with real-time concept-driven segmentation.
#'   * **Best For:** Advanced foundation segmentation across diverse biological structures.
#'
#' * **`"birefnet-lite"` (Bilateral Reference Network Lite)**:
#'   * **Size:** ~213.6 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** Bilateral Reference Network for High-Resolution Dichotomous Image Segmentation (BiRefNet Lite).
#'   * **Best For:** Extremely fine margins, hair-thin leaf serrations, complex petioles, translucent halos, and subtle lesions.
#'
#' * **`"isnet-general-use"` (High Precision Boundary Matting)**:
#'   * **Size:** ~170.4 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** IS-Net (Intermediate Supervision Network) specialized in high-resolution boundary matting.
#'   * **Best For:** Plant leaves, fine serrated margins, lesion borders, chlorotic halos, veins, and complex natural backgrounds.
#'
#' * **`"rmbg-1.4"` (BRIA AI Standard, alias `"rmbg"`)**:
#'   * **Size:** ~168.0 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** BRIA RMBG 1.4 trained on extensive curated foreground datasets.
#'   * **Best For:** General background removal under challenging lighting, shadows, reflections, and fine object edges.
#'
#' * **`"rmbg-2.0"` (BRIA AI Next-Gen Model)**:
#'   * **Size:** ~976.9 MB | **Input Resolution:** 1024 × 1024
#'   * **Architecture:** BRIA RMBG 2.0 based on a BiRefNet backbone trained on extensive curated commercial datasets.
#'   * **Best For:** Highest-fidelity boundary matting, transparent surfaces, fine plant hairs, and complex scenes.
#'
#' * **`"withoutbg"` (State-of-the-Art DepthAnythingV2 + ConvNeXt Matting)**:
#'   * **Size:** ~433.4 MB | **Input Resolution:** 448 × 448
#'   * **Architecture:** withoutBG Open Weights: DINOv3 ConvNeXt-fused U-Net guided by DepthAnythingV2 depth features.
#'   * **Best For:** Complex 3D depth-separated backgrounds, transparent materials, and fine details.
#'
#' * **`"silueta"`**:
#'   * **Size:** ~42.1 MB | **Input Resolution:** 320 × 320
#'   * **Architecture:** Compact silhouette extraction network.
#'   * **Best For:** Balanced speed-accuracy tradeoff on edge/CPU hardware.
#'
#' * **`"u2net"`**:
#'   * **Size:** ~167.8 MB | **Input Resolution:** 320 × 320
#'   * **Architecture:** Full-depth U2-Net salient object detection model.
#'   * **Best For:** Robust multi-scale salient object detection on standard resolution images.
#'
#' # 2. Zero-Python Architecture (Native C++ Runtime)
#'
#' Unlike standard Python-dependent deep learning packages, `pliman` executes ONNX models
#' natively in C++ via direct Microsoft ONNX Runtime C API dynamic bindings.
#'
#' * **Zero Python Requirement:** No Python, Conda, Virtualenv, or pip required.
#' * **Direct C++ Binding:** Connects directly to Microsoft's official ONNX Runtime C++ engine (`onnxruntime.dll` on Windows / `libonnxruntime.so` on Linux / `libonnxruntime.dylib` on macOS).
#' * **CPU & GPU Support:** Runs optimized CPU multi-threading (AVX2/AVX512) or CUDA GPU acceleration.
#'
#' # 3. Model Weight Storage & Cache
#'
#' Downloaded `.onnx` model weights are stored permanently in the user's R data directory:
#' `tools::R_user_dir("pliman", which = "data")/models`
#' (e.g., `C:/Users/<User>/AppData/Local/R/data/R/pliman/models` on Windows).
#' This follows CRAN policy and ensures models are preserved across R session restarts
#' and package upgrades.
#'
#' # 4. Quick Setup Guide
#'
#' To configure everything in a single call:
#' ```r
#' library(pliman)
#'
#' # 1. Set up the environment (downloads ONNX Runtime C++ library & default models)
#' pliman_configure_dl()
#'
#' # 2. Extract a binary mask
#' mask <- image_binary_dl(img, model = "u2netp")
#'
#' # 3. Remove background with transparent alpha channel (RGBA)
#' img_trans <- image_remove_bg_dl(img, model = "isnet-general-use")
#'
#' # 4. General segmentation
#' seg <- image_segment_dl(img, model = "u2netp")
#' ```
NULL

# ==============================================================================
# SECTION 1: MODEL MANAGEMENT & DIRECTORY MANAGEMENT
# ==============================================================================

#' Model Directory Management for pliman Deep Learning Models
#'
#' Gets or sets the local directory where pre-trained ONNX models are stored.
#' Defaults to the standard R user data directory:
#' `tools::R_user_dir("pliman", which = "data")/models`.
#'
#' @param dir Optional character string specifying a custom directory path.
#'   If provided, creates the directory if it does not exist.
#' @return A character string with the normalized absolute path to the model directory.
#' @export
#' @examples
#' pliman_model_dir()
pliman_model_dir <- function(dir = NULL) {
  if (!is.null(dir)) {
    dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(dir)
  }

  default_dir <- file.path(tools::R_user_dir("pliman", which = "data"), "models")
  default_dir <- normalizePath(default_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(default_dir)) {
    dir.create(default_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(default_dir)
}

# Internal helper to resolve user-friendly model aliases
.resolve_model_name <- function(model) {
  aliases <- c(
    "rmbg"                   = "rmbg-1.4",
    "bria"                   = "rmbg-1.4",
    "bria-rmbg"              = "rmbg-1.4",
    "bria_rmbg"              = "rmbg-1.4",
    "rmbg1.4"                = "rmbg-1.4",
    "rmbg_1.4"               = "rmbg-1.4",
    "rmbg2"                  = "rmbg-2.0",
    "rmbg-2"                 = "rmbg-2.0",
    "rmbg_2.0"               = "rmbg-2.0",
    "rmbg2.0"                = "rmbg-2.0",
    "bria-rmbg-2.0"          = "rmbg-2.0",
    "birefnet"               = "birefnet-lite",
    "birefnet_lite"          = "birefnet-lite",
    "without_bg"             = "withoutbg",
    "without-bg"             = "withoutbg",
    "withoutbg-open-weights" = "withoutbg",
    "u2net_p"                = "u2netp",
    "ben"                    = "ben2",
    "ben2"                   = "ben2",
    "ben-2"                  = "ben2",
    "boundary-aware"         = "ben2",
    "boundary_aware"         = "ben2",
    "boundary-aware-extraction-network" = "ben2",
    "sam"                    = "sam2.1",
    "sam2"                   = "sam2.1",
    "sam-2"                  = "sam2.1",
    "sam2.1"                 = "sam2.1",
    "sam-2.1"                = "sam2.1",
    "sam3"                   = "sam3.1",
    "sam-3"                  = "sam3.1",
    "sam3.1"                 = "sam3.1",
    "sam-3.1"                = "sam3.1",
    "segment-anything"       = "sam2.1",
    "segment_anything"       = "sam2.1",
    "segmentanything"        = "sam2.1",
    "grounded-sam"           = "grounded-sam",
    "grounded_sam"           = "grounded-sam",
    "groundedsam"            = "grounded-sam",
    "grounding-dino"         = "grounded-sam",
    "grounding_dino"         = "grounded-sam",
    "groundingdino"          = "grounded-sam",
    "persam"                 = "persam",
    "per-sam"                = "persam",
    "persam2"                = "persam",
    "persam2.1"              = "persam",
    "depth-anything"         = "depth-anything-v2",
    "depthanything"          = "depth-anything-v2",
    "depth_anything"         = "depth-anything-v2",
    "depth-anything-v2"      = "depth-anything-v2",
    "depthanythingv2"        = "depth-anything-v2",
    "depth_anything_v2"      = "depth-anything-v2",
    "depth"                  = "depth-anything-v2",
    "dinov2"                 = "dinov2",
    "dino-v2"                = "dinov2",
    "dino_v2"                = "dinov2",
    "dino"                   = "dinov2",
    # YOLO Object Detection
    "yolo"                   = "yolo26n",
    "yolo26"                 = "yolo26n",
    "yolo-26"                = "yolo26n",
    "yolo26n"                = "yolo26n",
    "yolo-26n"               = "yolo26n",
    "yolo26s"                = "yolo26s",
    "yolo-26s"               = "yolo26s",
    "yolo26m"                = "yolo26m",
    "yolo-26m"               = "yolo26m",
    "yolo26l"                = "yolo26l",
    "yolo-26l"               = "yolo26l",
    "yolo26x"                = "yolo26x",
    "yolo-26x"               = "yolo26x",
    "yolo11"                 = "yolo26n",
    "yolo-11"                = "yolo26n",
    "yolo11n"                = "yolo26n",
    "yolo-11n"               = "yolo26n",
    "yolov8"                 = "yolo26n",
    "yolo8"                  = "yolo26n",

    # YOLO Instance Segmentation
    "yolo-seg"               = "yolo26n-seg",
    "yoloseg"                = "yolo26n-seg",
    "yolo26-seg"             = "yolo26n-seg",
    "yolo-26-seg"            = "yolo26n-seg",
    "yolo26n-seg"            = "yolo26n-seg",
    "yolo26s-seg"            = "yolo26s-seg",
    "yolo26m-seg"            = "yolo26m-seg",
    "yolo26l-seg"            = "yolo26l-seg",
    "yolo26x-seg"            = "yolo26x-seg",
    "yolo11-seg"             = "yolo26n-seg",
    "yolo11n-seg"            = "yolo26n-seg",

    # YOLO Pose Estimation
    "yolo-pose"              = "yolo26n-pose",
    "yolopose"               = "yolo26n-pose",
    "yolo26-pose"            = "yolo26n-pose",
    "yolo-26-pose"           = "yolo26n-pose",
    "yolo26n-pose"           = "yolo26n-pose",
    "yolo26s-pose"           = "yolo26s-pose",
    "yolo26m-pose"           = "yolo26m-pose",
    "yolo26l-pose"           = "yolo26l-pose",
    "yolo26x-pose"           = "yolo26x-pose",

    # YOLO Classification
    "yolo-cls"               = "yolo26n-cls",
    "yolocls"                = "yolo26n-cls",
    "yolo26-cls"             = "yolo26n-cls",
    "yolo-26-cls"            = "yolo26n-cls",
    "yolo26n-cls"            = "yolo26n-cls",
    "yolo26s-cls"            = "yolo26s-cls",
    "yolo26m-cls"            = "yolo26m-cls",
    "yolo26l-cls"            = "yolo26l-cls",
    "stardist"               = "stardist",
    "star-dist"              = "stardist",
    "stardist-dsb2018"       = "stardist",
    "star_dist"              = "stardist",
    "realesrgan"             = "realesrgan-compact",
    "real-esrgan"            = "realesrgan-compact",
    "realesrgan-compact"     = "realesrgan-compact",
    "real_esrgan"            = "realesrgan-compact",
    "superres"               = "realesrgan-compact",
    "super-res"              = "realesrgan-compact"
  )
  vapply(model, function(m) {
    ml <- tolower(m)
    if (ml %in% names(aliases)) aliases[[ml]] else m
  }, character(1), USE.NAMES = FALSE)
}

.coco_classes <- c(
  "person", "bicycle", "car", "motorcycle", "airplane", "bus", "train", "truck", "boat", "traffic light",
  "fire hydrant", "stop sign", "parking meter", "bench", "bird", "cat", "dog", "horse", "sheep", "cow",
  "elephant", "bear", "zebra", "giraffe", "backpack", "umbrella", "handbag", "tie", "suitcase", "frisbee",
  "skis", "snowboard", "sports ball", "kite", "baseball bat", "baseball glove", "skateboard", "surfboard", "tennis racket", "bottle",
  "wine glass", "cup", "fork", "knife", "spoon", "bowl", "banana", "apple", "sandwich", "orange",
  "broccoli", "carrot", "hot dog", "pizza", "donut", "cake", "chair", "couch", "potted plant", "bed",
  "dining table", "toilet", "tv", "laptop", "mouse", "remote", "keyboard", "cell phone", "microwave", "oven",
  "toaster", "sink", "refrigerator", "book", "clock", "vase", "scissors", "teddy bear", "hair drier", "toothbrush"
)

.coco_keypoints <- c(
  "nose", "left_eye", "right_eye", "left_ear", "right_ear",
  "left_shoulder", "right_shoulder", "left_elbow", "right_elbow",
  "left_wrist", "right_wrist", "left_hip", "right_hip",
  "left_knee", "right_knee", "left_ankle", "right_ankle"
)

# Standard COCO 17-keypoint anatomical skeleton pairs (1-indexed)
.coco_skeleton_pairs <- matrix(c(
  1, 2,  1, 3,  2, 4,  3, 5,       # facial (nose-eye-ear)
  6, 7,  6, 12, 7, 13, 12, 13,     # torso (shoulders & hips)
  6, 8,  8, 10,                    # left arm (shoulder-elbow-wrist)
  7, 9,  9, 11,                    # right arm (shoulder-elbow-wrist)
  12, 14, 14, 16,                  # left leg (hip-knee-ankle)
  13, 15, 15, 17                   # right leg (hip-knee-ankle)
), ncol = 2, byrow = TRUE)

.imagenet_classes <- c(
  'tench', 'goldfish', 'great_white_shark', 'tiger_shark', 'hammerhead', 'electric_ray', 'stingray', 'cock', 'hen', 'ostrich',
  'brambling', 'goldfinch', 'house_finch', 'junco', 'indigo_bunting', 'robin', 'bulbul', 'jay', 'magpie', 'chickadee',
  'water_ouzel', 'kite', 'bald_eagle', 'vulture', 'great_grey_owl', 'European_fire_salamander', 'common_newt', 'eft', 'spotted_salamander', 'axolotl',
  'bullfrog', 'tree_frog', 'tailed_frog', 'loggerhead', 'leatherback_turtle', 'mud_turtle', 'terrapin', 'box_turtle', 'banded_gecko', 'common_iguana',
  'American_chameleon', 'whiptail', 'agama', 'frilled_lizard', 'alligator_lizard', 'Gila_monster', 'green_lizard', 'African_chameleon', 'Komodo_dragon', 'African_crocodile',
  'American_alligator', 'triceratops', 'thunder_snake', 'ringneck_snake', 'hognose_snake', 'green_snake', 'king_snake', 'garter_snake', 'water_snake', 'vine_snake',
  'night_snake', 'boa_constrictor', 'rock_python', 'Indian_cobra', 'green_mamba', 'sea_snake', 'horned_viper', 'diamondback', 'sidewinder', 'trilobite',
  'harvestman', 'scorpion', 'black_and_gold_garden_spider', 'barn_spider', 'garden_spider', 'black_widow', 'tarantula', 'wolf_spider', 'tick', 'centipede',
  'black_grouse', 'ptarmigan', 'ruffed_grouse', 'prairie_chicken', 'peacock', 'quail', 'partridge', 'African_grey', 'macaw', 'sulphur-crested_cockatoo',
  'lorikeet', 'coucal', 'bee_eater', 'hornbill', 'hummingbird', 'jacamar', 'toucan', 'drake', 'red-breasted_merganser', 'goose',
  'black_swan', 'tusker', 'echidna', 'platypus', 'wallaby', 'koala', 'wombat', 'jellyfish', 'sea_anemone', 'brain_coral',
  'flatworm', 'nematode', 'conch', 'snail', 'slug', 'sea_slug', 'chiton', 'chambered_nautilus', 'Dungeness_crab', 'rock_crab',
  'fiddler_crab', 'king_crab', 'American_lobster', 'spiny_lobster', 'crayfish', 'hermit_crab', 'isopod', 'white_stork', 'black_stork', 'spoonbill',
  'flamingo', 'little_blue_heron', 'American_egret', 'bittern', 'crane_(bird)', 'limpkin', 'European_gallinule', 'American_coot', 'bustard', 'ruddy_turnstone',
  'red-backed_sandpiper', 'redshank', 'dowitcher', 'oystercatcher', 'pelican', 'king_penguin', 'albatross', 'grey_whale', 'killer_whale', 'dugong',
  'sea_lion', 'Chihuahua', 'Japanese_spaniel', 'Maltese_dog', 'Pekinese', 'Shih-Tzu', 'Blenheim_spaniel', 'papillon', 'toy_terrier', 'Rhodesian_ridgeback',
  'Afghan_hound', 'basset', 'beagle', 'bloodhound', 'bluetick', 'black-and-tan_coonhound', 'Walker_hound', 'English_foxhound', 'redbone', 'borzoi',
  'Irish_wolfhound', 'Italian_greyhound', 'whippet', 'Ibizan_hound', 'Norwegian_elkhound', 'otterhound', 'Saluki', 'Scottish_deerhound', 'Weimaraner', 'Staffordshire_bullterrier',
  'American_Staffordshire_terrier', 'Bedlington_terrier', 'Border_terrier', 'Kerry_blue_terrier', 'Irish_terrier', 'Norfolk_terrier', 'Norwich_terrier', 'Yorkshire_terrier', 'wire-haired_fox_terrier', 'Lakeland_terrier',
  'Sealyham_terrier', 'Airedale', 'cairn', 'Australian_terrier', 'Dandie_Dinmont', 'Boston_bull', 'miniature_schnauzer', 'giant_schnauzer', 'standard_schnauzer', 'Scotch_terrier',
  'Tibetan_terrier', 'silky_terrier', 'soft-coated_wheaten_terrier', 'West_Highland_white_terrier', 'Lhasa', 'flat-coated_retriever', 'curly-coated_retriever', 'golden_retriever', 'Labrador_retriever', 'Chesapeake_Bay_retriever',
  'German_short-haired_pointer', 'vizsla', 'English_setter', 'Irish_setter', 'Gordon_setter', 'Brittany_spaniel', 'clumber', 'English_springer', 'Welsh_springer_spaniel', 'cocker_spaniel',
  'Sussex_spaniel', 'Irish_water_spaniel', 'kuvasz', 'schipperke', 'groenendael', 'malinois', 'briard', 'kelpie', 'komondor', 'Old_English_sheepdog',
  'Shetland_sheepdog', 'collie', 'Border_collie', 'Bouvier_des_Flandres', 'Rottweiler', 'German_shepherd', 'Doberman', 'miniature_pinscher', 'Greater_Swiss_Mountain_dog', 'Bernese_mountain_dog',
  'Appenzeller', 'EntleBucher', 'boxer', 'bull_mastiff', 'Tibetan_mastiff', 'French_bulldog', 'Great_Dane', 'Saint_Bernard', 'Eskimo_dog', 'malamute',
  'Siberian_husky', 'dalmatian', 'affenpinscher', 'basenji', 'pug', 'Leonberg', 'Newfoundland', 'Great_Pyrenees', 'Samoyed', 'Pomeranian',
  'chow', 'keeshond', 'Brabancon_griffon', 'Pembroke', 'Cardigan', 'toy_poodle', 'miniature_poodle', 'standard_poodle', 'Mexican_hairless', 'timber_wolf',
  'white_wolf', 'red_wolf', 'coyote', 'dingo', 'dhole', 'African_hunting_dog', 'hyena', 'red_fox', 'kit_fox', 'Arctic_fox',
  'grey_fox', 'tabby', 'tiger_cat', 'Persian_cat', 'Siamese_cat', 'Egyptian_cat', 'cougar', 'lynx', 'leopard', 'snow_leopard',
  'jaguar', 'lion', 'tiger', 'cheetah', 'brown_bear', 'American_black_bear', 'ice_bear', 'sloth_bear', 'mongoose', 'meerkat',
  'tiger_beetle', 'ladybug', 'ground_beetle', 'long-horned_beetle', 'leaf_beetle', 'dung_beetle', 'rhinoceros_beetle', 'weevil', 'fly', 'bee',
  'ant', 'grasshopper', 'cricket', 'walking_stick', 'cockroach', 'mantis', 'cicada', 'leafhopper', 'lacewing', 'dragonfly',
  'damselfly', 'admiral', 'ringlet', 'monarch', 'cabbage_butterfly', 'sulphur_butterfly', 'lycaenid', 'starfish', 'sea_urchin', 'sea_cucumber',
  'wood_rabbit', 'hare', 'Angora', 'hamster', 'porcupine', 'fox_squirrel', 'marmot', 'beaver', 'guinea_pig', 'sorrel',
  'zebra', 'hog', 'wild_boar', 'warthog', 'hippopotamus', 'ox', 'water_buffalo', 'bison', 'ram', 'bighorn',
  'ibex', 'hartebeest', 'impala', 'gazelle', 'Arabian_camel', 'llama', 'weasel', 'mink', 'polecat', 'black-footed_ferret',
  'otter', 'skunk', 'badger', 'armadillo', 'three-toed_sloth', 'orangutan', 'gorilla', 'chimpanzee', 'gibbon', 'siamang',
  'guenon', 'patas', 'baboon', 'macaque', 'langur', 'colobus', 'proboscis_monkey', 'marmoset', 'capuchin', 'howler_monkey',
  'titi', 'spider_monkey', 'squirrel_monkey', 'Madagascar_cat', 'indri', 'Indian_elephant', 'African_elephant', 'lesser_panda', 'giant_panda', 'barracouta',
  'eel', 'coho', 'rock_beauty', 'anemone_fish', 'sturgeon', 'gar', 'lionfish', 'puffer', 'abacus', 'abaya',
  'academic_gown', 'accordion', 'acoustic_guitar', 'aircraft_carrier', 'airliner', 'airship', 'altar', 'ambulance', 'amphibian', 'analog_clock',
  'apiary', 'apron', 'ashcan', 'assault_rifle', 'backpack', 'bakery', 'balance_beam', 'balloon', 'ballpoint', 'Band_Aid',
  'banjo', 'bannister', 'barbell', 'barber_chair', 'barbershop', 'barn', 'barometer', 'barrel', 'barrow', 'baseball',
  'basketball', 'bassinet', 'bassoon', 'bathing_cap', 'bath_towel', 'bathtub', 'beach_wagon', 'beacon', 'beaker', 'bearskin',
  'beer_bottle', 'beer_glass', 'bell_cote', 'bib', 'bicycle-built-for-two', 'bikini', 'binder', 'binoculars', 'birdhouse', 'boathouse',
  'bobsled', 'bolo_tie', 'bonnet', 'bookcase', 'bookshop', 'bottlecap', 'bow', 'bow_tie', 'brass', 'brassiere',
  'breakwater', 'breastplate', 'broom', 'bucket', 'buckle', 'bulletproof_vest', 'bullet_train', 'butcher_shop', 'cab', 'caldron',
  'candle', 'cannon', 'canoe', 'can_opener', 'cardigan', 'car_mirror', 'carousel', "carpenter's_kit", 'carton', 'car_wheel',
  'cash_machine', 'cassette', 'cassette_player', 'castle', 'catamaran', 'CD_player', 'cello', 'cellular_telephone', 'chain', 'chainlink_fence',
  'chain_mail', 'chain_saw', 'chest', 'chiffonier', 'chime', 'china_cabinet', 'Christmas_stocking', 'church', 'cinema', 'cleaver',
  'cliff_dwelling', 'cloak', 'clog', 'cocktail_shaker', 'coffee_mug', 'coffeepot', 'coil', 'combination_lock', 'computer_keyboard', 'confectionery',
  'container_ship', 'convertible', 'corkscrew', 'cornet', 'cowboy_boot', 'cowboy_hat', 'cradle', 'crane_(machine)', 'crash_helmet', 'crate',
  'crib', 'Crock_Pot', 'croquet_ball', 'crutch', 'cuirass', 'dam', 'desk', 'desktop_computer', 'dial_telephone', 'diaper',
  'digital_clock', 'digital_watch', 'dining_table', 'dishrag', 'dishwasher', 'disk_brake', 'dock', 'dogsled', 'dome', 'doormat',
  'drilling_platform', 'drum', 'drumstick', 'dumbbell', 'Dutch_oven', 'electric_fan', 'electric_guitar', 'electric_locomotive', 'entertainment_center', 'envelope',
  'espresso_maker', 'face_powder', 'feather_boa', 'file', 'fireboat', 'fire_engine', 'fire_screen', 'flagpole', 'flute', 'folding_chair',
  'football_helmet', 'forklift', 'fountain', 'fountain_pen', 'four-poster', 'freight_car', 'French_horn', 'frying_pan', 'fur_coat', 'garbage_truck',
  'gasmask', 'gas_pump', 'goblet', 'go-kart', 'golf_ball', 'golfcart', 'gondola', 'gong', 'gown', 'grand_piano',
  'greenhouse', 'grille', 'grocery_store', 'guillotine', 'hair_slide', 'hair_spray', 'half_track', 'hammer', 'hamper', 'hand_blower',
  'hand-held_computer', 'handkerchief', 'hard_disc', 'harmonica', 'harp', 'harvester', 'hatchet', 'holster', 'home_theater', 'honeycomb',
  'hook', 'hoopskirt', 'horizontal_bar', 'horse_cart', 'hourglass', 'iPod', 'iron', "jack-o'-lantern", 'jean', 'jeep',
  'jersey', 'jigsaw_puzzle', 'jinrikisha', 'joystick', 'kimono', 'knee_pad', 'knot', 'lab_coat', 'ladle', 'lampshade',
  'laptop', 'lawn_mower', 'lens_cap', 'letter_opener', 'library', 'lifeboat', 'lighter', 'limousine', 'liner', 'lipstick',
  'Loafer', 'lotion', 'loudspeaker', 'loupe', 'lumbermill', 'magnetic_compass', 'mailbag', 'mailbox', 'maillot_(tights)', 'maillot_(tank_suit)',
  'manhole_cover', 'maraca', 'marimba', 'mask', 'matchstick', 'maypole', 'maze', 'measuring_cup', 'medicine_chest', 'megalith',
  'microphone', 'microwave', 'military_uniform', 'milk_can', 'minibus', 'miniskirt', 'minivan', 'missile', 'mitten', 'mixing_bowl',
  'mobile_home', 'Model_T', 'modem', 'monastery', 'monitor', 'moped', 'mortar', 'mortarboard', 'mosque', 'mosquito_net',
  'motor_scooter', 'mountain_bike', 'mountain_tent', 'mouse', 'mousetrap', 'moving_van', 'muzzle', 'nail', 'neck_brace', 'necklace',
  'nipple', 'notebook', 'obelisk', 'oboe', 'ocarina', 'odometer', 'oil_filter', 'organ', 'oscilloscope', 'overskirt',
  'oxcart', 'oxygen_mask', 'packet', 'paddle', 'paddlewheel', 'padlock', 'paintbrush', 'pajama', 'palace', 'panpipe',
  'paper_towel', 'parachute', 'parallel_bars', 'park_bench', 'parking_meter', 'passenger_car', 'patio', 'pay-phone', 'pedestal', 'pencil_box',
  'pencil_sharpener', 'perfume', 'Petri_dish', 'photocopier', 'pick', 'pickelhaube', 'picket_fence', 'pickup', 'pier', 'piggy_bank',
  'pill_bottle', 'pillow', 'ping-pong_ball', 'pinwheel', 'pirate', 'pitcher', 'plane', 'planetarium', 'plastic_bag', 'plate_rack',
  'plow', 'plunger', 'Polaroid_camera', 'pole', 'police_van', 'poncho', 'pool_table', 'pop_bottle', 'pot', "potter's_wheel",
  'power_drill', 'prayer_rug', 'printer', 'prison', 'projectile', 'projector', 'puck', 'punching_bag', 'purse', 'quill',
  'quilt', 'racer', 'racket', 'radiator', 'radio', 'radio_telescope', 'rain_barrel', 'recreational_vehicle', 'reel', 'reflex_camera',
  'refrigerator', 'remote_control', 'restaurant', 'revolver', 'rifle', 'rocking_chair', 'rotisserie', 'rubber_eraser', 'rugby_ball', 'rule',
  'running_shoe', 'safe', 'safety_pin', 'saltshaker', 'sandal', 'sarong', 'sax', 'scabbard', 'scale', 'school_bus',
  'schooner', 'scoreboard', 'screen', 'screw', 'screwdriver', 'seat_belt', 'sewing_machine', 'shield', 'shoe_shop', 'shoji',
  'shopping_basket', 'shopping_cart', 'shovel', 'shower_cap', 'shower_curtain', 'ski', 'ski_mask', 'sleeping_bag', 'slide_rule', 'sliding_door',
  'slot', 'snorkel', 'snowmobile', 'snowplow', 'soap_dispenser', 'soccer_ball', 'sock', 'solar_dish', 'sombrero', 'soup_bowl',
  'space_bar', 'space_heater', 'space_shuttle', 'spatula', 'speedboat', 'spider_web', 'spindle', 'sports_car', 'spotlight', 'stage',
  'steam_locomotive', 'steel_arch_bridge', 'steel_drum', 'stethoscope', 'stole', 'stone_wall', 'stopwatch', 'stove', 'strainer', 'streetcar',
  'stretcher', 'studio_couch', 'stupa', 'submarine', 'suit', 'sundial', 'sunglass', 'sunglasses', 'sunscreen', 'suspension_bridge',
  'swab', 'sweatshirt', 'swimming_trunks', 'swing', 'switch', 'syringe', 'table_lamp', 'tank', 'tape_player', 'teapot',
  'teddy', 'television', 'tennis_ball', 'thatch', 'theater_curtain', 'thimble', 'thresher', 'throne', 'tile_roof', 'toaster',
  'tobacco_shop', 'toilet_seat', 'torch', 'totem_pole', 'tow_truck', 'toyshop', 'tractor', 'trailer_truck', 'tray', 'trench_coat',
  'tricycle', 'trimaran', 'tripod', 'triumphal_arch', 'trolleybus', 'trombone', 'tub', 'turnstile', 'typewriter_keyboard', 'umbrella',
  'unicycle', 'upright', 'vacuum', 'vase', 'vault', 'velvet', 'vending_machine', 'vestment', 'viaduct', 'violin',
  'volleyball', 'waffle_iron', 'wall_clock', 'wallet', 'wardrobe', 'warplane', 'washbasin', 'washer', 'water_bottle', 'water_jug',
  'water_tower', 'whiskey_jug', 'whistle', 'wig', 'window_screen', 'window_shade', 'Windsor_tie', 'wine_bottle', 'wing', 'wok',
  'wooden_spoon', 'wool', 'worm_fence', 'wreck', 'yawl', 'yurt', 'web_site', 'comic_book', 'crossword_puzzle', 'street_sign',
  'traffic_light', 'book_jacket', 'menu', 'plate', 'guacamole', 'consomme', 'hot_pot', 'trifle', 'ice_cream', 'ice_lolly',
  'French_loaf', 'bagel', 'pretzel', 'cheeseburger', 'hotdog', 'mashed_potato', 'head_cabbage', 'broccoli', 'cauliflower', 'zucchini',
  'spaghetti_squash', 'acorn_squash', 'butternut_squash', 'cucumber', 'artichoke', 'bell_pepper', 'cardoon', 'mushroom', 'Granny_Smith', 'strawberry',
  'orange', 'lemon', 'fig', 'pineapple', 'banana', 'jackfruit', 'custard_apple', 'pomegranate', 'hay', 'carbonara',
  'chocolate_sauce', 'dough', 'meat_loaf', 'pizza', 'potpie', 'burrito', 'red_wine', 'espresso', 'cup', 'eggnog',
  'alp', 'bubble', 'cliff', 'coral_reef', 'geyser', 'lakeside', 'promontory', 'sandbar', 'seashore', 'valley',
  'volcano', 'ballplayer', 'groom', 'scuba_diver', 'rapeseed', 'daisy', "yellow_lady's_slipper", 'corn', 'acorn', 'hip',
  'buckeye', 'coral_fungus', 'agaric', 'gyromitra', 'stinkhorn', 'earthstar', 'hen-of-the-woods', 'bolete', 'ear', 'toilet_tissue'
)

#' Get Base URL for Pre-Trained Neural Models
#'
#' Returns the base URL for downloading pre-trained ONNX models from GitHub Releases.
#' Can be customized globally via \code{options(pliman.models_base_url = "...")}.
#'
#' @return Character URL string ending in a slash.
#' @export
pliman_models_base_url <- function() {
  getOption("pliman.models_base_url", "https://github.com/NEPEM-UFSC/models/releases/download/v1.0.0/")
}

# Internal helper for robust downloads with SSL error recovery
.download_url_robust <- function(url, dest_file, min_size = 1000) {
  dl_success <- FALSE
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    tryCatch({
      curl::curl_download(url, destfile = dest_file, quiet = FALSE, mode = "wb", handle = h)
      dl_success <- file.exists(dest_file) && file.info(dest_file)$size >= min_size
    }, error = function(e) {
      if (grepl("certificate|schannel|SEC_E|SSL", e$message, ignore.case = TRUE)) {
        cli::cli_alert_warning("SSL handshake failed ({e$message}). Retrying with SSL verification bypass...")
        h_insecure <- curl::new_handle()
        curl::handle_setopt(h_insecure, ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)
        tryCatch({
          curl::curl_download(url, destfile = dest_file, quiet = FALSE, mode = "wb", handle = h_insecure)
          dl_success <<- file.exists(dest_file) && file.info(dest_file)$size >= min_size
        }, error = function(e2) NULL)
      }
    })
  }

  if (!dl_success) {
    tryCatch({
      utils::download.file(url, destfile = dest_file, mode = "wb", quiet = FALSE)
      dl_success <- file.exists(dest_file) && file.info(dest_file)$size >= min_size
    }, error = function(e) NULL)
  }

  if (!isTRUE(dl_success) || !file.exists(dest_file) || file.info(dest_file)$size < min_size) {
    if (file.exists(dest_file)) unlink(dest_file)
    return(FALSE)
  }
  return(TRUE)
}

#' List Available Pre-Trained Neural Models
#'
#' Lists all supported background removal and salient object detection models,
#' their disk sizes, input resolutions, description, and local download status.
#'
#' @param dir Directory where models are stored (default: `pliman_model_dir()`).
#' @return A data frame containing model metadata and local availability.
#' @export
#' @examples
#' pliman_available_models()
pliman_available_models <- function(dir = pliman_model_dir()) {
  base_url <- pliman_models_base_url()
  hf_yolo_url <- "https://huggingface.co/zwh20081/yolo26-onnx/resolve/main/"
  models_info <- list(
    list(
      name = "u2netp",
      description = "U2-Net Portable: Ultra-lightweight (~4.6 MB) and fast for CPU",
      size_mb = 4.6,
      input_size = 320,
      filename = "u2netp.onnx",
      url = paste0(base_url, "u2netp.onnx")
    ),
    list(
      name = "birefnet-lite",
      description = "BiRefNet Lite: Bilateral Reference Network for High-Resolution Dichotomous Segmentation (1024x1024)",
      size_mb = 213.6,
      input_size = 1024,
      filename = "birefnet-lite.onnx",
      url = paste0(base_url, "birefnet-lite.onnx")
    ),
    list(
      name = "isnet-general-use",
      description = "IS-Net: High-resolution (1024x1024) boundary matting (best for leaf contours & lesions)",
      size_mb = 178.6,
      input_size = 1024,
      filename = "isnet-general-use.onnx",
      url = paste0(base_url, "isnet-general-use.onnx")
    ),
    list(
      name = "rmbg-1.4",
      description = "BRIA RMBG 1.4: State-of-the-Art background removal (1024x1024, alias 'rmbg')",
      size_mb = 176.1,
      input_size = 1024,
      filename = "rmbg-1.4.onnx",
      url = paste0(base_url, "rmbg-1.4.onnx")
    ),
    list(
      name = "rmbg-2.0",
      description = "BRIA RMBG 2.0: Next-Gen BiRefNet-based background removal (1024x1024)",
      size_mb = 976.9,
      input_size = 1024,
      filename = "rmbg-2.0.onnx",
      url = paste0(base_url, "rmbg-2.0.onnx")
    ),
    list(
      name = "silueta",
      description = "Silueta: Compact model optimized for high accuracy and low latency (320x320)",
      size_mb = 44.2,
      input_size = 320,
      filename = "silueta.onnx",
      url = paste0(base_url, "silueta.onnx")
    ),
    list(
      name = "u2net",
      description = "U2-Net: Full salient object detection model (320x320)",
      size_mb = 176.3,
      input_size = 320,
      filename = "u2net.onnx",
      url = paste0(base_url, "u2net.onnx")
    ),
    list(
      name = "withoutbg",
      description = "withoutBG Open Weights: DepthAnythingV2 + ConvNeXt-fused U-Net matting (448x448)",
      size_mb = 433.4,
      input_size = 448,
      filename = "withoutbg.onnx",
      url = paste0(base_url, "withoutbg.onnx")
    ),
    list(
      name = "sam2.1",
      description = "Segment Anything Model 2.1 (Hiera-Tiny, Meta AI): Foundation image segmentation (1024x1024)",
      size_mb = 120.2,
      input_size = 1024,
      filename = "sam2.1.encoder.onnx",
      url = paste0(base_url, "sam2.1.encoder.onnx")
    ),
    list(
      name = "sam3.1",
      description = "Segment Anything Model 3.1 (Meta AI): Concept-driven foundation segmentation (1024x1024)",
      size_mb = 868.1,
      input_size = 1024,
      filename = "sam3.1.onnx",
      url = paste0(base_url, "sam3.1.onnx")
    ),
    list(
      name = "ben2",
      description = "Boundary-aware Extraction Network (BEN2 Base, Prama LLC): High-accuracy crisp boundary extraction (1024x1024)",
      size_mb = 212.6,
      input_size = 1024,
      filename = "ben2.onnx",
      url = paste0(base_url, "ben2.onnx")
    ),
    list(
      name = "grounded-sam",
      description = "Grounded-SAM: Grounding DINO Tiny + SAM 2.1 zero-shot text-prompted instance segmentation (800x800 & 1024x1024)",
      size_mb = 314.8,
      input_size = 800,
      filename = "groundingdino-tiny.onnx",
      url = paste0(base_url, "groundingdino-tiny.onnx")
    ),
    list(
      name = "persam",
      description = "PerSAM: One-shot visual exemplar instance segmentation and counting via SAM 2.1 embeddings (1024x1024)",
      size_mb = 120.2,
      input_size = 1024,
      filename = "sam2.1.encoder.onnx",
      url = paste0(base_url, "sam2.1.encoder.onnx")
    ),
    list(
      name = "depth-anything-v2",
      description = "Depth Anything V2 (Small): Monocular relative 3D depth estimation (518x518)",
      size_mb = 95.5,
      input_size = 518,
      filename = "depth-anything-v2-small.onnx",
      url = paste0(base_url, "depth-anything-v2-small.onnx")
    ),
    list(
      name = "dinov2",
      description = "DINOv2 (ViT-S/14, Meta AI): Self-supervised dense foundation vision transformer (518x518)",
      size_mb = 85.2,
      input_size = 518,
      filename = "dinov2-vits14.onnx",
      url = paste0(base_url, "dinov2-vits14.onnx")
    ),
    # YOLO26 Object Detection (Ultralytics)
    list(
      name = "yolo26n",
      description = "YOLO26 Nano Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 9.5,
      input_size = 640,
      filename = "yolo26n.onnx",
      url = paste0(hf_yolo_url, "yolo26n.onnx")
    ),
    list(
      name = "yolo26s",
      description = "YOLO26 Small Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 36.5,
      input_size = 640,
      filename = "yolo26s.onnx",
      url = paste0(hf_yolo_url, "yolo26s.onnx")
    ),
    list(
      name = "yolo26m",
      description = "YOLO26 Medium Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 78.2,
      input_size = 640,
      filename = "yolo26m.onnx",
      url = paste0(hf_yolo_url, "yolo26m.onnx")
    ),
    list(
      name = "yolo26l",
      description = "YOLO26 Large Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 95.0,
      input_size = 640,
      filename = "yolo26l.onnx",
      url = paste0(hf_yolo_url, "yolo26l.onnx")
    ),
    list(
      name = "yolo26x",
      description = "YOLO26 Extra-Large Detection: High-accuracy general object detection (640x640, 80 COCO classes)",
      size_mb = 212.9,
      input_size = 640,
      filename = "yolo26x.onnx",
      url = paste0(hf_yolo_url, "yolo26x.onnx")
    ),

    # YOLO26 Instance Segmentation (Ultralytics)
    list(
      name = "yolo26n-seg",
      description = "YOLO26 Nano Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 10.7,
      input_size = 640,
      filename = "yolo26n-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26n-seg.onnx")
    ),
    list(
      name = "yolo26s-seg",
      description = "YOLO26 Small Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 40.0,
      input_size = 640,
      filename = "yolo26s-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26s-seg.onnx")
    ),
    list(
      name = "yolo26m-seg",
      description = "YOLO26 Medium Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 90.2,
      input_size = 640,
      filename = "yolo26m-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26m-seg.onnx")
    ),
    list(
      name = "yolo26l-seg",
      description = "YOLO26 Large Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 107.1,
      input_size = 640,
      filename = "yolo26l-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26l-seg.onnx")
    ),
    list(
      name = "yolo26x-seg",
      description = "YOLO26 Extra-Large Segmentation: High-accuracy instance segmentation (640x640, 80 COCO classes)",
      size_mb = 240.0,
      input_size = 640,
      filename = "yolo26x-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26x-seg.onnx")
    ),

    # YOLO26 Pose Estimation (Ultralytics)
    list(
      name = "yolo26n-pose",
      description = "YOLO26 Nano Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 11.6,
      input_size = 640,
      filename = "yolo26n-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26n-pose.onnx")
    ),
    list(
      name = "yolo26s-pose",
      description = "YOLO26 Small Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 39.9,
      input_size = 640,
      filename = "yolo26s-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26s-pose.onnx")
    ),
    list(
      name = "yolo26m-pose",
      description = "YOLO26 Medium Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 82.6,
      input_size = 640,
      filename = "yolo26m-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26m-pose.onnx")
    ),
    list(
      name = "yolo26l-pose",
      description = "YOLO26 Large Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 99.4,
      input_size = 640,
      filename = "yolo26l-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26l-pose.onnx")
    ),
    list(
      name = "yolo26x-pose",
      description = "YOLO26 Extra-Large Pose: High-accuracy human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 220.0,
      input_size = 640,
      filename = "yolo26x-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26x-pose.onnx")
    ),

    # YOLO26 Image Classification (Ultralytics)
    list(
      name = "yolo26n-cls",
      description = "YOLO26 Nano Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 10.8,
      input_size = 640,
      filename = "yolo26n-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26n-cls.onnx")
    ),
    list(
      name = "yolo26s-cls",
      description = "YOLO26 Small Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 25.7,
      input_size = 640,
      filename = "yolo26s-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26s-cls.onnx")
    ),
    list(
      name = "yolo26m-cls",
      description = "YOLO26 Medium Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 44.4,
      input_size = 640,
      filename = "yolo26m-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26m-cls.onnx")
    ),
    list(
      name = "yolo26l-cls",
      description = "YOLO26 Large Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 53.9,
      input_size = 640,
      filename = "yolo26l-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26l-cls.onnx")
    ),
    list(
      name = "stardist",
      description = "StarDist (DSB 2018): Star-convex polygon detection for round/overlapping objects & cells",
      size_mb = 34.8,
      input_size = 256,
      filename = "stardist-dsb2018.onnx",
      url = paste0(base_url, "stardist-dsb2018.onnx")
    ),
    list(
      name = "realesrgan-compact",
      description = "Real-ESRGAN Compact: Fast 4x generative super-resolution with edge preservation",
      size_mb = 16.7,
      input_size = 256,
      filename = "realesrgan-compact.onnx",
      url = paste0(base_url, "realesrgan-compact.onnx")
    )
  )

  res <- do.call(rbind, lapply(models_info, function(m) {
    file_path <- file.path(dir, m$filename)
    downloaded <- file.exists(file_path) && file.info(file_path)$size > 1000000
    data.frame(
      model = m$name,
      size_mb = m$size_mb,
      input_size = paste0(m$input_size, "x", m$input_size),
      downloaded = downloaded,
      description = m$description,
      filename = m$filename,
      url = m$url,
      path = ifelse(downloaded, file_path, NA_character_),
      stringsAsFactors = FALSE
    )
  }))

  return(res)
}

#' Download Pre-Trained ONNX Models
#'
#' Downloads one or more pre-trained ONNX models into the user's model directory.
#' Supports downloading individual models, a subset of models, or all models at once.
#'
#' @param model Character vector naming the model(s) to download. Options include:
#'   * `"u2netp"` (default): Ultra-lightweight (~4.4 MB, 320x320).
#'   * `"ben2"`: Boundary-aware Extraction Network (~212.6 MB, 1024x1024).
#'   * `"grounded-sam"`: Grounding DINO Tiny + SAM 2.1 (~314.8 MB, 800x800 & 1024x1024).
#'   * `"sam2.1"`: Segment Anything 2.1 (~120.2 MB, 1024x1024).
#'   * `"sam3.1"`: Segment Anything 3.1 (~868.1 MB, 1024x1024).
#'   * `"rmbg-1.4"` (or `"rmbg"`): BRIA RMBG 1.4 (~168.0 MB, 1024x1024).
#'   * `"rmbg-2.0"`: BRIA RMBG 2.0 (~976.9 MB, 1024x1024).
#'   * `"withoutbg"`: withoutBG Open Weights (~433.4 MB, 448x448).
#'   * `"birefnet-lite"`: BiRefNet Lite (~213.6 MB, 1024x1024).
#'   * `"isnet-general-use"`: High-resolution boundary matting (~170.4 MB, 1024x1024).
#'   * `"silueta"`: Compact silhouette model (~42.1 MB, 320x320).
#'   * `"u2net"`: Full U2-Net (~167.8 MB, 320x320).
#'   * `"all"`: Downloads all 12 models above at once.
#'   Defaults to `"u2netp"`.
#' @param dir Directory to save the model files (default: `pliman_model_dir()`).
#' @param force Logical. If `TRUE`, forces re-download even if the file exists (default `FALSE`).
#' @return A character vector with the absolute path(s) to the downloaded `.onnx` model file(s).
#' @export
#' @examples
#' \dontrun{
#'   # Download a single model
#'   pliman_download_model("u2netp")
#'
#'   # Download BRIA RMBG
#'   pliman_download_model("rmbg")
#'
#'   # Download Segment Anything 2.1
#'   pliman_download_model("sam2.1")
#'
#'   # Download withoutBG Open Weights
#'   pliman_download_model("withoutbg")
#'
#'   # Download all models at once
#'   pliman_download_model("all")
#'
#'   # Download specific multiple models
#'   pliman_download_model(c("u2netp", "birefnet-lite"))
#' }
pliman_download_model <- function(model = "u2netp",
                                  dir = pliman_model_dir(),
                                  force = FALSE) {
  dir <- pliman_model_dir(dir)
  models_df <- pliman_available_models(dir = dir)
  valid_models <- models_df$model

  if (identical(model, "all") || ("all" %in% model)) {
    model <- valid_models
  } else {
    model <- .resolve_model_name(model)
  }

  # Validate model names
  invalid <- setdiff(model, valid_models)
  if (length(invalid) > 0) {
    stop("Unknown model(s): ", paste(invalid, collapse = ", "),
         "\nValid options are: 'all', ", paste(paste0("'", valid_models, "'"), collapse = ", "))
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(600, old_timeout))

  download_one <- function(m) {
    if (m %in% c("persam", "per-sam")) m <- "sam2.1"
    row <- models_df[models_df$model == m, ]
    base_url <- pliman_models_base_url()

    if (m == "sam2.1") {
      enc_file <- file.path(dir, "sam2.1.encoder.onnx")
      dec_file <- file.path(dir, "sam2.1.decoder.onnx")
      if (!isTRUE(force) && file.exists(enc_file) && file.info(enc_file)$size > 10000000 &&
          file.exists(dec_file) && file.info(dec_file)$size > 1000000) {
        return(enc_file)
      }
      cli::cli_alert_info("Downloading SAM 2.1 encoder and decoder ({row$size_mb} MB) to {.path {dir}}...")
      enc_url <- paste0(base_url, "sam2.1.encoder.onnx")
      dec_url <- paste0(base_url, "sam2.1.decoder.onnx")
      ok_enc <- .download_url_robust(enc_url, enc_file, min_size = 10000000)
      ok_dec <- .download_url_robust(dec_url, dec_file, min_size = 1000000)
      if (!ok_enc || !ok_dec) {
        stop("Failed to download SAM 2.1 model. Please check your internet connection.")
      }
      cli::cli_alert_success("Model {.val {m}} successfully downloaded!")
      return(enc_file)
    }

    if (m %in% c("grounded-sam", "grounding-dino")) {
      model_file <- file.path(dir, "groundingdino-tiny.onnx")
      vocab_file <- file.path(dir, "vocab.txt")
      if (!isTRUE(force) && file.exists(model_file) && file.info(model_file)$size > 100000000 &&
          file.exists(vocab_file) && file.info(vocab_file)$size > 50000) {
        return(model_file)
      }
      cli::cli_alert_info("Downloading Grounding DINO Tiny and BERT vocab ({row$size_mb} MB) to {.path {dir}}...")
      dino_url <- paste0(base_url, "groundingdino-tiny.onnx")
      vocab_url <- paste0(base_url, "vocab.txt")
      ok_dino <- .download_url_robust(dino_url, model_file, min_size = 100000000)
      ok_vocab <- .download_url_robust(vocab_url, vocab_file, min_size = 50000)
      if (!ok_dino || !ok_vocab) {
        stop("Failed to download Grounding DINO model/vocab. Please check your internet connection.")
      }
      cli::cli_alert_success("Grounding DINO model and vocab successfully downloaded!")
      return(model_file)
    }

    filename <- if (!is.null(row$filename) && length(row$filename) > 0 && !is.na(row$filename[1]) && nzchar(row$filename[1])) {
      row$filename[1]
    } else {
      paste0(m, ".onnx")
    }
    dest_file <- file.path(dir, filename)

    if (!isTRUE(force) && file.exists(dest_file) && file.info(dest_file)$size > 1000000) {
      return(dest_file)
    }

    url <- if (!is.null(row$url) && length(row$url) > 0 && !is.na(row$url[1]) && nzchar(row$url[1])) {
      row$url[1]
    } else {
      paste0(base_url, filename)
    }

    cli::cli_alert_info("Downloading model {.val {m}} ({row$size_mb} MB) to {.path {dir}}...")

    dl_success <- .download_url_robust(url, dest_file, min_size = 1000000)

    if (!isTRUE(dl_success) || !file.exists(dest_file) || file.info(dest_file)$size < 1000000) {
      if (file.exists(dest_file)) unlink(dest_file)
      stop("Failed to download model ", m, ". Please check your internet connection (e.g., Wi-Fi captive portal login).")
    }

    cli::cli_alert_success("Model {.val {m}} successfully downloaded!")
    return(dest_file)
  }

  res <- vapply(model, download_one, character(1))
  if (length(res) == 1L) {
    return(unname(res[1]))
  }
  return(res)
}

# ==============================================================================
# SECTION 2: ENVIRONMENT CONFIGURATION & DIRECT ONNX RUNTIME DOWNLOADER
# ==============================================================================

#' Directory for ONNX Runtime C++ Shared Libraries
#'
#' Gets or sets the local directory where the official Microsoft ONNX Runtime
#' C++ shared libraries (`onnxruntime.dll`, `libonnxruntime.so`, or `libonnxruntime.dylib`)
#' are saved.
#'
#' @param dir Optional custom directory path.
#' @return A character string with the path to the bin directory.
#' @export
#' @examples
#' pliman_onnx_dir()
pliman_onnx_dir <- function(dir = NULL) {
  if (!is.null(dir)) {
    dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(dir)
  }

  default_dir <- file.path(tools::R_user_dir("pliman", which = "data"), "bin")
  default_dir <- normalizePath(default_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(default_dir)) {
    dir.create(default_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(default_dir)
}

#' Download Microsoft ONNX Runtime C++ Library Directly
#'
#' Downloads and installs the official Microsoft ONNX Runtime C++ shared library
#' directly from Microsoft's GitHub releases, with **zero external package dependencies**
#' and **zero Python**.
#'
#' Automatically detects the host operating system (Windows, Linux, or macOS),
#' downloads the pre-compiled CPU binary (~15–60 MB), extracts the shared library
#' (`onnxruntime.dll` on Windows, `libonnxruntime.so` on Linux, `libonnxruntime.dylib` on macOS),
#' places it in `pliman`'s binary directory, and configures the environment.
#'
#' @param version Character string with the ONNX Runtime release version (default `"1.20.1"`).
#' @param dir Directory to store the library (default: `pliman_onnx_dir()`).
#' @param force Logical. If `TRUE`, forces re-download even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`. When `"gpu"` on Windows,
#'   downloads and configures `onnxruntime-directml.dll` and `DirectML.dll` from Microsoft's NuGet CDN.
#' @return A character string with the path to the installed shared library.
#' @export
#' @examples
#' \dontrun{
#'   pliman_download_onnx()
#'   pliman_download_onnx(engine = "gpu")
#' }
pliman_download_onnx <- function(version = "1.20.1",
                                 dir = pliman_onnx_dir(),
                                 force = FALSE,
                                 engine = c("cpu", "gpu")) {
  engine <- match.arg(engine)
  os <- tolower(Sys.info()[["sysname"]])

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  if (engine == "gpu") {
    if (os != "windows") {
      cli::cli_alert_warning("GPU acceleration via DirectML is currently only supported on Windows (DirectX 12). Falling back to CPU engine.")
      return(pliman_download_onnx(version = version, dir = dir, force = force, engine = "cpu"))
    }

    dest_dml_ort <- file.path(dir, "onnxruntime-directml.dll")
    dest_dml_core <- file.path(dir, "DirectML.dll")

    if (!isTRUE(force) &&
        file.exists(dest_dml_ort) && file.info(dest_dml_ort)$size > 1000000 &&
        file.exists(dest_dml_core) && file.info(dest_dml_core)$size > 1000000) {
      return(normalizePath(dest_dml_ort, winslash = "/"))
    }

    # Fast direct download of pre-extracted DLLs from repository
    base_url <- pliman_models_base_url()
    cli::cli_alert_info("Downloading DirectML (GPU) binaries from pliman repository...")
    ok_ort <- .download_url_robust(paste0(base_url, "onnxruntime-directml.dll"), dest_dml_ort, min_size = 1000000)
    ok_dml <- .download_url_robust(paste0(base_url, "DirectML.dll"), dest_dml_core, min_size = 1000000)

    if (!isTRUE(ok_ort) || !isTRUE(ok_dml)) {
      stop("Failed to download DirectML GPU libraries from ", base_url,
           ". Please check your connection or repository availability.")
    }

    cli::cli_alert_success("Microsoft ONNX Runtime DirectML (GPU) installed at {.path {dest_dml_ort}}!")
    return(normalizePath(dest_dml_ort, winslash = "/"))
  }

  lib_name <- switch(os,
    "windows" = "onnxruntime.dll",
    "linux"   = "libonnxruntime.so",
    "darwin"  = "libonnxruntime.dylib",
    stop("Unsupported operating system: ", os)
  )

  dest_lib <- file.path(dir, lib_name)
  onnxr_lib_dir <- file.path(tools::R_user_dir("onnxr", which = "data"), "lib")
  dest_onnxr_lib <- file.path(onnxr_lib_dir, lib_name)

  if (!isTRUE(force) && file.exists(dest_lib) && file.info(dest_lib)$size > 1000000) {
    if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(dest_onnxr_lib) || file.info(dest_onnxr_lib)$size < 1000000) {
      file.copy(dest_lib, dest_onnxr_lib, overwrite = TRUE)
    }
    Sys.setenv(ORT_ROOT = dir)
    return(normalizePath(dest_lib, winslash = "/"))
  }

  if (os == "windows") {
    base_url <- pliman_models_base_url()
    cli::cli_alert_info("Downloading ONNX Runtime CPU library from pliman repository...")
    ok_cpu <- .download_url_robust(paste0(base_url, "onnxruntime.dll"), dest_lib, min_size = 1000000)
    if (isTRUE(ok_cpu)) {
      if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(dest_lib, dest_onnxr_lib, overwrite = TRUE)
      Sys.setenv(ORT_ROOT = dir)
      cli::cli_alert_success("Microsoft ONNX Runtime CPU library installed at {.path {dest_lib}}!")
      return(normalizePath(dest_lib, winslash = "/"))
    }
  }

  cli::cli_alert_info("Downloading official Microsoft ONNX Runtime C++ library (v{version}) for {.val {os}}...")

  url <- switch(os,
    "windows" = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-win-x64-", version, ".zip"),
    "linux"   = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-linux-x64-", version, ".tgz"),
    "darwin"  = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-osx-universal-", version, ".tgz")
  )

  temp_archive <- file.path(tempdir(), paste0("ort_archive_", version, if (os == "windows") ".zip" else ".tgz"))
  old_timeout <- getOption("timeout")
  on.exit({
    options(timeout = old_timeout)
    if (file.exists(temp_archive)) unlink(temp_archive)
  }, add = TRUE)
  options(timeout = max(600, old_timeout))

  if (requireNamespace("curl", quietly = TRUE)) {
    curl::curl_download(url, destfile = temp_archive, quiet = FALSE, mode = "wb")
  } else {
    utils::download.file(url, destfile = temp_archive, mode = "wb", quiet = FALSE)
  }

  if (!file.exists(temp_archive) || file.info(temp_archive)$size < 100000) {
    stop("Failed to download ONNX Runtime from ", url, ". Please check your internet connection.")
  }

  cli::cli_alert_info("Extracting {.val {lib_name}}...")
  extract_tmp <- file.path(tempdir(), paste0("ort_extracted_", Sys.getpid()))
  dir.create(extract_tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(extract_tmp, recursive = TRUE), add = TRUE)

  if (os == "windows") {
    utils::unzip(temp_archive, exdir = extract_tmp)
  } else {
    utils::untar(temp_archive, exdir = extract_tmp)
  }

  found_files <- list.files(extract_tmp, pattern = paste0("^", gsub("\\.", "\\\\.", lib_name), "$"), recursive = TRUE, full.names = TRUE)
  if (length(found_files) == 0) {
    found_files <- list.files(extract_tmp, pattern = "onnxruntime", recursive = TRUE, full.names = TRUE)
    found_files <- found_files[grepl(paste0("\\", tools::file_ext(lib_name), "$"), found_files)]
  }

  if (length(found_files) == 0) {
    stop("Could not locate ", lib_name, " in the downloaded archive.")
  }

  src_lib <- found_files[1]
  file.copy(src_lib, dest_lib, overwrite = TRUE)

  tryCatch({
    if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(src_lib, dest_onnxr_lib, overwrite = TRUE)
  }, error = function(e) NULL, warning = function(w) NULL)

  Sys.setenv(ORT_ROOT = dir)
  cli::cli_alert_success("Microsoft ONNX Runtime C++ library installed at {.path {dest_lib}}!")
  return(normalizePath(dest_lib, winslash = "/"))
}

#' Proprietary ONNX Runtime Installer for pliman
#'
#' Downloads and configures the official Microsoft ONNX Runtime C++ shared library
#' directly into pliman's local bin directory, requiring **zero external packages**
#' and **zero Python**.
#'
#' @param version Character string with the ONNX Runtime release version (default `"1.20.1"`).
#' @param force Logical. If `TRUE`, forces re-download even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return A character string with the path to the installed shared library.
#' @export
#' @examples
#' \dontrun{
#'   onnx_install()
#'   onnx_install(engine = "gpu")
#' }
onnx_install <- function(version = "1.20.1", force = FALSE, engine = c("cpu", "gpu")) {
  pliman_download_onnx(version = version, force = force, engine = engine)
}

#' @rdname onnx_install
#' @export
pliman_install_onnx <- function(version = "1.20.1", force = FALSE, engine = c("cpu", "gpu")) {
  onnx_install(version = version, force = force, engine = engine)
}

#' Locate ONNX Runtime C++ Shared Library
#'
#' Finds the local path to the ONNX Runtime dynamic library (`onnxruntime.dll` on Windows,
#' `libonnxruntime.so` on Linux, `libonnxruntime.dylib` on macOS), or `onnxruntime-directml.dll`
#' when `engine = "gpu"` on Windows.
#'
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return Character path to the dynamic library or `NULL` if not found.
#' @export
pliman_onnx_lib_path <- function(engine = c("cpu", "gpu")) {
  engine <- match.arg(engine)
  os <- tolower(Sys.info()[["sysname"]])

  if (os == "windows") {
    p_gpu <- file.path(pliman_onnx_dir(), "onnxruntime-directml.dll")
    if (file.exists(p_gpu)) {
      return(normalizePath(p_gpu, winslash = "/"))
    }
    if (engine == "gpu") {
      return(NULL)
    }
  } else if (engine == "gpu") {
    cli::cli_alert_warning("DirectML GPU acceleration is only supported on Windows (DirectX 12). Falling back to 'cpu'.")
    engine <- "cpu"
  }

  lib_name <- switch(os,
    "windows" = "onnxruntime.dll",
    "linux"   = "libonnxruntime.so",
    "darwin"  = "libonnxruntime.dylib",
    "onnxruntime.dll"
  )

  p1 <- file.path(pliman_onnx_dir(), lib_name)
  if (file.exists(p1)) return(normalizePath(p1, winslash = "/"))

  ort_root <- Sys.getenv("ORT_ROOT", "")
  if (nzchar(ort_root)) {
    p2 <- file.path(ort_root, lib_name)
    if (file.exists(p2)) return(normalizePath(p2, winslash = "/"))
    p2b <- file.path(ort_root, "lib", lib_name)
    if (file.exists(p2b)) return(normalizePath(p2b, winslash = "/"))
  }

  p3 <- file.path(tools::R_user_dir("onnxr", which = "data"), "lib", lib_name)
  if (file.exists(p3)) return(normalizePath(p3, winslash = "/"))

  return(NULL)
}

pliman_onnx_library_path <- function(engine = c("cpu", "gpu")) {
  pliman_onnx_lib_path(engine = engine)
}


#' GPU Information and DirectML Status in pliman
#'
#' Queries the system's graphics hardware (via DirectX DXGI on Windows) to list
#' available GPU adapters, dedicated video memory (VRAM), and DirectML status.
#'
#' @return A data frame containing information on detected GPU adapters, or a message if no GPU is available.
#' @export
#' @examples
#' \dontrun{
#'   pliman_gpu_info()
#' }
pliman_gpu_info <- function() {
  info <- pliman_gpu_info_cpp()
  if (!isTRUE(info$available) || length(info$device_id) == 0) {
    cli::cli_alert_info("No DirectX 12 compatible GPU adapters detected or running on non-Windows OS.")
    return(invisible(NULL))
  }

  df <- data.frame(
    device_id = info$device_id,
    name = info$name,
    vram_mb = round(info$vram_mb, 1),
    is_dedicated = info$is_dedicated,
    is_default = (info$device_id == info$default_device_id),
    stringsAsFactors = FALSE
  )

  cli::cli_h1("Detected GPU Adapters (DirectX 12 / DirectML)")
  for (i in seq_len(nrow(df))) {
    prefix <- if (df$is_default[i]) "* (Default) " else "  "
    type_str <- if (df$is_dedicated[i]) "Dedicated" else "Integrated"
    cli::cli_alert_info("{prefix}[Device {df$device_id[i]}] {df$name[i]} ({type_str}, {df$vram_mb[i]} MB VRAM)")
  }

  dml_lib <- pliman_onnx_lib_path(engine = "gpu")
  if (!is.null(dml_lib) && file.exists(dml_lib)) {
    cli::cli_alert_success("DirectML GPU library: Installed at {.path {dml_lib}}")
  } else {
    cli::cli_alert_warning("DirectML GPU library: Not installed. Run {.code onnx_install(engine = 'gpu')} to configure.")
  }

  invisible(df)
}

#' Configure Complete Deep Learning Environment for pliman (100% Native C++)
#'
#' Sets up the entire zero-Python, zero-external-package Deep Learning environment
#' in a single call with interactive CLI feedback.
#'
#' This function:
#' 1. Creates and verifies the local model directory (`tools::R_user_dir("pliman", "data")/models`).
#' 2. Downloads and verifies the official Microsoft ONNX Runtime C++ shared library directly.
#' 3. Downloads the specified pre-trained models (`u2netp` and `isnet-general-use` by default).
#' 4. Runs a self-test inference using the native C++ engine to confirm 100% operational status.
#'
#' @param models Character vector of models to pre-download. Defaults to `"u2netp"`.
#'   Pass `"all"` to download all 5 models, or `NULL` to only configure the engine without downloading models.
#' @param dir Directory to store models (default: `pliman_model_dir()`).
#' @param force Logical. Force re-download and re-configuration even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return Logical `TRUE` invisibly on success.
#' @export
#' @examples
#' \dontrun{
#'   pliman_configure_dl()
#'   pliman_configure_dl(engine = "gpu")
#' }
pliman_configure_dl <- function(models = "u2netp",
                                dir = pliman_model_dir(),
                                force = FALSE,
                                engine = c("cpu", "gpu")) {
  engine <- match.arg(engine)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(600, old_timeout))

  cli::cli_h1("Configuring pliman Deep Learning Environment (100% Native C++)")

  # 1. Check and configure model directory
  cli::cli_h2("1. Model Storage Directory")
  model_path <- pliman_model_dir(dir)
  cli::cli_alert_success("Model directory: {.path {model_path}}")

  # 2. Check and configure ONNX Runtime C++ engine directly
  cli::cli_h2("2. Microsoft ONNX Runtime C++ Engine (Zero-Package Native)")
  lib_path <- onnx_install(force = force, engine = engine)
  cli::cli_alert_success("ONNX Runtime C++ library verified: {.path {lib_path}}")

  # 3. Download Requested Models
  if (!is.null(models)) {
    cli::cli_h2("3. Pre-trained Neural Models")
    if (identical(models, "all")) {
      models <- c("u2netp", "isnet-general-use", "rmbg-1.4", "silueta", "u2net")
    }

    for (m in models) {
      tryCatch({
        pliman_download_model(model = m, dir = model_path, force = force)
      }, error = function(e) {
        cli::cli_alert_danger("Failed to download model {.val {m}}: {e$message}")
      })
    }
  }

  # 4. Perform a quick verification test
  cli::cli_h2("4. Verification & Self-Test")
  test_model <- if ("u2netp" %in% models) "u2netp" else models[1]
  test_file <- file.path(model_path, paste0(test_model, ".onnx"))

  if (file.exists(test_file)) {
    cli::cli_alert_info("Running self-test with {.val {test_model}} in native C++...")
    test_ok <- tryCatch({
      dummy_tensor <- numeric(1 * 3 * 320 * 320)
      attr(dummy_tensor, "dims") <- as.integer(c(1, 3, 320, 320))
      res <- .run_onnx_inference(dummy_tensor, test_file, target_size = 320, engine = engine)
      if (is.matrix(res) && nrow(res) == 320 && ncol(res) == 320) {
        TRUE
      } else {
        FALSE
      }
    }, error = function(e) {
      cli::cli_alert_warning("Self-test error: {e$message}")
      FALSE
    })

    if (test_ok) {
      cli::cli_alert_success("Self-test passed! Native C++ neural inference is working at full speed.")
    }
  }

  cli::cli_h2("Configuration Complete!")
  cli::cli_alert_info("You can now run: {.code mask <- image_binary_dl(img, model = 'u2netp')}")
  cli::cli_alert_info("Or remove background: {.code img_tb <- image_remove_bg_dl(img, model = 'u2netp')}")
  return(invisible(TRUE))
}

# ==============================================================================
# SECTION 3: INTERNAL TENSOR PREPROCESSING & INFERENCE RUNNER
# ==============================================================================

# Bilinear interpolation resize for numeric/raw matrix via high-speed C++
.bilinear_resize_2d <- function(mat, out_w, out_h) {
  in_w <- nrow(mat)
  in_h <- ncol(mat)
  if (in_w == out_w && in_h == out_h) return(mat)
  image_resize_cpp(mat, as.integer(out_w), as.integer(out_h), filter = 1L)
}

# Fast polygon decimation for screen display (preserves full resolution in returned contours)
.decimate_poly <- function(p, max_pts = 2500L) {
  n <- nrow(p)
  if (n <= max_pts) return(p)
  step <- ceiling(n / max_pts)
  idx <- c(seq(1L, n, by = step), 1L)
  p[idx, , drop = FALSE]
}

.preprocess_nchw <- function(mat,
                             target_size = 320,
                             mean = c(0.485, 0.456, 0.406),
                             std = c(0.229, 0.224, 0.225),
                             letterbox = FALSE) {
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  # Extract R, G, B normalized [0, 1]
  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  dim(R) <- dim(G) <- dim(B) <- c(orig_w, orig_h)

  if (isTRUE(letterbox)) {
    scale <- target_size / max(orig_w, orig_h)
    new_w <- max(1L, round(orig_w * scale))
    new_h <- max(1L, round(orig_h * scale))

    R_res <- .bilinear_resize_2d(R, new_w, new_h)
    G_res <- .bilinear_resize_2d(G, new_w, new_h)
    B_res <- .bilinear_resize_2d(B, new_w, new_h)

    # Pad at top-left of canvas target_size x target_size [width, height]
    R_pad <- matrix(0.0, nrow = target_size, ncol = target_size)
    G_pad <- matrix(0.0, nrow = target_size, ncol = target_size)
    B_pad <- matrix(0.0, nrow = target_size, ncol = target_size)

    R_pad[1:new_w, 1:new_h] <- R_res
    G_pad[1:new_w, 1:new_h] <- G_res
    B_pad[1:new_w, 1:new_h] <- B_res

    R_norm <- (R_pad - mean[1]) / std[1]
    G_norm <- (G_pad - mean[2]) / std[2]
    B_norm <- (B_pad - mean[3]) / std[3]

    vec_R <- as.numeric(R_norm)
    vec_G <- as.numeric(G_norm)
    vec_B <- as.numeric(B_norm)

    flat_tensor <- c(vec_R, vec_G, vec_B)
    attr(flat_tensor, "dims") <- as.integer(c(1, 3, target_size, target_size))
    attr(flat_tensor, "letterbox_dims") <- c(new_w, new_h)
    return(flat_tensor)
  }

  # Resize each channel to target_size x target_size
  R_res <- .bilinear_resize_2d(R, target_size, target_size)
  G_res <- .bilinear_resize_2d(G, target_size, target_size)
  B_res <- .bilinear_resize_2d(B, target_size, target_size)

  # Normalize
  R_norm <- (R_res - mean[1]) / std[1]
  G_norm <- (G_res - mean[2]) / std[2]
  B_norm <- (B_res - mean[3]) / std[3]

  # In pliman, M[x, y] is column-major with x as row and y as col.
  # Reading column-by-column traverses y=1 (all x), then y=2 (all x),
  # which corresponds to the C-order NCHW row-major layout expected by ONNX.
  vec_R <- as.numeric(R_norm)
  vec_G <- as.numeric(G_norm)
  vec_B <- as.numeric(B_norm)

  # Combined NCHW flat vector [1, 3, target_size, target_size]
  flat_tensor <- c(vec_R, vec_G, vec_B)
  attr(flat_tensor, "dims") <- as.integer(c(1, 3, target_size, target_size))

  return(flat_tensor)
}

.run_onnx_inference <- function(tensor, model_path, target_size = 320, threads = 0,
                                engine = c("cpu", "gpu"), device_id = -1) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_path) || !file.exists(lib_path)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_path <- onnx_install(engine = engine)
  }

  tensor_dims <- attr(tensor, "dims")
  if (is.null(tensor_dims)) {
    tensor_dims <- as.integer(c(1, 3, target_size, target_size))
  }
  tensor_vec <- as.numeric(tensor)

  pred_mat <- run_onnx_inference_cpp(
    tensor_vec = tensor_vec,
    tensor_dims = tensor_dims,
    model_path = normalizePath(model_path, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_path, winslash = "/", mustWork = FALSE),
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  return(pred_mat)
}

# Internal core mask computation helper
.compute_dl_mask <- function(mat,
                             model = "u2netp",
                             threshold = 0.5,
                             fill_hull = TRUE,
                             filter = 0,
                             erode = 0,
                             dilate = 0,
                             opening = 0,
                             closing = 0,
                             min_area = 0,
                             invert = FALSE,
                             pick_object = FALSE,
                             prompt = NULL,
                             exemplar = FALSE,
                             threads = 0,
                             engine = c("cpu", "gpu"),
                             device_id = -1,
                             verbose = TRUE,
                             dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  model <- .resolve_model_name(model[1])
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Interactive picking of prompt point(s)
  is_text_prompt <- is.character(prompt) && length(prompt) >= 1 && !all(prompt %in% c("center", "box", "exemplar"))
  is_sam <- grepl("sam", model, ignore.case = TRUE)
  use_persam <- isTRUE(exemplar) || (model == "persam") || (is_sam && identical(prompt, "exemplar"))
  use_grounded_sam <- !use_persam && ((model == "grounded-sam") || (is_sam && is_text_prompt))

  if (use_persam) {
    p_res <- .run_persam(
      mat = mat,
      exemplar_points = if (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt)) prompt else NULL,
      threshold = threshold,
      threads = threads,
      engine = engine,
      device_id = device_id,
      fill_hull = fill_hull,
      filter = filter,
      erode = erode,
      dilate = dilate,
      opening = opening,
      closing = closing,
      min_area = min_area,
      invert = invert,
      mask = TRUE,
      verbose = verbose,
      dir = dir
    )
    if (isTRUE(verbose)) {
      .print_detection_summary(p_res$summary)
    }
    return(p_res$mask * 1.0)
  }

  if (use_grounded_sam) {
    if (is.null(prompt)) {
      cli::cli_abort("Model {.val grounded-sam} requires a text {.arg prompt} (e.g., {.code prompt = 'leaf'}).")
    }
    gs_res <- .run_grounded_sam(
      mat = mat,
      prompt = prompt,
      threshold = threshold,
      threads = threads,
      engine = engine,
      device_id = device_id,
      fill_hull = fill_hull,
      filter = filter,
      erode = erode,
      dilate = dilate,
      opening = opening,
      closing = closing,
      min_area = min_area,
      invert = invert,
      mask = TRUE,
      verbose = verbose,
      dir = dir
    )
    if (isTRUE(verbose)) {
      .print_detection_summary(gs_res$summary)
    }
    return(gs_res$mask * 1.0)
  }

  if (isTRUE(pick_object)) {
    cli::cli_alert_info("Click on the object(s) of interest in the plot window. Press <Esc> or right-click when finished.")
    plot(as_image(mat))
    pts <- tryCatch(graphics::locator(n = 512, type = "p", col = "red", pch = 19), error = function(e) NULL)
    if (!is.null(pts) && length(pts$x) > 0) {
      prompt <- cbind(pts$x, pts$y)
    } else {
      cli::cli_alert_warning("No point selected. Defaulting to center of image.")
      prompt <- c(orig_w / 2, orig_h / 2)
    }
  }

  if (is_sam) {
    if (is.null(prompt)) {
      prompt <- c(orig_w / 2, orig_h / 2)
    }

    # Ensure SAM models are downloaded
    model_file <- pliman_download_model(model = model, dir = dir)
    enc_path <- file.path(dir, "sam2.1.encoder.onnx")
    dec_path <- file.path(dir, "sam2.1.decoder.onnx")

    # SAM uses standard ImageNet normalization and 1024x1024
    tensor <- .preprocess_nchw(mat,
                               target_size = 1024L,
                               mean = c(0.485, 0.456, 0.406),
                               std = c(0.229, 0.224, 0.225),
                               letterbox = FALSE)

    # Convert prompt coordinates to 1024x1024 space
    if (is.character(prompt)) {
      if (prompt == "center") {
        pts_x <- 512.0
        pts_y <- 512.0
        pts_lbl <- 1L
      } else { # "box"
        pts_x <- c(10.0, 1014.0)
        pts_y <- c(10.0, 1014.0)
        pts_lbl <- c(2L, 3L)
      }
    } else if (is.matrix(prompt) || is.data.frame(prompt)) {
      pts_x <- (as.numeric(prompt[, 1]) / orig_w) * 1024.0
      pts_y <- (as.numeric(prompt[, 2]) / orig_h) * 1024.0
      pts_lbl <- rep(1L, length(pts_x))
    } else if (length(prompt) == 4) { # box: c(xmin, ymin, xmax, ymax)
      pts_x <- c(prompt[1] / orig_w * 1024.0, prompt[3] / orig_w * 1024.0)
      pts_y <- c(prompt[2] / orig_h * 1024.0, prompt[4] / orig_h * 1024.0)
      pts_lbl <- c(2L, 3L)
    } else { # single point: c(x, y)
      pts_x <- (as.numeric(prompt[1]) / orig_w) * 1024.0
      pts_y <- (as.numeric(prompt[2]) / orig_h) * 1024.0
      pts_lbl <- 1L
    }

    lib_path <- pliman_onnx_lib_path(engine = engine)
    if (is.null(lib_path) || !file.exists(lib_path)) {
      cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
      lib_path <- onnx_install(engine = engine)
    }

    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Running SAM 2.1 inference [{toupper(engine)}]...",
        msg_done = "SAM 2.1 inference complete"
      )
    }

    pred_raw <- run_sam2_inference_cpp(
      tensor_vec = as.numeric(tensor),
      points_x = pts_x,
      points_y = pts_y,
      point_labels = pts_lbl,
      encoder_path = normalizePath(enc_path, winslash = "/", mustWork = FALSE),
      decoder_path = normalizePath(dec_path, winslash = "/", mustWork = FALSE),
      lib_path = normalizePath(lib_path, winslash = "/", mustWork = FALSE),
      num_threads = as.integer(threads),
      use_gpu = use_gpu,
      device_id = as.integer(device_id)
    )

    mask_prob <- .bilinear_resize_2d(pred_raw, orig_w, orig_h)
  } else {
    # Ensure model is downloaded
    model_file <- pliman_download_model(model = model, dir = dir)

    # Determine target input size, normalization, and letterboxing
    is_wbg <- grepl("withoutbg", model)
    if (model %in% c("rmbg-1.4", "isnet-general-use")) {
      target_size <- 1024L
      mean_val <- c(0.5, 0.5, 0.5)
      std_val <- c(1.0, 1.0, 1.0)
      use_letterbox <- FALSE
    } else if (model == "ben2") {
      target_size <- 1024L
      mean_val <- c(0, 0, 0)
      std_val <- c(1, 1, 1)
      use_letterbox <- FALSE
    } else if (is_wbg) {
      target_size <- 448L
      mean_val <- c(0, 0, 0)
      std_val <- c(1, 1, 1)
      use_letterbox <- TRUE
    } else {
      models_df <- pliman_available_models(dir = dir)
      row <- models_df[models_df$model == model, ]
      target_size <- if (nrow(row) > 0) as.integer(strsplit(as.character(row$input_size), "x")[[1]][1]) else 320L
      mean_val <- c(0.485, 0.456, 0.406)
      std_val <- c(0.229, 0.224, 0.225)
      use_letterbox <- FALSE
    }

    # 1. Pre-process to NCHW normalized tensor (row-major flat vector)
    tensor <- .preprocess_nchw(mat,
                               target_size = target_size,
                               mean = mean_val,
                               std = std_val,
                               letterbox = use_letterbox)

    # 2. Run ONNX Inference
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Running neural inference with {.val {model}} [{toupper(engine)}]...",
        msg_done = "Inference with {.val {model}} complete"
      )
    }

    pred_raw <- .run_onnx_inference(tensor, model_file, target_size = target_size, threads = threads, engine = engine, device_id = device_id)

    # Normalize raw predictions to [0, 1] probability
    p_min <- min(pred_raw, na.rm = TRUE)
    p_max <- max(pred_raw, na.rm = TRUE)
    if (p_min < 0 || p_max > 1) {
      prob_map <- 1.0 / (1.0 + exp(-pred_raw))
    } else {
      prob_map <- pred_raw
    }

    # 3. Un-letterbox if applicable and resize probability map back to original image dimensions
    if (isTRUE(use_letterbox)) {
      lb_dims <- attr(tensor, "letterbox_dims")
      prob_crop <- prob_map[1:lb_dims[1], 1:lb_dims[2], drop = FALSE]
      mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
    } else {
      mask_prob <- .bilinear_resize_2d(prob_map, orig_w, orig_h)
    }
  }

  # 4. Threshold into initial binary mask
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing binary mask...",
      msg_done = "Mask post-processing complete"
    )
  }

  mask <- mask_prob >= threshold

  # 5. Morphological Post-Processing on the Mask
  # 5.1 Fill internal holes in the mask
  if (isTRUE(fill_hull)) {
    mask <- fill_holes_cpp(mask)
    mask <- (mask > 0)
  }

  # 5.2 Median filter smoothing to reduce noise/speckles
  if (filter > 0) {
    mask_arr <- as.array(mask)
    if (length(dim(mask_arr)) == 2) dim(mask_arr) <- c(orig_w, orig_h, 1L)
    f_res <- median_filter_binary_cpp(mask_arr, orig_w, orig_h, 1L, as.integer(filter))
    mask <- (f_res[, , 1] > 0)
  }

  # 5.3 Morphological erosion (shrinks mask boundaries)
  if (erode > 0) {
    mask <- erode_cpp(mask, raio = as.integer(erode))
    mask <- (mask > 0)
  }

  # 5.4 Morphological dilation (expands mask boundaries)
  if (dilate > 0) {
    mask <- dilate_cpp(mask, raio = as.integer(dilate))
    mask <- (mask > 0)
  }

  # 5.5 Morphological opening (erode then dilate)
  if (opening > 0) {
    mask <- erode_cpp(mask, raio = as.integer(opening))
    mask <- dilate_cpp(mask, raio = as.integer(opening))
    mask <- (mask > 0)
  }

  # 5.6 Morphological closing (dilate then erode)
  if (closing > 0) {
    mask <- dilate_cpp(mask, raio = as.integer(closing))
    mask <- erode_cpp(mask, raio = as.integer(closing))
    mask <- (mask > 0)
  }

  # 5.7 Filter small artifacts / retain largest object
  if (min_area <= 0) {
    lbls <- bwlabel_cpp(mask)
    if (max(lbls) > 0) {
      tbl <- table(lbls[lbls > 0])
      largest_id <- as.integer(names(which.max(tbl)))
      mask <- (lbls == largest_id)
    }
  } else {
    lbls <- bwlabel_cpp(mask)
    if (max(lbls) > 0) {
      tbl <- table(lbls[lbls > 0])
      keep_ids <- as.integer(names(tbl[tbl >= min_area]))
      mask <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
    }
  }

  # 5.8 Invert mask if requested
  if (isTRUE(invert)) {
    mask <- !mask
  }

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }

  return(mask)
}

# ==============================================================================
# SECTION 4: USER-FACING SEGMENTATION, BINARIZATION & BACKGROUND REMOVAL
# ==============================================================================

#' Deep Learning Binary Mask Extraction
#'
#' Extracts a pixel-accurate binary mask (foreground = 1, background = 0) from an image using
#' state-of-the-art Deep Learning neural networks executed in ONNX format via native C++ bindings.
#' Supports salient object detection, dichotomous image segmentation, high-resolution boundary matting,
#' foundation segmentation models (SAM 2.1 / SAM 3.1), zero-shot text-prompted instance segmentation (Grounded-SAM),
#' and one-shot visual exemplar segmentation (PerSAM).
#'
#' @param img An `image` object, 2D grayscale matrix, 3D color array (RGB/RGBA), or a `list`
#'   of images/arrays. When a list is passed, each element is automatically processed in batch,
#'   returning a corresponding list of binary masks.
#' @param model Character string specifying the pre-trained neural network architecture to use:
#'   * **General Salient Object Detection & Background Extraction:**
#'     - `"u2netp"` (default): Ultra-lightweight U2-Net Portable (~4.6 MB, 320x320 input). Highly recommended
#'       for rapid CPU batch processing on standard laptops without requiring a dedicated GPU.
#'     - `"ben2"` (or `"ben"`): Boundary-aware Extraction Network Base (~212.6 MB, 1024x1024 input). High-accuracy
#'       crisp boundary extraction, exceptional for fine leaf contours, intricate roots, thin petioles, and plant serrations.
#'     - `"rmbg-1.4"` (or `"rmbg"`): BRIA AI RMBG 1.4 (~176.1 MB, 1024x1024 input). Enterprise-grade salient foreground
#'       extractor trained on diverse commercial datasets, exceptionally robust against studio lighting, reflections, and shadows.
#'     - `"rmbg-2.0"`: BRIA AI RMBG 2.0 (~976.9 MB, 1024x1024 input). Next-generation BiRefNet-based model for maximum boundary
#'       detail, fine textures, and translucent edge matting.
#'     - `"isnet-general-use"`: Intermediate Supervision Network (~178.6 MB, 1024x1024 input). High-precision boundary matting,
#'       ideal for subtle leaf lesions, chlorotic halos, and fine botanical structures.
#'     - `"withoutbg"`: withoutBG Open Weights (~433.4 MB, 448x448 input). Combines DepthAnythingV2 3D depth features with
#'       ConvNeXt-fused U-Net matting, excelling in depth-cluttered backgrounds.
#'     - `"birefnet-lite"` (or `"birefnet"`): Bilateral Reference Network Lite (~213.6 MB, 1024x1024 input). Bilateral reference
#'       architecture specialized in dichotomous image segmentation and complex shapes.
#'     - `"silueta"`: Compact silhouette extraction model (~44.2 MB, 320x320 input). Fast inference with solid accuracy on edge/CPU hardware.
#'     - `"u2net"`: Full-depth U2-Net architecture (~176.3 MB, 320x320 input). Deep multi-scale salient object detection model.
#'   * **Zero-Shot & Interactive Foundation Models:**
#'     - `"grounded-sam"`: Zero-shot open-vocabulary text-prompted instance segmentation (~314.8 MB: Grounding DINO Tiny + SAM 2.1).
#'       Detects and segments specific concepts described in natural language (e.g., `prompt = "leaf"`, `prompt = "fruit"`, `prompt = "insect"`).
#'     - `"persam"`: Personalize Segment Anything via SAM 2.1 embeddings (~120.2 MB, 1024x1024 input). Segments all instances
#'       visually similar to one or more user-clicked exemplar objects using Multi-Prototype Cosine Similarity Max-Pooling.
#'     - `"sam2.1"` (or `"sam2"`): Meta AI Segment Anything 2.1 (Hiera-Tiny, ~120.2 MB, 1024x1024 input). Foundation vision transformer
#'       supporting point prompts, box prompts, or interactive clicking (`pick_object = TRUE`).
#'     - `"sam3.1"` (or `"sam3"`): Meta AI Segment Anything 3.1 (~868.1 MB, 1024x1024 input). High-capacity concept-driven
#'       foundation segmentation model.
#' @param threshold Numeric value in `[0, 1]` specifying the cutoff threshold applied to the model's
#'   continuous probability/logit output map. Default is `0.5`. Increasing this threshold (e.g., `0.7`) makes
#'   segmentation more conservative (avoiding background noise), whereas decreasing it (e.g., `0.3`) includes
#'   faint, thin, or translucent object edges.
#' @param fill_hull Logical. If `TRUE` (default), automatically fills internal holes and voids
#'   (such as specular glare or reflective spots) within the segmented foreground object using morphological
#'   hole reconstruction (`fill_holes_cpp()`).
#' @param filter Integer specifying the window radius for binary median filtering (`median_filter_binary_cpp()`).
#'   Default is `0` (disabled). A positive value (e.g., `2` or `3`) eliminates salt-and-pepper noise and
#'   smoothes jagged boundaries without blurring sharp object corners.
#' @param erode Integer specifying the radius for morphological erosion (`erode_cpp()`). Default is `0` (disabled).
#'   Shrinks foreground boundaries inward by the given pixel radius, useful for severing thin touching bridges
#'   between adjacent objects or removing outer halo artifacts.
#' @param dilate Integer specifying the radius for morphological dilation (`dilate_cpp()`). Default is `0` (disabled).
#'   Expands foreground boundaries outward by the given pixel radius, useful for restoring peripheral details
#'   or slightly undersized masks.
#' @param opening Integer specifying the radius for morphological opening (erosion followed by dilation).
#'   Default is `0` (disabled). Removes small stray islands and thin protrusions while preserving overall
#'   object geometry.
#' @param closing Integer specifying the radius for morphological closing (dilation followed by erosion).
#'   Default is `0` (disabled). Fuses narrow internal cracks, gaps, and small indentations along the perimeter.
#' @param min_area Integer minimum object size in pixels evaluated via connected component labeling (`bwlabel_cpp()`).
#'   If positive, discards all foreground connected components smaller than this area. If `0` (default),
#'   automatically isolates and retains only the single largest connected component in the image.
#' @param invert Logical. If `TRUE`, inverts the output binary mask such that foreground becomes `0` and
#'   background becomes `1`. Default is `FALSE`.
#' @param pick_object Logical. If `TRUE`, activates an interactive graphics prompt via `graphics::locator()`,
#'   allowing the user to directly click on the object(s) of interest in the plot window to guide SAM models.
#'   Press `<Esc>` or right-click when finished. Default is `FALSE`.
#' @param prompt Flexible prompt specification for SAM, Grounded-SAM, and PerSAM models:
#'   * `NULL` (default): Defaults to the center of the image `c(width / 2, height / 2)` for SAM models.
#'   * Character string (e.g., `"leaf"`, `"fruit"`, `"weed"`): Triggers Grounded-SAM open-vocabulary object detection.
#'   * Numeric vector `c(x, y)`: Coordinates of a single foreground point prompt or single visual exemplar.
#'   * Two-column numeric matrix or data frame `cbind(x, y)`: Specifies multiple positive point prompts or multiple exemplar prototypes.
#'   * Numeric vector of length 4 `c(xmin, ymin, xmax, ymax)`: Specifies a bounding box prompt in pixel coordinates.
#'   * `"center"`: Uses the image center coordinates `c(width / 2, height / 2)`.
#'   * `"box"`: Uses a bounding box spanning the full image canvas.
#' @param exemplar Logical. If `TRUE` (or if `model = "persam"`), activates one-shot visual exemplar segmentation (PerSAM)
#'   using SAM 2.1 embeddings. Default is `FALSE`.
#' @param threads Integer specifying the number of threads for parallel ONNX Runtime execution. Default is `0`,
#'   which automatically detects and utilizes all available logical CPU cores on the system.
#' @param engine Character specifying the backend execution engine: `"cpu"` (default) or `"gpu"`. On Windows,
#'   `"gpu"` utilizes Microsoft DirectML over DirectX 12 hardware (compatible with NVIDIA, AMD, and Intel GPUs
#'   without requiring CUDA toolkit installation).
#' @param device_id Integer GPU device adapter index. Default is `-1`, which automatically selects the high-performance
#'   discrete GPU adapter with the maximum dedicated video memory (VRAM).
#' @param verbose Logical. If `TRUE` (default), displays step-by-step progress feedback and timing indicators
#'   via `cli::cli_progress_step()`. Set to `FALSE` for silent execution in automated pipelines.
#' @param plot Logical. If `TRUE` (default), renders the resulting binary mask to the active graphics device.
#' @param dir Character string specifying the local directory where ONNX model files are stored and cached.
#'   Defaults to `pliman_model_dir()`.
#' @param ... Additional arguments passed down to internal plotting methods.
#'
#' @return A binary `image` object (Grayscale, values 0 and 1) or a list of binary `image` objects.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # Fast binarization on CPU (ultra-lightweight, 4.6 MB)
#'   mask <- image_binary_dl(img, model = "u2netp")
#'
#'   # Sub-pixel boundary segmentation for fine veins and serrations
#'   mask_ben <- image_binary_dl(img, model = "ben2")
#'
#'   # Interactive prompt: click directly on objects of interest
#'   mask_interactive <- image_binary_dl(img, model = "sam2.1", pick_object = TRUE)
#'
#'   # Open-vocabulary text-prompted binarization
#'   mask_weed <- image_binary_dl(img, model = "grounded-sam", prompt = "weed")
#'
#'   # One-shot visual exemplar segmentation (PerSAM)
#'   mask_exemplar <- image_binary_dl(img, model = "persam", prompt = c(150, 200))
#' }
image_binary_dl <- function(img,
                            model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                            threshold = 0.5,
                            fill_hull = TRUE,
                            filter = 0,
                            erode = 0,
                            dilate = 0,
                            opening = 0,
                            closing = 0,
                            min_area = 0,
                            invert = FALSE,
                            pick_object = FALSE,
                            prompt = NULL,
                            exemplar = FALSE,
                            threads = 0,
                            engine = c("cpu", "gpu"),
                            device_id = -1,
                            verbose = TRUE,
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_binary_dl(x, model = model, threshold = threshold,
                      fill_hull = fill_hull, filter = filter,
                      erode = erode, dilate = dilate,
                      opening = opening, closing = closing,
                      min_area = min_area, invert = invert,
                      pick_object = pick_object, prompt = prompt,
                      exemplar = exemplar,
                      threads = threads,
                      engine = engine, device_id = device_id,
                      verbose = verbose,
                      plot = FALSE, dir = dir)
    })
    return(res)
  }

  if (is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.null(dim(mat)) || length(dim(mat)) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  mask <- .compute_dl_mask(mat = mat,
                           model = model,
                           threshold = threshold,
                           fill_hull = fill_hull,
                           filter = filter,
                           erode = erode,
                           dilate = dilate,
                           opening = opening,
                           closing = closing,
                           min_area = min_area,
                           invert = invert,
                           pick_object = pick_object,
                           prompt = prompt,
                           exemplar = exemplar,
                           threads = threads,
                           engine = engine,
                           device_id = device_id,
                           verbose = verbose,
                           dir = dir)

  mask_img <- as_image(mask, colormode = "Grayscale")

  if (isTRUE(plot)) {
    plot(mask_img, ...)
  }

  invisible(mask_img)
}

#' Deep Learning Background Removal (Transparent / RGBA)
#'
#' Removes the background of an image using pre-trained Deep Learning models,
#' producing a 4-channel **RGBA** image with a transparent background by default,
#' or replacing the background with a solid color.
#'
#' @inheritParams image_binary_dl
#' @param transparent Logical. If `TRUE` (default), returns a 4-channel RGBA `image` object
#'   where background pixels have an alpha channel value of `0` (100% transparent).
#' @param bg_color Optional background color to fill when `transparent = FALSE`. Accepts standard
#'   R color names (e.g., `"white"`, `"black"`) or hex strings (e.g., `"#FFFFFF"`). Defaults to `"black"`.
#' @return An `image` object (4-channel RGBA if `transparent = TRUE`, 3-channel RGB otherwise) or a list of images.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # Fast background removal with transparent alpha channel (CPU friendly)
#'   img_trans <- image_remove_bg_dl(img, model = "u2netp")
#'
#'   # Razor-sharp plant boundary cutout (leaf margins, serrations, petioles)
#'   img_ben <- image_remove_bg_dl(img, model = "ben2")
#'
#'   # Robust background removal under complex studio shadows and reflections
#'   img_rmbg <- image_remove_bg_dl(img, model = "rmbg-1.4")
#'
#'   # Replace background with solid white
#'   img_white <- image_remove_bg_dl(img, model = "u2netp", transparent = FALSE, bg_color = "white")
#'
#'   # Extract specific object via natural language prompt (Grounded-SAM)
#'   img_fruit <- image_remove_bg_dl(img, model = "grounded-sam", prompt = "red fruit")
#'
#'   # Extract instances matching a visual exemplar click (PerSAM)
#'   img_exemplar <- image_remove_bg_dl(img, model = "persam", prompt = c(120, 150))
#' }
image_remove_bg_dl <- function(img,
                               model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                               threshold = 0.5,
                               fill_hull = TRUE,
                               filter = 0,
                               erode = 0,
                               dilate = 0,
                               opening = 0,
                               closing = 0,
                               min_area = 0,
                               invert = FALSE,
                               pick_object = FALSE,
                               prompt = NULL,
                               exemplar = FALSE,
                               threads = 0,
                               engine = c("cpu", "gpu"),
                               device_id = -1,
                               verbose = TRUE,
                               transparent = TRUE,
                               bg_color = "black",
                               plot = TRUE,
                               dir = pliman_model_dir(),
                               ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_remove_bg_dl(x, model = model, threshold = threshold,
                         fill_hull = fill_hull, filter = filter,
                         erode = erode, dilate = dilate,
                         opening = opening, closing = closing,
                         min_area = min_area, invert = invert,
                         pick_object = pick_object, prompt = prompt,
                         exemplar = exemplar,
                         threads = threads,
                         engine = engine, device_id = device_id,
                         verbose = verbose,
                         transparent = transparent, bg_color = bg_color,
                         plot = FALSE, dir = dir)
    })
    return(res)
  }

  if (is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  mask <- .compute_dl_mask(mat = mat,
                           model = model,
                           threshold = threshold,
                           fill_hull = fill_hull,
                           filter = filter,
                           erode = erode,
                           dilate = dilate,
                           opening = opening,
                           closing = closing,
                           min_area = min_area,
                           invert = invert,
                           pick_object = pick_object,
                           prompt = prompt,
                           exemplar = exemplar,
                           threads = threads,
                           engine = engine,
                           device_id = device_id,
                           verbose = verbose,
                           dir = dir)

  if (isTRUE(transparent)) {
    # Generate 4-channel RGBA image with alpha layer
    is_raw_img <- is.raw(mat)
    nch <- if (length(dims) >= 3) dims[3] else 1

    if (nch >= 3) {
      R <- mat[, , 1]
      G <- mat[, , 2]
      B <- mat[, , 3]
    } else {
      R <- G <- B <- mat[, , 1]
    }

    if (is_raw_img) {
      alpha_layer <- as.raw(ifelse(mask, 255, 0))
      arr <- array(c(R, G, B, alpha_layer), dim = c(orig_w, orig_h, 4))
    } else {
      alpha_layer <- ifelse(mask, 1.0, 0.0)
      arr <- array(c(R, G, B, alpha_layer), dim = c(orig_w, orig_h, 4))
    }
    out <- as_image(arr, colormode = "Color")
  } else {
    out <- as_image(mat)
    # Fill background with bg_color
    col_rgb <- tryCatch(col2rgb(bg_color) / 255.0, error = function(e) c(0, 0, 0))
    if (length(dims) == 3) {
      for (k in seq_len(dims[3])) {
        fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
        out[, , k][!mask] <- fill_val
      }
    } else {
      fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
      out[!mask] <- fill_val
    }
  }

  if (isTRUE(plot)) {
    plot(out, ...)
  }

  invisible(out)
}

# ==============================================================================
# SECTION 3.5: GROUNDED-SAM (GROUNDING DINO + SAM 2.1) INSTANCE SEGMENTATION
# ==============================================================================

# Global cache for BERT vocabulary
.pliman_vocab_cache <- new.env(parent = emptyenv())

# Native WordPiece tokenizer for BERT / Grounding DINO (Zero external dependencies)
.bert_tokenize <- function(text, vocab_path) {
  if (!file.exists(vocab_path)) {
    stop("BERT vocab file not found: ", vocab_path)
  }
  if (!exists(vocab_path, envir = .pliman_vocab_cache, inherits = FALSE)) {
    vocab_lines <- readLines(vocab_path, warn = FALSE, encoding = "UTF-8")
    vocab_env <- new.env(hash = TRUE, parent = emptyenv(), size = length(vocab_lines))
    for (i in seq_along(vocab_lines)) {
      vocab_env[[vocab_lines[i]]] <- i - 1L
    }
    assign(vocab_path, vocab_env, envir = .pliman_vocab_cache)
  } else {
    vocab_env <- get(vocab_path, envir = .pliman_vocab_cache)
  }

  if (length(text) > 1) {
    text <- paste(trimws(text), collapse = " and ")
  }
  text_clean <- tolower(trimws(text))
  # Support commas, semicolons, and normalize delimiters into natural conjunctions
  text_clean <- gsub("[,;]+", " and ", text_clean)
  # Normalize internal periods between class names (e.g. "people . book ." -> "people and book .")
  text_clean <- gsub("\\s*\\.\\s*(?=[a-zA-Z0-9])", " and ", text_clean, perl = TRUE)
  # Collapse redundant 'and's or spaces
  text_clean <- gsub("\\b(and\\s*)+", "and ", text_clean)
  text_clean <- gsub("\\s+", " ", text_clean)
  text_clean <- trimws(text_clean)
  # Ensure clean trailing period expected by Grounding DINO
  text_clean <- gsub("\\s*\\.*$", " .", text_clean)
  if (!grepl("\\.$", text_clean)) {
    text_clean <- paste0(text_clean, " .")
  }

  # Insert spaces around punctuation
  text_spaced <- gsub("([[:punct:]])", " \\1 ", text_clean)
  words <- strsplit(text_spaced, "\\s+")[[1]]
  words <- words[nzchar(words)]

  token_ids <- 101L # [CLS]
  token_words <- "[CLS]"
  unk_id <- 100L

  for (w in words) {
    start <- 1L
    w_len <- nchar(w)
    cur_ids <- integer(0)
    is_unk <- FALSE

    while (start <= w_len) {
      end <- w_len
      found_sub <- FALSE
      while (start <= end) {
        sub <- substr(w, start, end)
        if (start > 1L) sub <- paste0("##", sub)
        if (exists(sub, envir = vocab_env, inherits = FALSE)) {
          cur_ids <- c(cur_ids, as.integer(get(sub, envir = vocab_env)))
          start <- end + 1L
          found_sub <- TRUE
          break
        }
        end <- end - 1L
      }
      if (!found_sub) {
        is_unk <- TRUE
        break
      }
    }

    if (is_unk) {
      token_ids <- c(token_ids, unk_id)
      token_words <- c(token_words, w)
    } else {
      token_ids <- c(token_ids, cur_ids)
      token_words <- c(token_words, rep(w, length(cur_ids)))
    }
  }

  token_ids <- c(token_ids, 102L) # [SEP]
  token_words <- c(token_words, "[SEP]")

  return(list(
    input_ids = token_ids,
    attention_mask = rep(1L, length(token_ids)),
    token_type_ids = rep(0L, length(token_ids)),
    token_words = token_words
  ))
}

# Helper to generate rounded rectangle polygon vertices
.rounded_rect_pts <- function(xleft, ybottom, xright, ytop, r = NULL, n = 8) {
  x0 <- min(xleft, xright)
  x1 <- max(xleft, xright)
  y0 <- min(ybottom, ytop)
  y1 <- max(ybottom, ytop)

  w <- x1 - x0
  h <- y1 - y0

  if (is.null(r) || is.na(r) || r <= 0) {
    r <- min(w * 0.25, h * 0.35, 4)
  } else {
    r <- min(r, w / 2, h / 2)
  }

  theta <- seq(0, pi / 2, length.out = n)
  arc_tr <- cbind((x1 - r) + r * cos(theta), (y0 + r) - r * sin(theta))
  arc_tl <- cbind((x0 + r) - r * sin(theta), (y0 + r) - r * cos(theta))
  arc_bl <- cbind((x0 + r) - r * cos(theta), (y1 - r) + r * sin(theta))
  arc_br <- cbind((x1 - r) + r * sin(theta), (y1 - r) + r * cos(theta))

  rbind(arc_tr, arc_tl, arc_bl, arc_br)
}

# Helper to draw YOLO-style bounding boxes with class label and confidence score
.plot_yolo_bboxes <- function(boxes, palette_colors, lwd = 2) {
  if (is.null(boxes) || !is.data.frame(boxes) || nrow(boxes) == 0) return(invisible(NULL))

  u <- graphics::par("usr")
  img_h <- abs(u[3] - u[4])
  base_th <- abs(graphics::strheight("Ag", units = "user", cex = 1))
  target_th <- max(7, min(14, img_h * 0.012))
  cex_val <- if (base_th > 0) max(0.45, min(0.65, target_th / base_th)) else 0.55

  num_inst <- nrow(boxes)
  for (i in seq_len(num_inst)) {
    bx <- boxes[i, ]
    k_col <- palette_colors[((i - 1) %% length(palette_colors)) + 1]
    x1 <- bx$xmin; y1 <- bx$ymin; x2 <- bx$xmax; y2 <- bx$ymax

    # Draw box outline
    graphics::rect(xleft = x1, ybottom = y1, xright = x2, ytop = y2, border = k_col, lwd = lwd)

    # Label text: e.g. "person 0.85"
    has_lbl <- !is.null(bx$label) && !is.na(bx$label) && nzchar(as.character(bx$label))
    has_scr <- !is.null(bx$score) && !is.na(bx$score)
    lbl_txt <- if (has_lbl && has_scr) {
      sprintf("%s %.2f", bx$label, bx$score)
    } else if (has_lbl) {
      as.character(bx$label)
    } else if (has_scr) {
      sprintf("%.2f", bx$score)
    } else {
      as.character(bx$id)
    }

    tw <- abs(graphics::strwidth(lbl_txt, units = "user", cex = cex_val, font = 2))
    th <- abs(graphics::strheight(lbl_txt, units = "user", cex = cex_val, font = 2))
    pad_x <- max(3, th * 0.3)
    pad_y <- max(2, th * 0.2)

    badge_w <- tw + 2 * pad_x
    badge_h <- th + 2 * pad_y

    outside <- (y1 - badge_h) >= min(u[3], u[4])
    by0 <- if (outside) y1 - badge_h else y1
    by1 <- if (outside) y1 else y1 + badge_h
    bx0 <- x1
    bx1 <- x1 + badge_w

    # Ensure badge does not overflow right canvas boundary
    x_max_plot <- max(u[1], u[2])
    if (bx1 > x_max_plot) {
      shift <- bx1 - x_max_plot
      bx0 <- max(min(u[1], u[2]), bx0 - shift)
      bx1 <- bx0 + badge_w
    }

    # Text contrast color (black on bright badges, white on dark badges)
    rgb_vals <- tryCatch(grDevices::col2rgb(k_col)[, 1] / 255.0, error = function(e) c(0, 1, 0))
    lum <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
    txt_col <- if (lum > 0.55) "#111111" else "#ffffff"

    # YOLO rounded badge tag
    r_val <- min(badge_h * 0.35, 4)
    pts <- .rounded_rect_pts(bx0, by0, bx1, by1, r = r_val)
    graphics::polygon(pts[, 1], pts[, 2], col = k_col, border = NA)
    graphics::text(bx0 + pad_x, (by0 + by1) / 2, labels = lbl_txt, col = txt_col, font = 2, adj = c(0, 0.5), cex = cex_val)
  }
}

# Helper to draw YOLO pose keypoints and skeleton connections
.plot_yolo_keypoints <- function(keypoints_list,
                                 palette_colors,
                                 kpt_threshold = 0.3,
                                 lwd = 2,
                                 kpt_radius = 4) {
  if (length(keypoints_list) == 0) return(invisible(NULL))

  # Define colors for limbs
  limb_colors <- c(
    rep("#FF4B4B", 4), # facial: red
    rep("#FFA500", 4), # torso: orange
    rep("#00CC66", 2), # left arm: green
    rep("#0099FF", 2), # right arm: blue
    rep("#9933FF", 2), # left leg: purple
    rep("#FF00FF", 2)  # right leg: magenta
  )

  for (i in seq_along(keypoints_list)) {
    kp_df <- keypoints_list[[i]]
    if (is.null(kp_df) || nrow(kp_df) < 17) next
    k_col <- palette_colors[((i - 1) %% length(palette_colors)) + 1]

    # Draw skeleton limbs
    for (b in seq_len(nrow(.coco_skeleton_pairs))) {
      p1 <- .coco_skeleton_pairs[b, 1]
      p2 <- .coco_skeleton_pairs[b, 2]
      if (kp_df$conf[p1] >= kpt_threshold && kp_df$conf[p2] >= kpt_threshold) {
        graphics::lines(
          c(kp_df$x[p1], kp_df$x[p2]),
          c(kp_df$y[p1], kp_df$y[p2]),
          col = limb_colors[min(b, length(limb_colors))],
          lwd = lwd
        )
      }
    }

    # Draw keypoint dots
    valid_kpts <- which(kp_df$conf >= kpt_threshold)
    if (length(valid_kpts) > 0) {
      graphics::points(
        kp_df$x[valid_kpts],
        kp_df$y[valid_kpts],
        pch = 21,
        bg = k_col,
        col = "white",
        cex = max(0.6, kpt_radius / 3),
        lwd = 1.5
      )
    }
  }
  invisible(NULL)
}

# Helper to format human-readable instance detection summary
.format_detection_summary <- function(labels) {
  if (length(labels) == 0) return("0 objects")
  counts <- table(factor(labels, levels = unique(labels)))
  parts <- vapply(names(counts), function(cls) {
    cnt <- counts[[cls]]
    cls_lower <- tolower(cls)
    if (cnt == 1) {
      paste0(cnt, " ", cls)
    } else {
      # Inherent or already plural nouns
      if (cls_lower %in% c("shorts", "pants", "jeans", "trousers", "glasses", "sunglasses", "scissors", "children", "people")) {
        paste0(cnt, " ", cls)
      } else if (cls_lower == "child") {
        paste0(cnt, " children")
      } else if (cls_lower == "leaf") {
        paste0(cnt, " leaves")
      } else if (cls_lower == "person") {
        paste0(cnt, " people")
      } else if (grepl("[^aeiou]y$", cls_lower)) {
        # Words ending in consonant + y (e.g., butterfly -> butterflies)
        stem <- substr(cls, 1, nchar(cls) - 1)
        paste0(cnt, " ", stem, "ies")
      } else if (grepl("(ss|sh|ch|x|z)$", cls_lower)) {
        # Words ending in -ss, -sh, -ch, -x, -z (e.g. dress -> dresses, sunglass -> sunglasses)
        paste0(cnt, " ", cls, "es")
      } else if (grepl("s$", cls_lower)) {
        # Words ending in single s (already plural, e.g. cars, shoes, cows)
        paste0(cnt, " ", cls)
      } else {
        paste0(cnt, " ", cls, "s")
      }
    }
  }, character(1), USE.NAMES = FALSE)

  if (length(parts) == 1) {
    parts[1]
  } else if (length(parts) == 2) {
    paste(parts[1], "and", parts[2])
  } else {
    paste0(paste(parts[-length(parts)], collapse = ", "), ", and ", parts[length(parts)])
  }
}

# Helper to print highlighted detection summary block
.print_detection_summary <- function(summary_str, title = "Detection Summary") {
  if (is.null(summary_str)) return(invisible(NULL))

  if (is.data.frame(summary_str)) {
    if (nrow(summary_str) == 0) return(invisible(NULL))
    lbl_col <- intersect(c("class", "label"), tolower(names(summary_str)))
    cnt_col <- intersect(c("count", "freq"), tolower(names(summary_str)))
    if (length(lbl_col) > 0 && length(cnt_col) > 0) {
      lbls <- rep(as.character(summary_str[[lbl_col[1]]]), as.integer(summary_str[[cnt_col[1]]]))
      summary_text <- .format_detection_summary(lbls)
    } else {
      summary_text <- paste(apply(summary_str, 1, paste, collapse = ": "), collapse = ", ")
    }
  } else if (is.table(summary_str)) {
    if (length(summary_str) == 0) return(invisible(NULL))
    lbls <- rep(names(summary_str), as.integer(summary_str))
    summary_text <- .format_detection_summary(lbls)
  } else if (is.character(summary_str)) {
    if (length(summary_str) == 0 || !any(nzchar(summary_str))) return(invisible(NULL))
    if (length(summary_str) > 1) {
      summary_text <- .format_detection_summary(summary_str)
    } else {
      summary_text <- summary_str
    }
  } else {
    summary_text <- as.character(summary_str)
  }

  if (is.null(summary_text) || !nzchar(summary_text) || summary_text == "0 objects") return(invisible(NULL))

  cat("\n")
  cli::cli_rule(left = paste0("{.bold ", title, "}"))
  cli::cli_alert_success("pliman detected {.bold {summary_text}} in the image.")
  cli::cli_rule()
  cat("\n")
}



# Grounded-SAM runner: detects instances from text prompt and segments each with SAM 2.1
.run_grounded_sam <- function(mat,
                              prompt,
                              threshold = 0.5,
                              box_threshold = 0.25,
                              text_threshold = 0.25,
                              iou_threshold = 0.5,
                              threads = 0,
                              engine = c("cpu", "gpu"),
                              device_id = -1,
                              fill_hull = TRUE,
                              filter = 0,
                              erode = 0,
                              dilate = 0,
                              opening = 0,
                              closing = 0,
                              min_area = 0,
                              invert = FALSE,
                              mask = TRUE,
                              verbose = TRUE,
                              dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  dir <- pliman_model_dir(dir)
  lib_file <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_file) || !file.exists(lib_file)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_file <- onnx_install(engine = engine)
  }

  # Ensure Grounding DINO and SAM 2.1 (if masking) are downloaded
  pliman_download_model("grounded-sam", dir = dir)
  if (isTRUE(mask)) {
    pliman_download_model("sam2.1", dir = dir)
  }

  dino_model <- file.path(dir, "groundingdino-tiny.onnx")
  vocab_file <- file.path(dir, "vocab.txt")
  sam_enc <- file.path(dir, "sam2.1.encoder.onnx")
  sam_dec <- file.path(dir, "sam2.1.decoder.onnx")

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # 1. Preprocess for Grounding DINO (800x800, ImageNet norm)
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Tokenizing prompt and preparing image...",
      msg_done = "Prompt tokenized"
    )
  }

  dino_tensor <- .preprocess_nchw(mat, target_size = 800L,
                                  mean = c(0.485, 0.456, 0.406),
                                  std = c(0.229, 0.224, 0.225),
                                  letterbox = FALSE)

  # 2. Tokenize prompt
  tok <- .bert_tokenize(prompt, vocab_file)

  # 3. Run Grounding DINO
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Detecting objects with Grounding DINO [{toupper(engine)}]...",
      msg_done = "Grounding DINO object detection complete"
    )
  }

  dino_res <- run_grounding_dino_cpp(
    pixel_values = as.numeric(dino_tensor),
    input_ids = as.integer(tok$input_ids),
    token_type_ids = as.integer(tok$token_type_ids),
    attention_mask = as.integer(tok$attention_mask),
    model_path = normalizePath(dino_model, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
    box_threshold = box_threshold,
    text_threshold = text_threshold,
    iou_threshold = iou_threshold,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  num_boxes <- length(dino_res$scores)
  if (num_boxes == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
      cli::cli_alert_warning("Grounding DINO detected 0 instances for prompt: {.val {prompt}}.")
    }
    return(list(
      boxes = data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                         xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                         label = character(0), stringsAsFactors = FALSE),
      labels = as_image(matrix(0L, nrow = orig_w, ncol = orig_h), colormode = "Grayscale", storage = "integer"),
      mask = as_image(matrix(FALSE, nrow = orig_w, ncol = orig_h))
    ))
  }

  # Scale normalized box coordinates to pixel coordinates on original image
  raw_boxes <- dino_res$boxes
  boxes_px <- matrix(0.0, nrow = num_boxes, ncol = 4)
  boxes_px[, 1] <- pmax(1.0, raw_boxes[, 1] * orig_w)
  boxes_px[, 2] <- pmax(1.0, raw_boxes[, 2] * orig_h)
  boxes_px[, 3] <- pmin(as.double(orig_w), raw_boxes[, 3] * orig_w)
  boxes_px[, 4] <- pmin(as.double(orig_h), raw_boxes[, 4] * orig_h)

  # Map token indices to labels (avoiding conjunctions or punctuation)
  labels_str <- character(num_boxes)
  stop_tokens <- c("[cls]", "[sep]", "and", ".", ",", ";")
  for (i in seq_len(num_boxes)) {
    t_idx <- dino_res$token_indices[i] + 1L
    w_cand <- if (t_idx >= 1 && t_idx <= length(tok$token_words)) tok$token_words[t_idx] else ""
    if (tolower(w_cand) %in% stop_tokens || !nzchar(w_cand)) {
      dists <- abs(seq_along(tok$token_words) - t_idx)
      valid_mask <- !(tolower(tok$token_words) %in% stop_tokens)
      if (any(valid_mask)) {
        dists[!valid_mask] <- Inf
        w_cand <- tok$token_words[which.min(dists)]
      } else {
        w_cand <- prompt
      }
    }
    labels_str[i] <- w_cand
  }

  boxes_df <- data.frame(
    id = seq_len(num_boxes),
    xmin = round(boxes_px[, 1], 1),
    ymin = round(boxes_px[, 2], 1),
    xmax = round(boxes_px[, 3], 1),
    ymax = round(boxes_px[, 4], 1),
    score = round(dino_res$scores, 3),
    label = labels_str,
    stringsAsFactors = FALSE
  )

  counts_df <- as.data.frame(table(factor(labels_str, levels = unique(labels_str))), stringsAsFactors = FALSE)
  colnames(counts_df) <- c("label", "count")
  summary_str <- .format_detection_summary(labels_str)

  # If no boxes found, return early
  if (nrow(boxes_df) == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = counts_df,
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = NULL
    ))
  }

  # 4. If mask = FALSE, skip SAM 2.1 completely (pure object detection mode)
  if (!isTRUE(mask)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = counts_df,
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = matrix(FALSE, nrow = orig_w, ncol = orig_h),
      features = data.frame()
    ))
  }

  # 4. Run SAM 2.1 on all detected boxes
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Segmenting {num_boxes} instance{?s} with SAM 2.1 [{toupper(engine)}]...",
      msg_done = "SAM 2.1 segmented {num_boxes} instance{?s}"
    )
  }

  sam_tensor <- .preprocess_nchw(mat, target_size = 1024L,
                                mean = c(0.485, 0.456, 0.406),
                                std = c(0.229, 0.224, 0.225),
                                letterbox = FALSE)

  raw_masks <- run_sam2_instances_cpp(
    tensor_vec = as.numeric(sam_tensor),
    boxes = boxes_px,
    orig_w = as.double(orig_w),
    orig_h = as.double(orig_h),
    encoder_path = normalizePath(sam_enc, winslash = "/", mustWork = FALSE),
    decoder_path = normalizePath(sam_dec, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  # 5. Post-process masks directly into a single multi-label matrix
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing multi-label instance masks...",
      msg_done = "Instance masks ready"
    )
  }

  combined_labels <- matrix(0L, nrow = orig_w, ncol = orig_h)

  for (k in seq_len(num_boxes)) {
    prob_crop <- raw_masks[[k]]
    # Resize continuous float logits to full resolution for silky-smooth sub-pixel boundaries
    mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
    m <- mask_prob >= threshold

    bx <- boxes_px[k, ]
    x1 <- max(1L, floor(bx[1]) - 2L)
    x2 <- min(orig_w, ceiling(bx[3]) + 2L)
    y1 <- max(1L, floor(bx[2]) - 2L)
    y2 <- min(orig_h, ceiling(bx[4]) + 2L)

    if (isTRUE(fill_hull)) {
      sub_m <- m[x1:x2, y1:y2]
      m[x1:x2, y1:y2] <- (fill_holes_cpp(sub_m) > 0)
    }
    if (filter > 0) {
      m_arr <- as.array(m)
      if (length(dim(m_arr)) == 2) dim(m_arr) <- c(orig_w, orig_h, 1L)
      f_res <- median_filter_binary_cpp(m_arr, orig_w, orig_h, 1L, as.integer(filter))
      m <- (f_res[, , 1] > 0)
    }
    if (erode > 0) {
      m <- (erode_cpp(m, raio = as.integer(erode)) > 0)
    }
    if (dilate > 0) {
      m <- (dilate_cpp(m, raio = as.integer(dilate)) > 0)
    }
    if (opening > 0) {
      m <- (erode_cpp(m, raio = as.integer(opening)) > 0)
      m <- (dilate_cpp(m, raio = as.integer(opening)) > 0)
    }
    if (closing > 0) {
      m <- (dilate_cpp(m, raio = as.integer(closing)) > 0)
      m <- (erode_cpp(m, raio = as.integer(closing)) > 0)
    }
    if (min_area > 0) {
      lbls <- bwlabel_cpp(m)
      if (max(lbls) > 0) {
        tbl <- table(lbls[lbls > 0])
        keep_ids <- as.integer(names(tbl[tbl >= min_area]))
        m <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
      }
    }

    combined_labels[m] <- as.integer(k)
  }

  labels_img <- as_image(combined_labels, colormode = "Grayscale", storage = "integer")

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }
  contornos <- contour(labels_img)
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  return(list(
    boxes = boxes_df,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (combined_labels > 0L),
    features = poly_measures(contornos)
  ))
}

# PerSAM (Personalize Segment Anything) runner: segments all instances similar to visual exemplar(s)
.run_persam <- function(mat,
                        exemplar_points = NULL,
                        sim_threshold = 0.5,
                        min_dist = 16,
                        iou_threshold = 0.5,
                        max_objects = 300,
                        feat_res = 256,
                        threshold = 0.5,
                        threads = 0,
                        engine = c("cpu", "gpu"),
                        device_id = -1,
                        fill_hull = TRUE,
                        filter = 0,
                        erode = 0,
                        dilate = 0,
                        opening = 0,
                        closing = 0,
                        min_area = 0,
                        invert = FALSE,
                        mask = TRUE,
                        verbose = TRUE,
                        dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  if (is.character(feat_res)) {
    feat_res <- switch(tolower(feat_res[1]),
      "fast" = 64L,
      "vit" = 64L,
      "low" = 64L,
      "medium" = 256L,
      "high" = 1024L,
      "dense" = 1024L,
      "conv" = 1024L,
      "grain" = 1024L,
      "grains" = 1024L,
      "ultra" = 1024L,
      1024L
    )
  }
  feat_res <- as.integer(feat_res)
  if (is.na(feat_res) || feat_res <= 64L) {
    feat_res <- 64L
  } else if (feat_res <= 256L) {
    feat_res <- 256L
  } else {
    feat_res <- 1024L
  }
  use_gpu <- (engine == "gpu")
  dir <- pliman_model_dir(dir)
  lib_file <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_file) || !file.exists(lib_file)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_file <- onnx_install(engine = engine)
  }

  # Ensure SAM 2.1 model is downloaded
  pliman_download_model("sam2.1", dir = dir)
  sam_enc <- file.path(dir, "sam2.1.encoder.onnx")
  sam_dec <- file.path(dir, "sam2.1.decoder.onnx")

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Interactive picking of exemplar points if not supplied
  if (is.null(exemplar_points)) {
    cli::cli_alert_info("Click on 1 or more exemplar object(s) in the plot window. Press <Esc> or right-click when finished.")
    plot(as_image(mat))
    pts <- tryCatch(graphics::locator(n = 512, type = "p", col = "cyan", pch = 19), error = function(e) NULL)
    if (!is.null(pts) && length(pts$x) > 0) {
      ex_x <- pts$x
      ex_y <- pts$y
    } else {
      cli::cli_abort("No point selected.")
    }
  } else if (is.matrix(exemplar_points) || is.data.frame(exemplar_points)) {
    ex_x <- as.numeric(exemplar_points[, 1])
    ex_y <- as.numeric(exemplar_points[, 2])
    if (length(ex_x) == 0) {
      cli::cli_abort("No exemplar point provided.")
    }
  } else if (is.numeric(exemplar_points)) {
    if (length(exemplar_points) == 2) {
      ex_x <- exemplar_points[1]
      ex_y <- exemplar_points[2]
    } else if (length(exemplar_points) > 0 && length(exemplar_points) %% 2 == 0) {
      half <- length(exemplar_points) / 2
      ex_x <- exemplar_points[seq_len(half)]
      ex_y <- exemplar_points[(half + 1L):length(exemplar_points)]
    } else {
      cli::cli_abort("Invalid {.arg exemplar_points}. Expected numeric coordinates or a matrix/data.frame.")
    }
  } else {
    cli::cli_abort("Invalid {.arg exemplar_points}. Expected numeric coordinates or a matrix/data.frame.")
  }

  # 1. Preprocess for SAM 2.1 (1024x1024, ImageNet norm)
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Computing visual embeddings with SAM 2.1 [{toupper(engine)}]...",
      msg_done = "SAM 2.1 visual embeddings ready"
    )
  }

  sam_tensor <- .preprocess_nchw(mat, target_size = 1024L,
                                mean = c(0.485, 0.456, 0.406),
                                std = c(0.229, 0.224, 0.225),
                                letterbox = FALSE)

  # 2. Run PerSAM C++ routine
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Calculating PerSAM feature similarity map & segmenting matches...",
      msg_done = "PerSAM exemplar matching complete"
    )
  }

  persam_res <- run_sam2_persam_cpp(
    tensor_vec = as.numeric(sam_tensor),
    exemplar_x = as.numeric(ex_x),
    exemplar_y = as.numeric(ex_y),
    orig_w = as.double(orig_w),
    orig_h = as.double(orig_h),
    encoder_path = normalizePath(sam_enc, winslash = "/", mustWork = FALSE),
    decoder_path = normalizePath(sam_dec, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
    sim_threshold = as.double(sim_threshold),
    min_dist = as.double(min_dist),
    iou_threshold = as.double(iou_threshold),
    max_objects = as.integer(max_objects),
    feat_res = as.integer(feat_res),
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  # Process similarity map: normalize to [0, 1] as an image object (resolution governed by feat_res)
  sim_raw <- persam_res$similarity_map
  sim_min <- min(sim_raw)
  sim_max <- max(sim_raw)
  sim_norm <- if (sim_max > sim_min) (sim_raw - sim_min) / (sim_max - sim_min) else sim_raw
  sim_img <- as_image(sim_norm, colormode = "Grayscale", storage = "double")
  attr(sim_img, "raw_matrix") <- sim_raw
  attr(sim_img, "raw_range") <- c(sim_min, sim_max)

  num_inst <- length(persam_res$scores)
  if (num_inst == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
      cli::cli_alert_warning("PerSAM found 0 instances matching the exemplar (sim_threshold = {sim_threshold}).")
    }
    return(list(
      boxes = data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                         xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                         label = character(0), stringsAsFactors = FALSE),
      counts = data.frame(label = "exemplar", count = 0L, stringsAsFactors = FALSE),
      summary = "0 objects similar to the exemplar",
      contours = list(),
      labels = as_image(matrix(0L, nrow = orig_w, ncol = orig_h), colormode = "Grayscale", storage = "integer"),
      mask = as_image(matrix(FALSE, nrow = orig_w, ncol = orig_h)),
      similarity_map = sim_img
    ))
  }

  raw_boxes <- persam_res$boxes
  boxes_df <- data.frame(
    id = seq_len(num_inst),
    xmin = round(raw_boxes[, 1], 1),
    ymin = round(raw_boxes[, 2], 1),
    xmax = round(raw_boxes[, 3], 1),
    ymax = round(raw_boxes[, 4], 1),
    score = round(persam_res$scores, 3),
    label = "exemplar",
    stringsAsFactors = FALSE
  )

  summary_str <- paste0(num_inst, if (num_inst == 1) " object similar to the exemplar" else " objects similar to the exemplar")

  # If mask = FALSE, skip mask post-processing and contour extraction
  if (!isTRUE(mask)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = data.frame(label = "exemplar", count = num_inst, stringsAsFactors = FALSE),
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = matrix(FALSE, nrow = orig_w, ncol = orig_h),
      similarity_map = sim_img,
      features = data.frame()
    ))
  }

  # 3. Post-process masks into multi-label matrix
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing multi-label instance masks...",
      msg_done = "Instance masks ready"
    )
  }

  combined_labels <- matrix(0L, nrow = orig_w, ncol = orig_h)
  for (k in seq_len(num_inst)) {
    prob_crop <- persam_res$masks[[k]]
    mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
    m <- mask_prob >= threshold

    bx <- raw_boxes[k, ]
    x1 <- max(1L, floor(bx[1]) - 2L)
    x2 <- min(orig_w, ceiling(bx[3]) + 2L)
    y1 <- max(1L, floor(bx[2]) - 2L)
    y2 <- min(orig_h, ceiling(bx[4]) + 2L)

    if (isTRUE(fill_hull)) {
      sub_m <- m[x1:x2, y1:y2]
      m[x1:x2, y1:y2] <- (fill_holes_cpp(sub_m) > 0)
    }
    if (filter > 0) {
      m_arr <- as.array(m)
      if (length(dim(m_arr)) == 2) dim(m_arr) <- c(orig_w, orig_h, 1L)
      f_res <- median_filter_binary_cpp(m_arr, orig_w, orig_h, 1L, as.integer(filter))
      m <- (f_res[, , 1] > 0)
    }
    if (erode > 0) {
      m <- (erode_cpp(m, raio = as.integer(erode)) > 0)
    }
    if (dilate > 0) {
      m <- (dilate_cpp(m, raio = as.integer(dilate)) > 0)
    }
    if (opening > 0) {
      m <- (erode_cpp(m, raio = as.integer(opening)) > 0)
      m <- (dilate_cpp(m, raio = as.integer(opening)) > 0)
    }
    if (closing > 0) {
      m <- (dilate_cpp(m, raio = as.integer(closing)) > 0)
      m <- (erode_cpp(m, raio = as.integer(closing)) > 0)
    }
    if (min_area > 0) {
      lbls <- bwlabel_cpp(m)
      if (max(lbls) > 0) {
        tbl <- table(lbls[lbls > 0])
        keep_ids <- as.integer(names(tbl[tbl >= min_area]))
        m <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
      }
    }

    combined_labels[m] <- as.integer(k)
  }

  labels_img <- as_image(combined_labels, colormode = "Grayscale", storage = "integer")

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }

  summary_str <- paste0(num_inst, if (num_inst == 1) " object similar to the exemplar" else " objects similar to the exemplar")
  contornos <- contour(labels_img)
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  return(list(
    boxes = boxes_df,
    counts = data.frame(label = "exemplar", count = num_inst, stringsAsFactors = FALSE),
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (combined_labels > 0L),
    similarity_map = sim_img,
    features = poly_measures(contornos)
  ))
}

# Internal helper to preprocess an image for YOLO models with letterboxing to 640x640
.preprocess_yolo <- function(mat, target_size = 640L) {
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  dim(R) <- dim(G) <- dim(B) <- c(orig_w, orig_h)

  gain <- min(target_size / orig_w, target_size / orig_h)
  new_w <- max(1L, round(orig_w * gain))
  new_h <- max(1L, round(orig_h * gain))

  pad_x <- (target_size - new_w) / 2.0
  pad_y <- (target_size - new_h) / 2.0
  x1 <- floor(pad_x) + 1L
  y1 <- floor(pad_y) + 1L
  x2 <- x1 + new_w - 1L
  y2 <- y1 + new_h - 1L

  R_res <- .bilinear_resize_2d(R, new_w, new_h)
  G_res <- .bilinear_resize_2d(G, new_w, new_h)
  B_res <- .bilinear_resize_2d(B, new_w, new_h)

  fill_val <- 114.0 / 255.0
  canvas_R <- matrix(fill_val, nrow = target_size, ncol = target_size)
  canvas_G <- matrix(fill_val, nrow = target_size, ncol = target_size)
  canvas_B <- matrix(fill_val, nrow = target_size, ncol = target_size)

  canvas_R[x1:x2, y1:y2] <- R_res
  canvas_G[x1:x2, y1:y2] <- G_res
  canvas_B[x1:x2, y1:y2] <- B_res

  vec <- c(as.numeric(canvas_R), as.numeric(canvas_G), as.numeric(canvas_B))
  attr(vec, "gain") <- gain
  attr(vec, "pad_x") <- pad_x
  attr(vec, "pad_y") <- pad_y
  vec
}

# Internal helper to run YOLO detection, segmentation, or pose estimation
.run_yolo <- function(mat,
                      model = "yolo26n",
                      conf_threshold = 0.25,
                      iou_threshold = 0.45,
                      labels = NULL,
                      threads = 0,
                      engine = c("cpu", "gpu"),
                      device_id = -1,
                      fill_hull = TRUE,
                      filter = 0,
                      erode = 0,
                      dilate = 0,
                      opening = 0,
                      closing = 0,
                      min_area = 0,
                      dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (file.exists(model)) {
    model_file <- normalizePath(model, winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model, dir = dir)
  }

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  tensor <- .preprocess_yolo(mat, target_size = 640L)

  raw_res <- run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = orig_w,
    orig_h = orig_h,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  class_names <- if (!is.null(labels)) labels else .coco_classes

  num_boxes <- nrow(raw_res$boxes)
  if (num_boxes > 0) {
    c_ids <- raw_res$class_ids
    lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))
    df_boxes <- data.frame(
      id = seq_len(num_boxes),
      xmin = raw_res$boxes[, 1],
      ymin = raw_res$boxes[, 2],
      xmax = raw_res$boxes[, 3],
      ymax = raw_res$boxes[, 4],
      label = lbls,
      score = round(raw_res$scores, 4),
      class_id = c_ids,
      stringsAsFactors = FALSE
    )
    counts <- table(factor(df_boxes$label, levels = unique(df_boxes$label)))
    counts_df <- data.frame(
      class = names(counts),
      count = as.integer(counts),
      stringsAsFactors = FALSE
    )
    summary_str <- .format_detection_summary(df_boxes$label)
  } else {
    df_boxes <- data.frame(
      id = integer(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      label = character(0),
      score = numeric(0),
      class_id = integer(0),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(class = character(0), count = integer(0), stringsAsFactors = FALSE)
    summary_str <- "0 objects"
  }

  lbl_mat <- raw_res$labels
  if (num_boxes > 0 && any(lbl_mat > 0L)) {
    for (k in seq_len(num_boxes)) {
      bx <- df_boxes[k, ]
      x1 <- max(1L, floor(bx$xmin) - 2L)
      x2 <- min(orig_w, ceiling(bx$xmax) + 2L)
      y1 <- max(1L, floor(bx$ymin) - 2L)
      y2 <- min(orig_h, ceiling(bx$ymax) + 2L)
      sub_m <- (lbl_mat[x1:x2, y1:y2] == k)
      if (!any(sub_m)) next

      if (isTRUE(fill_hull)) {
        sub_m <- (fill_holes_cpp(sub_m) > 0)
      }
      if (filter > 0) {
        m_arr <- as.array(sub_m)
        if (length(dim(m_arr)) == 2) dim(m_arr) <- c(nrow(sub_m), ncol(sub_m), 1L)
        sub_m <- (median_filter_binary_cpp(m_arr, nrow(sub_m), ncol(sub_m), 1L, as.integer(filter))[, , 1] > 0)
      }
      if (erode > 0) {
        sub_m <- (erode_cpp(sub_m, raio = as.integer(erode)) > 0)
      }
      if (dilate > 0) {
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(dilate)) > 0)
      }
      if (opening > 0) {
        sub_m <- (erode_cpp(sub_m, raio = as.integer(opening)) > 0)
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(opening)) > 0)
      }
      if (closing > 0) {
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(closing)) > 0)
        sub_m <- (erode_cpp(sub_m, raio = as.integer(closing)) > 0)
      }
      if (min_area > 0) {
        lbls <- bwlabel_cpp(sub_m)
        if (max(lbls) > 0) {
          tbl <- table(lbls[lbls > 0])
          keep_ids <- as.integer(names(tbl[tbl >= min_area]))
          sub_m <- matrix(lbls %in% keep_ids, nrow = nrow(sub_m), ncol = ncol(sub_m))
        }
      }
      lbl_mat[x1:x2, y1:y2][lbl_mat[x1:x2, y1:y2] == k] <- 0L
      lbl_mat[x1:x2, y1:y2][sub_m] <- as.integer(k)
    }
  }

  labels_img <- as_image(lbl_mat, colormode = "Grayscale", storage = "integer")
  contornos <- if (num_boxes > 0 && any(lbl_mat > 0L)) contour(labels_img) else list()
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  feats <- if (length(contornos) > 0) {
    tryCatch(poly_measures(contornos), error = function(e) data.frame())
  } else {
    data.frame()
  }

  has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && num_boxes > 0
  kpts_list <- if (has_kpts) {
    lapply(seq_len(num_boxes), function(k) {
      vals <- raw_res$keypoints[k, ]
      data.frame(
        id = k,
        keypoint = .coco_keypoints,
        x = vals[seq(1, 51, by = 3)],
        y = vals[seq(2, 51, by = 3)],
        conf = round(vals[seq(3, 51, by = 3)], 4),
        stringsAsFactors = FALSE
      )
    })
  } else {
    list()
  }

  list(
    boxes = df_boxes,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (lbl_mat > 0L),
    features = feats,
    keypoints = kpts_list,
    raw_keypoints = if (has_kpts) raw_res$keypoints else NULL
  )
}

# Internal helper to run StarDist polygon detection
.run_stardist <- function(mat,
                          model = "stardist",
                          prob_threshold = 0.5,
                          nms_threshold = 0.3,
                          threads = 0,
                          engine = c("cpu", "gpu"),
                          device_id = -1,
                          dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (file.exists(model)) {
    model_file <- normalizePath(model, winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model, dir = dir)
  }

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  target_size <- 256L
  tensor <- .preprocess_nchw(
    mat,
    target_size = target_size,
    mean = c(0, 0, 0),
    std = c(1, 1, 1),
    letterbox = FALSE
  )

  raw_res <- run_stardist_cpp(
    tensor_vec = tensor,
    in_w = target_size,
    in_h = target_size,
    orig_w = orig_w,
    orig_h = orig_h,
    prob_threshold = prob_threshold,
    nms_threshold = nms_threshold,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  num_objs <- nrow(raw_res$boxes)
  if (num_objs > 0) {
    df_boxes <- data.frame(
      id = seq_len(num_objs),
      xmin = raw_res$boxes[, 1],
      ymin = raw_res$boxes[, 2],
      xmax = raw_res$boxes[, 3],
      ymax = raw_res$boxes[, 4],
      center_x = raw_res$centers_x,
      center_y = raw_res$centers_y,
      label = rep("object", num_objs),
      score = round(raw_res$scores, 4),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(
      class = "object",
      count = num_objs,
      stringsAsFactors = FALSE
    )
    summary_str <- .format_detection_summary(df_boxes$label)
  } else {
    df_boxes <- data.frame(
      id = integer(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      center_x = numeric(0),
      center_y = numeric(0),
      label = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(class = character(0), count = integer(0), stringsAsFactors = FALSE)
    summary_str <- "0 objects"
  }

  lbl_mat <- raw_res$labels
  labels_img <- as_image(lbl_mat, colormode = "Grayscale", storage = "integer")
  contornos <- if (num_objs > 0 && any(lbl_mat > 0L)) contour(labels_img) else list()
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  feats <- if (length(contornos) > 0) {
    tryCatch(poly_measures(contornos), error = function(e) data.frame())
  } else {
    data.frame()
  }

  list(
    boxes = df_boxes,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (lbl_mat > 0L),
    features = feats,
    polygons_x = raw_res$polygons_x,
    polygons_y = raw_res$polygons_y,
    centers_x = raw_res$centers_x,
    centers_y = raw_res$centers_y
  )
}

#' Deep Learning Semantic & Instance Segmentation, Object Detection, and Visual Exemplar Counting
#'
#' Segments foreground objects (such as leaves, fruits, seeds, grains, roots,
#' animals, or plants) from complex, textured, or out-of-focus backgrounds using pre-trained
#' Deep Learning models in **ONNX** format.
#'
#' This unified function supports five major computer vision workflows:
#' 1. **Salient Foreground Segmentation & Background Removal:** Semantic cutout using architectures
#'    such as U2-Net, BEN2, BRIA RMBG 1.4/2.0, IS-Net, withoutBG, BiRefNet, and Silueta.
#' 2. **Real-Time YOLO26 Instance Segmentation (`yolo26n-seg`):** Real-time instance segmentation
#'    combining bounding box detection with 32 prototype mask channels and continuous bilinear sub-pixel
#'    interpolation. Returns a single multi-label mask image (`type = "mask"`), smooth polygon overlays
#'    (`type = "highlight"`), or background cutouts (`type = "segment"`), accompanied by per-instance
#'    morphological measurements (area, perimeter, circularity, etc.).
#' 3. **Star-Convex Object Detection & Segmentation (`stardist`):** Predicts 48 radial star-convex
#'    distance rays and object probability, specifically tailored for densely packed, touching, or convex
#'    structures such as seeds, grains, nuclei, and cells.
#' 4. **Zero-Shot Text-Prompted Detection & Instance Segmentation (Grounded-SAM):** Combines Grounding DINO
#'    with SAM 2.1 to locate, outline, and count arbitrary concepts described in natural language (e.g., `prompt = "leaf, fruit, insect"`).
#' 5. **One-Shot Visual Exemplar Segmentation & Counting (PerSAM):** Personalize Segment Anything via SAM 2.1.
#'    The user clicks on 1 or more exemplar objects (or provides coordinates), and the algorithm automatically searches,
#'    segments, and counts all visually similar objects across the entire image using Multi-Prototype Max-Pooling.
#'
#' @inheritParams image_binary_dl
#' @param type Character specifying the visual and output modality:
#'   * `"segment"` (default): Returns an `image` object where detected foreground is preserved
#'     and background is replaced by a solid color specified by `col_background`.
#'   * `"mask"`: Returns the binary mask (for semantic models) or single multi-label integer instance map
#'     (for YOLO-seg, StarDist, Grounded-SAM, and PerSAM) as an `image` object (`colormode = "Grayscale"`).
#'   * `"highlight"`: Overlays the original image with semi-transparent colored polygons covering
#'     each detected instance, ideal for visual inspection, quality control, and publication figures.
#'     When `mask = FALSE`, overlays bounding boxes directly.
#'   * `"boxes"`: Overlays bounding boxes on the original image and returns a data frame of detected boxes.
#' @param model Model architecture or path to a custom ONNX file. Built-in choices:
#'   * `"yolo26n-seg"` (or `"yolo26s-seg"`, `"yolo26m-seg"`, `"yolo26l-seg"`, `"yolo26x-seg"`): YOLO26 Instance Segmentation (real-time polygon mask & box detection, 80 COCO categories).
#'   * `"stardist"` (or `"stardist-dsb2018"`): StarDist star-convex radial polygon segmentation for overlapping grains/cells.
#'   * `"grounded-sam"`: Zero-shot open-vocabulary instance segmentation prompted by text.
#'   * `"persam"`: One-shot visual exemplar segmentation guided by reference click points.
#'   * `"u2netp"`, `"u2net"`: U2-Net salient foreground cutout (fast CPU-friendly default).
#'   * `"ben2"`: Boundary-aware Enhanced Network (BEN2) for razor-sharp leaf margins, serrations, and roots.
#'   * `"rmbg-1.4"`, `"rmbg-2.0"`: BRIA state-of-the-art commercial-grade background removal.
#'   * `"withoutbg"`, `"birefnet-lite"`, `"isnet-general-use"`, `"silueta"`: Specialized salient cutout backends.
#' @param col_background Character string or hex code specifying the solid background color when
#'   `type = "segment"`. Defaults to `"white"`. Accepts any valid R color name (e.g., `"black"`, `"transparent"`)
#'   or hexadecimal code (e.g., `"#FFFFFF"`).
#' @param col_highlight Fill color for the semi-transparent overlay polygon when `type = "highlight"`.
#'   Defaults to `"salmon"`. When multiple instances are detected, distinct categorical
#'   colors are automatically generated from a rainbow palette.
#' @param alpha Numeric transparency value in `[0, 1]` for the highlight polygon overlay when `type = "highlight"`.
#'   Default is `0.4` (40% opacity).
#' @param border Color for the border stroke of the highlight polygon when `type = "highlight"`.
#'   Default is `NA` (no border stroke).
#' @param lwd Numeric line width for the bounding boxes and highlight polygon contours. Default is `1`.
#' @param bbox Logical. If `TRUE`, draws YOLO-style bounding boxes with class labels and confidence scores
#'   around detected instances. Defaults to `FALSE`.
#' @param mask Logical. If `TRUE` (default), computes pixel-accurate instance or semantic masks.
#'   When `mask = FALSE` and `bbox = TRUE`, skips the mask decoder step and returns only the detected
#'   bounding box coordinates (fast object detection). Defaults to `TRUE`.
#' @param exemplar Logical. If `TRUE` (or if `model = "persam"`), activates one-shot visual exemplar segmentation (PerSAM)
#'   using SAM 2.1 embeddings. The user can click on 1 or more representative exemplar objects in the plot window
#'   (or pass coordinates via `prompt`), and the model will automatically compute visual feature embeddings,
#'   search the entire image for all similar instances via cosine similarity, and segment all matching objects.
#'   **Multi-Prototype Max-Pooling:** When multiple exemplars are clicked (e.g., one bright log, one dark log, one with moss),
#'   each click forms an independent 256-D prototype vector. The similarity map calculates the maximum cosine similarity
#'   against *any* prototype, preventing vector dilution and capturing wide visual diversity without prototype cancellation. Default is `FALSE`.
#' @param sim_threshold Numeric cosine similarity threshold in `[0, 1]` for PerSAM visual exemplar matching.
#'   Higher values (e.g. `0.55` - `0.65`) require closer visual resemblance, while lower values (e.g. `0.40` - `0.50`)
#'   capture more diverse, shadowed, or smaller instances. Default is `0.5`.
#' @param min_dist Numeric minimum spatial distance in pixels (within the 1024×1024 embedding space) between candidate
#'   peaks to prevent duplicate instance detections of the same object. Default is `16`.
#' @param box_threshold Numeric confidence threshold in `[0, 1]` for Grounded-SAM and YOLO detection bounding
#'   box proposals. Default is `0.25`.
#' @param text_threshold Numeric token logit confidence threshold in `[0, 1]` for Grounded-SAM text-to-box
#'   association. Default is `0.25`.
#' @param iou_threshold Numeric Non-Maximum Suppression (NMS) Intersection over Union (IoU) threshold in `[0, 1]`
#'   for deduplicating overlapping bounding box proposals. Default is `0.5`.
#' @param feat_res Spatial resolution / algorithm for exemplar matching. Options:
#'   * `64` (or `"fast"`, `"vit"`, `"low"`): Fast 64×64 ViT patch token matching (stride 16). Extremely fast
#'     and ideal for large objects (e.g. people, leaves, fruits, animals).
#'   * `256` (or `"medium"`): 256×256 multi-scale CNN feature fusion (stride 4).
#'   * `1024` (or `"dense"`, `"conv"`, `"high"`, `"grain"`): Dense Stride-1 Convolutional Exemplar matching
#'     at full 1024×1024 pixel resolution with single-pass Zero-Mean Normalized Cross-Correlation (ZNCC) and
#'     ViT semantic guidance. Essential for resolving small touching objects, seeds, and grains with sharp pixel boundaries.
#'   Defaults to `256`.
#' @param max_objects Maximum number of instances to detect and segment. Default is `300`.
#' @param show_mask Deprecated. If `TRUE`, equivalent to `type = "mask"`.
#'
#' @return
#'   * If `mask = FALSE` and `bbox = TRUE`: A data frame containing the bounding box coordinates
#'     (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `score`, `label`).
#'   * If `type = "segment"`: An `image` object containing the segmented image.
#'   * If `type = "mask"`: A single multi-label or binary `image` object (`colormode = "Grayscale"`).
#'   * If `type = "highlight"`: Invisibly returns a list with:
#'     - `"boxes"`: Data frame of detected bounding boxes and scores.
#'     - `"labels"`: Single multi-label integer mask `image` where each pixel value corresponds to instance `id`.
#'     - `"summary"`: Text summary of detected objects and counts.
#'     - `"counts"`: Data frame of instance counts per class label.
#'     - `"contours"`: List of polygon contour coordinate matrices.
#'     - `"features"`: Morphological measures (`area`, `perimeter`, `solidity`, `convexity`, etc.) computed via `poly_measures()`.
#'     - `"similarity_map"`: Continuous cosine similarity `image` object (when `exemplar = TRUE`).
#'   * If `img` is a list: A list of the corresponding outputs.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # 1. Basic salient foreground segmentation with white background (CPU-friendly)
#'   seg <- image_segment_dl(img, model = "u2netp")
#'
#'   # 2. Razor-sharp plant boundary extraction (leaf margins, serrations, roots)
#'   seg_ben <- image_segment_dl(img, model = "ben2")
#'
#'   # 3. Semi-transparent polygon overlay for visual inspection
#'   image_segment_dl(img, model = "u2netp", type = "highlight")
#'
#'   # 4. Zero-shot text-prompted instance segmentation (Grounded-SAM)
#'   # Segments leaves and lesions, displaying bounding boxes and summary
#'   res_gs <- image_segment_dl(
#'     img,
#'     model = "grounded-sam",
#'     prompt = "plant leaf, lesion",
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#'   # Inspect counts and morphological measures
#'   print(attr(res_gs, "counts"))
#'   print(head(attr(res_gs, "features")))
#'
#'   # 5. Ultra-fast pure object detection without mask generation
#'   # (Skips SAM 2.1 completely for 10x faster execution)
#'   boxes <- image_segment_dl(
#'     img,
#'     model = "grounded-sam",
#'     prompt = "insect, leaf",
#'     mask = FALSE,
#'     bbox = TRUE
#'   )
#'
#'   # 6. One-shot visual exemplar segmentation (PerSAM):
#'   # Click on 1 or more representative items (e.g. seeds, logs, grains, cells)
#'   # to automatically find, outline, and count all matching instances:
#'   res_persam <- image_segment_dl(
#'     img,
#'     exemplar = TRUE,
#'     sim_threshold = 0.48,
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#'
#'   # 7. Multi-prototype visual exemplar segmentation with predefined coordinates:
#'   # Passing coordinates of distinct log types (light, dark, weathered)
#'   pts <- rbind(c(65, 380), c(425, 380), c(630, 380))
#'   res_logs <- image_segment_dl(
#'     img,
#'     exemplar = TRUE,
#'     prompt = pts,
#'     sim_threshold = 0.45,
#'     type = "highlight",
#'     engine = "gpu"
#'   )
#'   # 8. YOLO26 instance segmentation (real-time polygon mask & box detection)
#'   res_yolo <- image_segment_dl(
#'     img,
#'     model = "yolo26n-seg",
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#' }
image_segment_dl <- function(img,
                             type = c("segment", "mask", "highlight"),
                             col_background = "white",
                             col_highlight = "salmon",
                             alpha = 0.4,
                             border = NA,
                             lwd = 1,
                             bbox = FALSE,
                             mask = TRUE,
                             model = c("u2netp", "ben2", "yolo26n-seg", "stardist", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                             threshold = 0.5,
                             box_threshold = 0.25,
                             text_threshold = 0.25,
                             iou_threshold = 0.5,
                             exemplar = FALSE,
                             sim_threshold = 0.5,
                             min_dist = 16,
                             feat_res = 256,
                             max_objects = 300,
                             threads = 0,
                             engine = c("cpu", "gpu"),
                             device_id = -1,
                             fill_hull = TRUE,
                             filter = 0,
                             erode = 0,
                             dilate = 0,
                             opening = 0,
                             closing = 0,
                             min_area = 0,
                             invert = FALSE,
                             pick_object = FALSE,
                             prompt = NULL,
                             show_mask = FALSE,
                             verbose = TRUE,
                             plot = TRUE,
                             dir = pliman_model_dir(),
                             ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_segment_dl(x,
                       type = type,
                       col_background = col_background,
                       col_highlight = col_highlight,
                       alpha = alpha,
                       border = border,
                       lwd = lwd,
                       bbox = bbox,
                       mask = mask,
                       model = model,
                       threshold = threshold,
                       box_threshold = box_threshold,
                       text_threshold = text_threshold,
                       iou_threshold = iou_threshold,
                       exemplar = exemplar,
                       sim_threshold = sim_threshold,
                       min_dist = min_dist,
                       feat_res = feat_res,
                       max_objects = max_objects,
                       threads = threads,
                       engine = engine,
                       device_id = device_id,
                       fill_hull = fill_hull,
                       filter = filter,
                       erode = erode,
                       dilate = dilate,
                       opening = opening,
                       closing = closing,
                       min_area = min_area,
                       invert = invert,
                       pick_object = pick_object,
                       prompt = prompt,
                       show_mask = show_mask,
                       verbose = verbose,
                       plot = FALSE,
                       dir = dir,
                       ...)
    })
    return(res)
  }

  if (!isTRUE(mask) && !isTRUE(bbox)) {
    cli::cli_abort("At least one of {.arg mask} or {.arg bbox} must be TRUE.")
  }

  if (isTRUE(show_mask)) {
    type <- "mask"
  } else {
    type_str <- tolower(type[1])
    if (grepl("^high", type_str)) {
      type <- "highlight"
    } else {
      type <- match.arg(type, c("segment", "mask", "highlight"))
    }
  }

  if (is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Grounded-SAM, PerSAM, YOLO-Seg, and StarDist instance segmentation
  is_text_prompt <- is.character(prompt) && length(prompt) >= 1 && !all(prompt %in% c("center", "box", "exemplar"))
  is_sam_model <- grepl("sam", model, ignore.case = TRUE)
  use_persam <- isTRUE(exemplar) || (model == "persam") || (is_sam_model && identical(prompt, "exemplar"))
  use_grounded_sam <- !use_persam && ((model == "grounded-sam") || (is_sam_model && is_text_prompt))
  use_yolo_seg <- !use_persam && !use_grounded_sam && (model %in% c("yolo26n-seg", "yolo11n-seg", "yolo-seg") || grepl("-seg", model))
  use_stardist <- !use_persam && !use_grounded_sam && !use_yolo_seg && (model %in% c("stardist", "star-dist"))

  if (use_persam || use_grounded_sam || use_yolo_seg || use_stardist) {
    if (use_persam) {
      if (!missing(box_threshold) && missing(sim_threshold)) {
        sim_threshold <- box_threshold
      }
      dots <- list(...)
      ex_pts <- if (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt)) {
        prompt
      } else if (!is.null(dots$exemplar_points)) {
        dots$exemplar_points
      } else {
        NULL
      }
      gs_res <- .run_persam(
        mat = mat,
        exemplar_points = ex_pts,
        sim_threshold = sim_threshold,
        min_dist = min_dist,
        iou_threshold = iou_threshold,
        max_objects = max_objects,
        feat_res = feat_res,
        threshold = threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        invert = invert,
        mask = mask,
        verbose = verbose,
        dir = dir
      )
    } else if (use_grounded_sam) {
      gs_res <- .run_grounded_sam(
        mat = mat,
        prompt = prompt,
        threshold = threshold,
        box_threshold = box_threshold,
        text_threshold = text_threshold,
        iou_threshold = iou_threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        invert = invert,
        mask = mask,
        verbose = verbose,
        dir = dir
      )
    } else if (use_yolo_seg) {
      gs_res <- .run_yolo(
        mat = mat,
        model = model,
        conf_threshold = box_threshold,
        iou_threshold = iou_threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        dir = dir
      )
    } else if (use_stardist) {
      gs_res <- .run_stardist(
        mat = mat,
        model = model,
        prob_threshold = threshold,
        nms_threshold = iou_threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        dir = dir
      )
    }

    sum_title <- if (use_persam) {
      "PerSAM Detection Summary"
    } else if (use_yolo_seg) {
      "YOLO Detection Summary"
    } else if (use_stardist) {
      "StarDist Detection Summary"
    } else {
      "Grounded-SAM Detection Summary"
    }

    if (!isTRUE(mask) && isTRUE(bbox)) {
      if (isTRUE(plot)) {
        plot(as_image(mat), ...)
        num_inst <- nrow(gs_res$boxes)
        if (num_inst > 0) {
          palette_colors <- if (num_inst == 1) col_highlight else grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
          .plot_yolo_bboxes(
            boxes = gs_res$boxes,
            palette_colors = palette_colors,
            lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd
          )
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(gs_res$boxes, "summary") <- gs_res$summary
      attr(gs_res$boxes, "counts") <- gs_res$counts
      if (!is.null(gs_res$similarity_map)) attr(gs_res$boxes, "similarity_map") <- gs_res$similarity_map
      return(invisible(gs_res$boxes))
    }

    if (type == "mask") {
      if (!isTRUE(mask) || is.null(gs_res$labels)) {
        cli::cli_alert_warning("type = 'mask' requested but mask = FALSE. Returning bounding boxes instead.")
        return(invisible(gs_res$boxes))
      }
      out_img <- gs_res$labels
      if (isTRUE(plot)) {
        plot(out_img, ...)
        if (isTRUE(bbox) && nrow(gs_res$boxes) > 0) {
          num_inst <- nrow(gs_res$boxes)
          palette_colors <- if (num_inst == 1) col_highlight else grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
          for (k in seq_len(num_inst)) {
            bx <- gs_res$boxes[k, ]
            k_col <- palette_colors[((k - 1) %% length(palette_colors)) + 1]
            graphics::rect(xleft = bx$xmin, ybottom = bx$ymin, xright = bx$xmax, ytop = bx$ymax,
                           border = k_col, lwd = lwd)
          }
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(out_img, "boxes") <- gs_res$boxes
      attr(out_img, "summary") <- gs_res$summary
      attr(out_img, "counts") <- gs_res$counts
      if (!is.null(gs_res$similarity_map)) attr(out_img, "similarity_map") <- gs_res$similarity_map
      return(invisible(out_img))
    }

    if (type == "segment") {
      out <- as_image(mat)
      if (isTRUE(mask) && !is.null(gs_res$mask)) {
        col_rgb <- tryCatch(grDevices::col2rgb(col_background) / 255.0, error = function(e) c(1, 1, 1))
        m_bg <- gs_res$mask
        if (length(dims) == 3) {
          for (k in seq_len(dims[3])) {
            fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
            out[, , k][!m_bg] <- fill_val
          }
        } else {
          fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
          out[!m_bg] <- fill_val
        }
      }

      if (isTRUE(plot)) {
        plot(out, ...)
        if ((isTRUE(bbox) || !isTRUE(mask)) && nrow(gs_res$boxes) > 0) {
          num_inst <- nrow(gs_res$boxes)
          palette_colors <- if (num_inst == 1) col_highlight else grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
          .plot_yolo_bboxes(
            boxes = gs_res$boxes,
            palette_colors = palette_colors,
            lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd
          )
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(out, "boxes") <- gs_res$boxes
      attr(out, "labels") <- gs_res$labels
      attr(out, "summary") <- gs_res$summary
      attr(out, "counts") <- gs_res$counts
      if (!is.null(gs_res$similarity_map)) attr(out, "similarity_map") <- gs_res$similarity_map
      return(invisible(out))
    }

    if (type == "highlight") {
      num_inst <- nrow(gs_res$boxes)
      conts_all <- vector("list", num_inst)

      if (isTRUE(mask) && !is.null(gs_res$labels) && num_inst > 0) {
        lbl_mat <- image_data(gs_res$labels)
        has_boxes <- !is.null(gs_res$boxes) && is.data.frame(gs_res$boxes) &&
                     nrow(gs_res$boxes) == num_inst && all(c("xmin", "ymin", "xmax", "ymax") %in% names(gs_res$boxes))

        pad <- 2L
        for (k in seq_len(num_inst)) {
          if (has_boxes) {
            bx <- gs_res$boxes[k, ]
            x1 <- max(1L, floor(bx$xmin) - pad)
            x2 <- min(orig_w, ceiling(bx$xmax) + pad)
            y1 <- max(1L, floor(bx$ymin) - pad)
            y2 <- min(orig_h, ceiling(bx$ymax) + pad)
            sub_lbl <- lbl_mat[x1:x2, y1:y2, drop = FALSE]
            sub_m <- (sub_lbl == k)
            if (!any(sub_m)) next
            lbls <- bwlabel_cpp(sub_m)
            conts <- list()
            if (max(lbls) > 0) {
              conts <- contour(lbls)
              for (ci in seq_along(conts)) {
                if (is.matrix(conts[[ci]])) {
                  conts[[ci]][, 1] <- conts[[ci]][, 1] + (x1 - 1L)
                  conts[[ci]][, 2] <- conts[[ci]][, 2] + (y1 - 1L)
                }
              }
            }
          } else {
            inst_m <- (lbl_mat == k)
            if (!any(inst_m)) next
            lbls <- bwlabel_cpp(inst_m)
            conts <- if (max(lbls) > 0) contour(lbls) else list()
          }
          conts_all[[k]] <- conts
        }
      }

      if (isTRUE(plot)) {
        if (isTRUE(verbose)) {
          cli::cli_progress_step(
            msg = if (isTRUE(mask)) "Rendering highlight overlay..." else "Rendering detection overlay...",
            msg_done = if (isTRUE(mask)) "Highlight overlay ready" else "Detection overlay ready"
          )
        }
        plot(as_image(mat), ...)
        if (num_inst > 0) {
          palette_colors <- if (num_inst == 1) {
            col_highlight
          } else {
            grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
          }

          if (isTRUE(mask) && !is.null(gs_res$labels)) {
            for (k in seq_len(num_inst)) {
              conts <- conts_all[[k]]
              if (is.null(conts) || length(conts) == 0) next

              k_col <- palette_colors[((k - 1) %% length(palette_colors)) + 1]
              c_rgb <- tryCatch(grDevices::col2rgb(k_col) / 255.0, error = function(e) c(0, 1, 0))
              poly_col <- grDevices::rgb(c_rgb[1], c_rgb[2], c_rgb[3], alpha = alpha)
              cur_border <- if (is.na(border)) NA else k_col

              for (p in conts) {
                if (is.matrix(p) && nrow(p) >= 3) {
                  p_plot <- .decimate_poly(p)
                  graphics::polygon(p_plot[, 1], p_plot[, 2], col = poly_col, border = cur_border, lwd = lwd)
                }
              }
            }
          }

          if (isTRUE(bbox) || !isTRUE(mask)) {
            .plot_yolo_bboxes(
              boxes = gs_res$boxes,
              palette_colors = palette_colors,
              lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd
            )
          } else {
            cx <- (gs_res$boxes$xmin + gs_res$boxes$xmax) / 2
            cy <- (gs_res$boxes$ymin + gs_res$boxes$ymax) / 2
            graphics::points(cx, cy, pch = 21, bg = "black", col = "white", cex = 2.2)
            graphics::text(cx, cy, labels = as.character(seq_len(num_inst)), col = "white", font = 2, cex = 0.9)
          }
        }
        if (isTRUE(verbose)) {
          cli::cli_progress_done()
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      flat_conts <- unlist(conts_all, recursive = FALSE)
      if (is.null(flat_conts)) flat_conts <- list()
      attr(gs_res, "contours") <- flat_conts
      attr(gs_res, "summary") <- gs_res$summary
      attr(gs_res, "counts") <- gs_res$counts
      if (!is.null(gs_res$similarity_map)) attr(gs_res, "similarity_map") <- gs_res$similarity_map
      return(invisible(gs_res))
    }
  }

  # 1. Compute refined mask with all morphological operations applied
  mask_mat <- .compute_dl_mask(mat = mat,
                               model = model,
                               threshold = threshold,
                               fill_hull = fill_hull,
                               filter = filter,
                               erode = erode,
                               dilate = dilate,
                               opening = opening,
                               closing = closing,
                               min_area = min_area,
                               invert = invert,
                               pick_object = pick_object,
                               prompt = prompt,
                               exemplar = exemplar,
                               threads = threads,
                               engine = engine,
                               device_id = device_id,
                               verbose = verbose,
                               dir = dir)

  if (!isTRUE(mask) && isTRUE(bbox)) {
    lbls <- bwlabel_cpp(mask_mat > 0)
    conts <- if (max(lbls) > 0) contour(lbls) else list()
    num_inst <- length(conts)
    if (num_inst > 0) {
      boxes_list <- lapply(seq_along(conts), function(i) {
        p <- conts[[i]]
        data.frame(
          id = i,
          xmin = round(min(p[, 1]), 1),
          ymin = round(min(p[, 2]), 1),
          xmax = round(max(p[, 1]), 1),
          ymax = round(max(p[, 2]), 1),
          score = 1.0,
          label = "object",
          stringsAsFactors = FALSE
        )
      })
      boxes_df <- do.call(rbind, boxes_list)
    } else {
      boxes_df <- data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                             xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                             label = character(0), stringsAsFactors = FALSE)
    }
    if (isTRUE(plot)) {
      plot(as_image(mat), ...)
      if (nrow(boxes_df) > 0) {
        palette_colors <- if (nrow(boxes_df) == 1) col_highlight else grDevices::rainbow(nrow(boxes_df), s = 0.85, v = 0.95)
        .plot_yolo_bboxes(
          boxes = boxes_df,
          palette_colors = palette_colors,
          lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd
        )
      }
    }
    return(invisible(boxes_df))
  }

  # 2. Output according to selected type
  if (type == "mask") {
    mask_img <- as_image(mask_mat, colormode = "Grayscale")
    if (isTRUE(plot)) {
      plot(mask_img, ...)
    }
    return(invisible(mask_img))
  }

  if (type == "segment") {
    out <- as_image(mat)
    col_rgb <- tryCatch(grDevices::col2rgb(col_background) / 255.0, error = function(e) c(1, 1, 1))
    if (length(dims) == 3) {
      for (k in seq_len(dims[3])) {
        fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
        out[, , k][!mask_mat] <- fill_val
      }
    } else {
      fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
      out[!mask_mat] <- fill_val
    }

    if (isTRUE(plot)) {
      plot(out, ...)
    }
    return(invisible(out))
  }

  if (type == "highlight") {
    lbls <- bwlabel_cpp(mask_mat > 0)
    conts <- if (max(lbls) > 0) contour(lbls) else list()

    if (isTRUE(plot)) {
      if (isTRUE(verbose)) {
        cli::cli_progress_step(
          msg = "Rendering highlight overlay...",
          msg_done = "Highlight overlay ready"
        )
      }
      plot(as_image(mat), ...)
      col_rgb <- tryCatch(grDevices::col2rgb(col_highlight) / 255.0, error = function(e) c(0, 1, 0))
      poly_col <- grDevices::rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
      for (p in conts) {
        if (is.matrix(p) && nrow(p) >= 3) {
          p_plot <- .decimate_poly(p)
          graphics::polygon(p_plot[, 1], p_plot[, 2], col = poly_col, border = border, lwd = lwd)
        }
        if (isTRUE(bbox) && is.matrix(p) && nrow(p) >= 1) {
          graphics::rect(xleft = min(p[, 1]), ybottom = min(p[, 2]),
                         xright = max(p[, 1]), ytop = max(p[, 2]),
                         border = col_highlight, lwd = lwd)
        }
      }
      if (isTRUE(verbose)) {
        cli::cli_progress_done()
      }
    }
    return(invisible(conts))
  }
}

#' Clear Cached ONNX Runtime Sessions
#'
#' Releases all in-memory compiled neural network sessions (e.g., Grounded-SAM,
#' SAM 2.1, U2-Net) and frees allocated RAM.
#'
#' @return Invisible \code{NULL}.
#' @export
pliman_clear_sessions <- function() {
  clear_onnx_sessions_cpp()
  if (exists(".pliman_vocab_cache")) {
    rm(list = ls(envir = .pliman_vocab_cache), envir = .pliman_vocab_cache)
  }
  invisible(NULL)
}

#' @rdname image_remove_bg_dl
#' @export
image_remove_bg <- function(img,
                            model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                            threshold = 0.5,
                            fill_hull = TRUE,
                            filter = 0,
                            erode = 0,
                            dilate = 0,
                            opening = 0,
                            closing = 0,
                            min_area = 0,
                            invert = FALSE,
                            pick_object = FALSE,
                            prompt = NULL,
                            exemplar = FALSE,
                            threads = 0,
                            engine = c("cpu", "gpu"),
                            device_id = -1,
                            verbose = TRUE,
                            transparent = TRUE,
                            bg_color = "black",
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  image_remove_bg_dl(img = img,
                     model = model,
                     threshold = threshold,
                     fill_hull = fill_hull,
                     filter = filter,
                     erode = erode,
                     dilate = dilate,
                     opening = opening,
                     closing = closing,
                     min_area = min_area,
                     invert = invert,
                     pick_object = pick_object,
                     prompt = prompt,
                     exemplar = exemplar,
                     threads = threads,
                     engine = engine,
                     device_id = device_id,
                     verbose = verbose,
                     transparent = transparent,
                     bg_color = bg_color,
                     plot = plot,
                     dir = dir,
                     ...)
}

#' Estimate Relative 3D Depth with Depth Anything V2
#'
#' Estimates a metric-aligned relative 3D depth map from any 2D RGB image using
#' the state-of-the-art Depth Anything V2 foundation model running purely via ONNX Runtime.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"depth-anything-v2"`.
#' @param col_palette Color palette used to visualize the depth map when `plot = TRUE`.
#'   Options include `"viridis"`, `"magma"`, `"inferno"`, `"plasma"`, `"spectral"`, or a custom vector of colors.
#' @param invert Logical. If `TRUE`, inverts the depth map so near objects are darker. Defaults to `FALSE`.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the resulting depth heatmap.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return An `Image` object representing the colorized depth map, with the raw normalized
#'   depth matrix stored in the `"depth"` attribute (`attr(res, "depth")`).
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'   depth <- image_depth_dl(img)
#' }
image_depth_dl <- function(img,
                           model = "depth-anything-v2",
                           col_palette = "magma",
                           invert = FALSE,
                           threads = 0,
                           engine = c("cpu", "gpu"),
                           device_id = -1,
                           verbose = TRUE,
                           plot = TRUE,
                           dir = pliman_model_dir(),
                           ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_depth_dl(x, model = model, col_palette = col_palette, invert = invert,
                     threads = threads, engine = engine, device_id = device_id,
                     verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "depth-anything-v2"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running Depth Anything V2 [{toupper(engine)}]...",
      msg_done = "Depth estimation complete"
    )
  }

  in_size <- 518L
  tensor <- .preprocess_nchw(
    mat,
    target_size = in_size,
    mean = c(0.485, 0.456, 0.406),
    std = c(0.229, 0.224, 0.225),
    letterbox = FALSE
  )

  depth_mat <- run_depth_anything_cpp(
    tensor_vec = tensor,
    in_w = in_size,
    in_h = in_size,
    orig_w = orig_w,
    orig_h = orig_h,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  min_d <- min(depth_mat, na.rm = TRUE)
  max_d <- max(depth_mat, na.rm = TRUE)
  rng <- if (max_d > min_d) (max_d - min_d) else 1.0
  norm_d <- (depth_mat - min_d) / rng
  if (isTRUE(invert)) norm_d <- 1.0 - norm_d

  pal_colors <- switch(col_palette,
    "viridis"  = grDevices::hcl.colors(256, "Viridis"),
    "magma"    = grDevices::hcl.colors(256, "Inferno"),
    "inferno"  = grDevices::hcl.colors(256, "Inferno"),
    "plasma"   = grDevices::hcl.colors(256, "Plasma"),
    "spectral" = grDevices::hcl.colors(256, "Spectral"),
    if (is.character(col_palette) && length(col_palette) > 1) {
      grDevices::colorRampPalette(col_palette)(256)
    } else {
      grDevices::hcl.colors(256, "Viridis")
    }
  )

  rgb_mat <- grDevices::col2rgb(pal_colors) / 255.0
  idx <- pmax(1L, pmin(256L, round(norm_d * 255.0) + 1L))

  out_arr <- array(0.0, dim = c(orig_w, orig_h, 3L))
  out_arr[, , 1] <- matrix(rgb_mat[1, idx], nrow = orig_w, ncol = orig_h)
  out_arr[, , 2] <- matrix(rgb_mat[2, idx], nrow = orig_w, ncol = orig_h)
  out_arr[, , 3] <- matrix(rgb_mat[3, idx], nrow = orig_w, ncol = orig_h)

  res_img <- as_image(out_arr)
  attr(res_img, "depth") <- depth_mat

  if (isTRUE(plot)) {
    plot(res_img, ...)
  }

  invisible(res_img)
}

#' Extract Foundation Features and Semantic Visualization with DINOv2
#'
#' Extracts dense patch embeddings and global representation vectors from DINOv2
#' (Vision Transformer ViT-S/14), and computes top-3 PCA false-color semantic segmentation
#' without any labels or supervised training.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"dinov2"`.
#' @param patch_size Size of ViT patches in pixels (default 14).
#' @param return_pca Logical. If `TRUE`, computes top-3 principal components of patch embeddings.
#' @param interpolate Logical. If `TRUE`, interpolates the PCA map from patch grid (37x37) to original image resolution.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the PCA false-color semantic image.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list containing:
#'   * `cls_token`: 384-dimensional global semantic representation vector.
#'   * `pca_image`: An `Image` object representing the 3-component PCA false-color semantic map.
#'   * `pca_r`, `pca_g`, `pca_b`: Matrices of normalized PCA coordinates.
#'   * `num_patches`: Total number of spatial patches.
#'   * `embed_dim`: Dimensionality of patch embeddings.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'   feats <- image_features_dl(img)
#' }
image_features_dl <- function(img,
                              model = "dinov2",
                              patch_size = 14,
                              return_pca = TRUE,
                              interpolate = TRUE,
                              threads = 0,
                              engine = c("cpu", "gpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_features_dl(x, model = model, patch_size = patch_size, return_pca = return_pca,
                        interpolate = interpolate, threads = threads, engine = engine,
                        device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "dinov2"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running DINOv2 foundation ViT [{toupper(engine)}]...",
      msg_done = "DINOv2 feature extraction complete"
    )
  }

  in_size <- 518L
  tensor <- .preprocess_nchw(
    mat,
    target_size = in_size,
    mean = c(0.485, 0.456, 0.406),
    std = c(0.229, 0.224, 0.225),
    letterbox = FALSE
  )

  dino_res <- run_dinov2_cpp(
    tensor_vec = tensor,
    in_w = in_size,
    in_h = in_size,
    patch_size = as.integer(patch_size),
    return_pca = isTRUE(return_pca),
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  pca_img <- NULL
  if (isTRUE(return_pca)) {
    wp <- dino_res$wp
    hp <- dino_res$hp
    r_mat <- dino_res$pca_r
    g_mat <- dino_res$pca_g
    b_mat <- dino_res$pca_b

    if (isTRUE(interpolate)) {
      r_mat <- .bilinear_resize_2d(r_mat, orig_w, orig_h)
      g_mat <- .bilinear_resize_2d(g_mat, orig_w, orig_h)
      b_mat <- .bilinear_resize_2d(b_mat, orig_w, orig_h)
      out_w <- orig_w
      out_h <- orig_h
    } else {
      out_w <- wp
      out_h <- hp
    }

    arr <- array(0.0, dim = c(out_w, out_h, 3L))
    arr[, , 1] <- r_mat
    arr[, , 2] <- g_mat
    arr[, , 3] <- b_mat
    pca_img <- as_image(arr)

    if (isTRUE(plot)) {
      plot(pca_img, ...)
    }
  }

  out <- list(
    cls_token = dino_res$cls_token,
    pca_image = pca_img,
    pca_r = dino_res$pca_r,
    pca_g = dino_res$pca_g,
    pca_b = dino_res$pca_b,
    num_patches = dino_res$num_patches,
    embed_dim = dino_res$embed_dim
  )
  invisible(out)
}

#' Universal Object Detection with YOLO26
#'
#' Performs real-time object detection using YOLO26 (or custom YOLO ONNX models)
#' with bounding box regression and End-to-End detection.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n"`.
#'   Options include `"yolo26n"`, `"yolo26s"`, `"yolo26m"`, `"yolo26l"`, `"yolo26x"`.
#' @param conf_threshold Minimum confidence score for candidate boxes (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param labels Optional character vector of class names (defaults to standard 80 COCO classes).
#' @param col Optional color or palette for bounding boxes.
#' @param lwd Line width for bounding boxes (default 2).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the original image with bounding boxes.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A data frame containing detected bounding boxes (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`),
#'   with `summary` and `counts` stored as attributes.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("objects.png")
#'   boxes <- image_detect_dl(img)
#' }
image_detect_dl <- function(img,
                            model = "yolo26n",
                            conf_threshold = 0.25,
                            iou_threshold = 0.45,
                            labels = NULL,
                            col = NULL,
                            lwd = 2,
                            threads = 0,
                            engine = c("cpu", "gpu"),
                            device_id = -1,
                            verbose = TRUE,
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_detect_dl(x, model = model, conf_threshold = conf_threshold, iou_threshold = iou_threshold,
                      labels = labels, col = col, lwd = lwd, threads = threads, engine = engine,
                      device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "yolo26n"

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running YOLO object detection [{toupper(engine)}]...",
      msg_done = "Object detection complete"
    )
  }

  res <- .run_yolo(
    mat = mat,
    model = model_str,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    labels = labels,
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  df_boxes <- res$boxes

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    num_inst <- nrow(df_boxes)
    if (num_inst > 0) {
      palette_colors <- if (!is.null(col)) {
        rep(col, length.out = num_inst)
      } else if (num_inst == 1) {
        "salmon"
      } else {
        grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
      }
      .plot_yolo_bboxes(
        boxes = df_boxes,
        palette_colors = palette_colors,
        lwd = lwd
      )
    }
  }

  if (isTRUE(verbose)) {
    .print_detection_summary(res$summary, title = "YOLO Object Detection Summary")
  }

  attr(df_boxes, "summary") <- res$summary
  attr(df_boxes, "counts") <- res$counts
  invisible(df_boxes)
}

#' Universal Human Pose Estimation with YOLO26
#'
#' Performs real-time human pose estimation and 17 keypoint detection using YOLO26
#' (e.g., `yolo26n-pose`, `yolo26s-pose`, etc.) with bounding box regression and
#' anatomical skeleton visualization.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n-pose"`.
#'   Options include `"yolo26n-pose"`, `"yolo26s-pose"`, `"yolo26m-pose"`, `"yolo26l-pose"`, `"yolo26x-pose"`.
#' @param conf_threshold Minimum confidence score for candidate person detections (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param kpt_threshold Minimum confidence threshold for rendering individual keypoints and skeleton limbs (default 0.3).
#' @param col Optional color or palette for bounding boxes and keypoint markers.
#' @param lwd Line width for bounding boxes and skeleton limbs (default 2).
#' @param kpt_radius Radius/size for keypoint markers (default 4).
#' @param bbox Logical. Whether to draw bounding boxes around detected persons (default `TRUE`).
#' @param skeleton Logical. Whether to render the 17-keypoint anatomical skeleton (default `TRUE`).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the image with pose overlays.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list with class `pliman_pose` containing:
#'   * `boxes`: Data frame containing detected person bounding boxes (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`).
#'   * `keypoints`: List of data frames (one per person) containing the 17 COCO keypoints (`keypoint`, `x`, `y`, `conf`).
#'   * `summary`: Human-readable detection summary string.
#'   * `counts`: Data frame of instance counts.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("person.png")
#'   res <- image_pose_dl(img)
#'   res$boxes
#'   res$keypoints[[1]]
#' }
image_pose_dl <- function(img,
                          model = "yolo26n-pose",
                          conf_threshold = 0.25,
                          iou_threshold = 0.45,
                          kpt_threshold = 0.3,
                          col = NULL,
                          lwd = 2,
                          kpt_radius = 4,
                          bbox = TRUE,
                          skeleton = TRUE,
                          threads = 0,
                          engine = c("cpu", "gpu"),
                          device_id = -1,
                          verbose = TRUE,
                          plot = TRUE,
                          dir = pliman_model_dir(),
                          ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_pose_dl(x, model = model, conf_threshold = conf_threshold, iou_threshold = iou_threshold,
                    kpt_threshold = kpt_threshold, col = col, lwd = lwd, kpt_radius = kpt_radius,
                    bbox = bbox, skeleton = skeleton, threads = threads, engine = engine,
                    device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "yolo26n-pose"

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running YOLO pose estimation [{toupper(engine)}]...",
      msg_done = "Pose estimation complete"
    )
  }

  res <- .run_yolo(
    mat = mat,
    model = model_str,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    labels = "person",
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  df_boxes <- res$boxes
  num_inst <- nrow(df_boxes)

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    if (num_inst > 0) {
      palette_colors <- if (!is.null(col)) {
        rep(col, length.out = num_inst)
      } else if (num_inst == 1) {
        "salmon"
      } else {
        grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
      }

      if (isTRUE(bbox)) {
        .plot_yolo_bboxes(
          boxes = df_boxes,
          palette_colors = palette_colors,
          lwd = lwd
        )
      }

      if (isTRUE(skeleton) && length(res$keypoints) > 0) {
        .plot_yolo_keypoints(
          keypoints_list = res$keypoints,
          palette_colors = palette_colors,
          kpt_threshold = kpt_threshold,
          lwd = lwd,
          kpt_radius = kpt_radius
        )
      }
    }
  }

  if (isTRUE(verbose)) {
    .print_detection_summary(res$summary, title = "YOLO Pose Estimation Summary")
  }

  out <- list(
    boxes = df_boxes,
    keypoints = res$keypoints,
    summary = res$summary,
    counts = res$counts
  )
  class(out) <- c("pliman_pose", "list")
  invisible(out)
}

#' Real-Time Image Classification with YOLO26
#'
#' Classifies images using pre-trained YOLO26 classification models across 1,000 ImageNet categories.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n-cls"`.
#'   Options include `"yolo26n-cls"`, `"yolo26s-cls"`, `"yolo26m-cls"`, `"yolo26l-cls"`.
#' @param top_k Integer. Number of top class predictions to return (default 5).
#' @param labels Optional character vector of class names (defaults to 1,000 standard ImageNet classes).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to display the image with top predicted class labels.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A data frame containing top predictions (`rank`, `class`, `probability`).
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("sample.png")
#'   image_classify_dl(img)
#' }
image_classify_dl <- function(img,
                              model = "yolo26n-cls",
                              top_k = 5,
                              labels = NULL,
                              threads = 0,
                              engine = c("cpu", "gpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = FALSE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_classify_dl(x, model = model, top_k = top_k, labels = labels,
                        threads = threads, engine = engine, device_id = device_id,
                        verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "yolo26n-cls"

  if (file.exists(model_str)) {
    model_file <- normalizePath(model_str, winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  tensor <- .preprocess_yolo(mat, target_size = 640L)

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running YOLO image classification [{toupper(engine)}]...",
      msg_done = "Image classification complete"
    )
  }

  probs <- run_yolo_cls_cpp(
    tensor_vec = tensor,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  class_names <- if (!is.null(labels)) labels else .imagenet_classes
  num_classes <- length(probs)
  k <- min(as.integer(top_k), num_classes)

  top_idx <- order(probs, decreasing = TRUE)[seq_len(k)]
  top_probs <- round(probs[top_idx], 4)
  top_labels <- ifelse(top_idx <= length(class_names), class_names[top_idx], paste0("class_", top_idx - 1L))

  df_res <- data.frame(
    rank = seq_len(k),
    class = top_labels,
    probability = top_probs,
    class_id = top_idx - 1L,
    stringsAsFactors = FALSE
  )

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    top1 <- df_res[1, ]
    graphics::title(sub = sprintf("Top 1: %s (%.1f%%)", top1$class, top1$probability * 100), col.sub = "darkgreen", font.sub = 2)
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Top-{k} Classification Results:")
    for (i in seq_len(nrow(df_res))) {
      cli::cli_bullets(c("*" = sprintf("%d. {.val %s}: %.2f%%", df_res$rank[i], df_res$class[i], df_res$probability[i] * 100)))
    }
  }

  invisible(df_res)
}

#' Star-Convex Object Detection with StarDist
#'
#' Detects round, convex, and overlapping objects (e.g. seeds, cells, nuclei, spores)
#' by predicting radial star-convex polygons and object probability using StarDist.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"stardist"`.
#' @param prob_threshold Probability threshold for detecting object centers (default 0.5).
#' @param nms_threshold Polygon IoU threshold for Non-Maximum Suppression (default 0.3).
#' @param type Output type: `"segment"` (highlight overlay), `"mask"` (integer labels), or `"polygons"`.
#' @param col_highlight Fill color for overlay polygons (default `"salmon"`).
#' @param border Border color for polygon contours (default `"white"`).
#' @param lwd Line width for polygon contours (default 2).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the results.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list containing:
#'   * `boxes`: Bounding boxes data frame with centers and scores.
#'   * `polygons_x`, `polygons_y`: Coordinates of radial polygon vertices for each object.
#'   * `labels`: Integer matrix with instance IDs.
#'   * `mask`: Logical foreground mask.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("grains.png")
#'   res <- image_stardist_dl(img)
#' }
image_stardist_dl <- function(img,
                              model = "stardist",
                              prob_threshold = 0.5,
                              nms_threshold = 0.3,
                              type = c("segment", "mask", "polygons"),
                              col_highlight = "salmon",
                              border = "white",
                              lwd = 2,
                              threads = 0,
                              engine = c("cpu", "gpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  type <- match.arg(type)
  engine <- match.arg(engine)

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_stardist_dl(x, model = model, prob_threshold = prob_threshold, nms_threshold = nms_threshold,
                        type = type, col_highlight = col_highlight, border = border, lwd = lwd,
                        threads = threads, engine = engine, device_id = device_id, verbose = verbose,
                        plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "stardist"

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running StarDist polygon detection [{toupper(engine)}]...",
      msg_done = "StarDist detection complete"
    )
  }

  res <- .run_stardist(
    mat = mat,
    model = model_str,
    prob_threshold = prob_threshold,
    nms_threshold = nms_threshold,
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  num_inst <- nrow(res$boxes)

  if (isTRUE(plot)) {
    if (type == "mask") {
      plot(as_image(res$labels, colormode = "Grayscale"), ...)
    } else {
      plot(as_image(mat), ...)
      if (num_inst > 0) {
        palette_colors <- if (num_inst == 1) col_highlight else grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
        for (i in seq_len(num_inst)) {
          px <- res$polygons_x[[i]]
          py <- res$polygons_y[[i]]
          k_col <- palette_colors[((i - 1) %% length(palette_colors)) + 1]
          col_poly <- grDevices::adjustcolor(k_col, alpha.f = 0.4)
          graphics::polygon(x = c(px, px[1]), y = c(py, py[1]), col = col_poly, border = border, lwd = lwd)
        }
      }
    }
  }

  if (isTRUE(verbose)) {
    .print_detection_summary(res$summary, title = "StarDist Object Detection Summary")
  }

  invisible(res)
}

#' 4x Super-Resolution with Real-ESRGAN Compact
#'
#' Upscales images by 4x using generative convolutional super-resolution with
#' intelligent tiled overlapping to prevent GPU/CPU memory exhaustion.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"realesrgan-compact"`.
#' @param scale Upscaling factor (default 4).
#' @param tile_size Processing tile size in pixels (default 256).
#' @param tile_pad Overlap margin in pixels for seamless blending (default 16).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the super-resolved image.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return An `Image` object with dimensions `(orig_w * 4) x (orig_h * 4) x 3`.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("low_res.png")
#'   img_hr <- image_superres_dl(img)
#' }
image_superres_dl <- function(img,
                              model = "realesrgan-compact",
                              scale = 4,
                              tile_size = 256,
                              tile_pad = 16,
                              threads = 0,
                              engine = c("cpu", "gpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_superres_dl(x, model = model, scale = scale, tile_size = tile_size,
                        tile_pad = tile_pad, threads = threads, engine = engine,
                        device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "realesrgan-compact"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running Real-ESRGAN {scale}x super-resolution [{toupper(engine)}]...",
      msg_done = "Super-resolution complete"
    )
  }

  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  tensor <- c(R, G, B)

  raw_vec <- run_super_resolution_cpp(
    tensor_vec = tensor,
    in_w = orig_w,
    in_h = orig_h,
    scale = as.integer(scale),
    tile_size = as.integer(tile_size),
    tile_pad = as.integer(tile_pad),
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  out_w <- orig_w * as.integer(scale)
  out_h <- orig_h * as.integer(scale)
  plane <- out_w * out_h

  out_arr <- array(0.0, dim = c(out_w, out_h, 3L))
  out_arr[, , 1] <- matrix(raw_vec[1:plane], nrow = out_w, ncol = out_h)
  out_arr[, , 2] <- matrix(raw_vec[(plane + 1):(2 * plane)], nrow = out_w, ncol = out_h)
  out_arr[, , 3] <- matrix(raw_vec[(2 * plane + 1):(3 * plane)], nrow = out_w, ncol = out_h)

  out_img <- as_image(out_arr)

  if (isTRUE(plot)) {
    plot(out_img, ...)
  }

  invisible(out_img)
}

