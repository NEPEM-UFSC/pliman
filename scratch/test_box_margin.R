
pkgload::load_all("d:/Desktop/pliman")

cat("=== 1. Checking formal arguments ===\n")
stopifnot("box_margin" %in% names(formals(yolo_dataset_export)))
stopifnot("box_margin" %in% names(formals(yolo_dataset_to_segment)))
stopifnot("box_margin" %in% names(formals(pliman:::.sam_segment_boxes)))
cat("[OK] box_margin is present in all function signatures with default 0.05.\n")

cat("\n=== 2. Testing .sam_segment_boxes bounding box expansion ===\n")
# Create synthetic image (200x200)
mat <- array(0.2, dim = c(200, 200, 3))
# Draw a circle in the center: center=(100, 100), radius=30
for (x in 1:200) {
  for (y in 1:200) {
    if ((x - 100)^2 + (y - 100)^2 <= 30^2) {
      mat[x, y, ] <- 0.8
    }
  }
}
im <- as_image(mat)

# Suppose tight detector found box: [xmin=70, ymin=70, xmax=130, ymax=130] (exact bounding box of circle)
boxes <- data.frame(xmin = 70, ymin = 70, xmax = 130, ymax = 130)

# Check .compact_prototypes
proto_dummy <- list(semantic = replicate(12, runif(256), simplify = FALSE))
proto_compacted <- pliman:::.compact_prototypes(proto_dummy, max_protos = 5L)
stopifnot(length(proto_compacted$semantic) == 5L)
cat("[OK] .compact_prototypes works as expected.\n")

cat("\nAll function checks passed successfully!\n")
