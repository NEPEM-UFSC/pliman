
pkgload::load_all("d:/Desktop/pliman")

# Create synthetic image (200x200)
mat <- array(0.2, dim = c(200, 200, 3))
for (x in 1:200) {
  for (y in 1:200) {
    if ((x - 100)^2 + (y - 100)^2 <= 30^2) {
      mat[x, y, ] <- 0.8
    }
  }
}
im <- as_image(mat)

# Bounding box
boxes <- data.frame(xmin = 70, ymin = 70, xmax = 130, ymax = 130)

cat("Calling .sam_segment_boxes with box_margin = 0.05...\n")
res <- pliman:::.sam_segment_boxes(im, boxes = boxes, box_margin = 0.05, engine = "cpu", verbose = FALSE)

cat("Extracted contours length:", length(res$contours), "\n")
if (length(res$contours) > 0) {
  cnt <- res$contours[[1]]
  cat(sprintf("Contour points: %d, X range: [%.1f, %.1f], Y range: [%.1f, %.1f]\n",
              nrow(cnt), min(cnt[, 1]), max(cnt[, 1]), min(cnt[, 2]), max(cnt[, 2])))
}
cat("[OK] SAM 2.1 segment boxes executed successfully with box_margin!\n")
