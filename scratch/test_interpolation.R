cat("=== Testing Box and Mask Interpolation in Tracker ===\n")
devtools::load_all(".")

# 1. Test interpolate_box = TRUE (default)
cat("\n--- Test 1: interpolate_box = TRUE ---\n")
state <- list()

# Frame 1: initial detection
res1 <- update_tracker_cpp(
  xmin = 100, ymin = 100, xmax = 200, ymax = 200,
  labels = "seed", scores = 0.95, state = state,
  class_ids = 0L, mask_ids = 1L,
  max_lost = 3L
)
state <- res1$state
cat(sprintf("Frame 1: Track ID = %d, lost = %d, boxes = %d\n", res1$track_ids[1], res1$lost[1], length(res1$track_ids)))

# Frame 2: motion towards right-down (vx = 10, vy = 10)
res2 <- update_tracker_cpp(
  xmin = 110, ymin = 110, xmax = 210, ymax = 210,
  labels = "seed", scores = 0.95, state = state,
  class_ids = 0L, mask_ids = 1L,
  max_lost = 3L
)
state <- res2$state
cat(sprintf("Frame 2: Track ID = %d, lost = %d, boxes = %d\n", res2$track_ids[1], res2$lost[1], length(res2$track_ids)))

# Frame 3: object lost (0 detections) with interpolate_box = TRUE
res3 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state,
  max_lost = 3L, interpolate_box = TRUE, interpolate_mask = TRUE
)
state <- res3$state
cat(sprintf("Frame 3 (lost 1): boxes = %d, Track ID = %d, lost = %d, xmin = %.1f, shift_x = %.1f\n",
            length(res3$track_ids), res3$track_ids[1], res3$lost[1], res3$xmin[1], res3$mask_shift_x[1]))
stopifnot(length(res3$track_ids) == 1)
stopifnot(res3$lost[1] == 1)
stopifnot(res3$xmin[1] > 110) # Position propagated forward

# Frame 4: object lost again (lost = 2)
res4 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state,
  max_lost = 3L, interpolate_box = TRUE
)
state <- res4$state
cat(sprintf("Frame 4 (lost 2): boxes = %d, lost = %d, xmin = %.1f, shift_x = %.1f\n",
            length(res4$track_ids), res4$lost[1], res4$xmin[1], res4$mask_shift_x[1]))

# Frame 5: object lost again (lost = 3, max_lost reached)
res5 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state,
  max_lost = 3L, interpolate_box = TRUE
)
state <- res5$state
cat(sprintf("Frame 5 (lost 3 = max_lost): boxes = %d, lost = %d\n",
            length(res5$track_ids), res5$lost[1]))

# Frame 6: lost > max_lost -> object disappears!
res6 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state,
  max_lost = 3L, interpolate_box = TRUE
)
state <- res6$state
cat(sprintf("Frame 6 (lost > max_lost): boxes = %d (expected 0), surviving tracks in state = %d (expected 0)\n",
            length(res6$track_ids), length(state$tracks)))
stopifnot(length(res6$track_ids) == 0)
stopifnot(length(state$tracks) == 0)


# 2. Test interpolate_box = FALSE
cat("\n--- Test 2: interpolate_box = FALSE ---\n")
state2 <- list()
res2_1 <- update_tracker_cpp(
  xmin = 100, ymin = 100, xmax = 200, ymax = 200,
  labels = "seed", scores = 0.95, state = state2,
  class_ids = 0L, mask_ids = 1L, max_lost = 3L
)
state2 <- res2_1$state

res2_2 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state2,
  max_lost = 3L, interpolate_box = FALSE
)
state2 <- res2_2$state
cat(sprintf("Frame 2 (lost 1, interpolate_box=FALSE): boxes returned = %d (expected 0)\n", length(res2_2$track_ids)))
cat(sprintf("Surviving tracks in state = %d (expected 1, object kept in tracker)\n", length(state2$tracks)))
stopifnot(length(res2_2$track_ids) == 0)
stopifnot(length(state2$tracks) == 1)

# Frame 3: Detection returns - track ID preserved!
res2_3 <- update_tracker_cpp(
  xmin = 105, ymin = 105, xmax = 205, ymax = 205,
  labels = "seed", scores = 0.95, state = state2,
  class_ids = 0L, mask_ids = 1L, max_lost = 3L
)
cat(sprintf("Frame 3 (detection returns): Track ID = %d (expected %d), lost = %d\n",
            res2_3$track_ids[1], res2_1$track_ids[1], res2_3$lost[1]))
stopifnot(res2_3$track_ids[1] == res2_1$track_ids[1])
stopifnot(res2_3$lost[1] == 0)


# 3. Test interpolate_mask = FALSE
cat("\n--- Test 3: interpolate_mask = FALSE ---\n")
state3 <- list()
res3_1 <- update_tracker_cpp(
  xmin = 100, ymin = 100, xmax = 200, ymax = 200,
  labels = "seed", scores = 0.95, state = state3,
  class_ids = 0L, mask_ids = 1L, max_lost = 3L
)
state3 <- res3_1$state

res3_2 <- update_tracker_cpp(
  xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
  labels = character(0), scores = numeric(0), state = state3,
  max_lost = 3L, interpolate_box = TRUE, interpolate_mask = FALSE
)
cat(sprintf("Frame 2 (lost 1, interpolate_box=TRUE, interpolate_mask=FALSE): boxes = %d, mask_ids = %d (expected -1)\n",
            length(res3_2$track_ids), res3_2$mask_ids[1]))
stopifnot(length(res3_2$track_ids) == 1)
stopifnot(res3_2$mask_ids[1] == -1)

cat("\nALL INTERPOLATION TESTS PASSED SUCCESSFULLY!\n")
