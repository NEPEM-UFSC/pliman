# ==============================================================================
# FEW-SHOT EXEMPLAR ANNOTATION & AUTO-LABELING FOR YOLO DATASETS
# pliman: Plant Image Analysis in R
# ==============================================================================


.point_in_polygon <- function(px, py, pol_x, pol_y) {
  n <- length(pol_x)
  if (n < 3L) return(FALSE)
  inside <- FALSE
  j <- n
  for (i in seq_len(n)) {
    xi <- pol_x[i]; yi <- pol_y[i]
    xj <- pol_x[j]; yj <- pol_y[j]
    if (((yi > py) != (yj > py)) &&
        (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) {
      inside <- !inside
    }
    j <- i
  }
  inside
}

.render_fewshot_review_plot <- function(im,
                                        boxes = data.frame(),
                                        polygons = list(),
                                        exemplar_points = NULL,
                                        task = "segment",
                                        label = "object",
                                        conf_threshold = 0.5,
                                        rainbow = TRUE,
                                        show_class = TRUE,
                                        show_conf = TRUE,
                                        badge = TRUE,
                                        img_name = "",
                                        img_idx = 1L,
                                        total_imgs = 1L) {
  w <- dim(im)[1]
  h <- dim(im)[2]
  n_b <- if (is.data.frame(boxes)) nrow(boxes) else 0L
  n_p <- if (is.list(polygons)) length(polygons) else 0L
  n_inst <- if (task == "segment" && n_p > 0L) n_p else n_b

  plot(im)

  if (isTRUE(rainbow) && n_inst > 0L) {
    cols_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
    cols_poly <- grDevices::adjustcolor(cols_border, alpha.f = 0.40)
  } else {
    cols_border <- rep("#00FFCC", max(1L, n_inst))
    cols_poly <- rep("#00FFCC55", max(1L, n_inst))
  }

  # 1. Draw segmentation polygons
  if (identical(task, "segment") && n_p > 0L) {
    for (i in seq_along(polygons)) {
      poly <- polygons[[i]]
      if (!is.null(poly) && length(poly) >= 6L) {
        px <- poly[c(TRUE, FALSE)] * w
        py <- poly[c(FALSE, TRUE)] * h
        cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
        cp <- cols_poly[((i - 1L) %% length(cols_poly)) + 1L]
        graphics::polygon(px, py, col = cp, border = cb, lwd = 2)
      }
    }
  }

  # 2. Draw bounding boxes (if detect, or if segment has boxes and no polygons)
  if (identical(task, "detect") && n_b > 0L) {
    for (b in seq_len(n_b)) {
      bx <- boxes[b, ]
      cb <- cols_border[((b - 1L) %% length(cols_border)) + 1L]
      graphics::rect(bx$xmin, bx$ymin, bx$xmax, bx$ymax, border = cb, lwd = 2)
    }
  }

  # 3. Draw labels and badges (for detect, or for segment when boxes are available)
  if (n_b > 0L && (isTRUE(show_class) || isTRUE(show_conf))) {
    for (b in seq_len(n_b)) {
      bx <- boxes[b, ]
      cb <- cols_border[((b - 1L) %% length(cols_border)) + 1L]
      lbl_parts <- character(0)
      if (isTRUE(show_class)) {
        lbl_txt <- if (!is.null(bx$class_name) && nzchar(bx$class_name)) bx$class_name else label
        if (nzchar(lbl_txt)) lbl_parts <- c(lbl_parts, lbl_txt)
      }
      if (isTRUE(show_conf) && !is.null(bx$conf) && !is.na(bx$conf)) {
        lbl_parts <- c(lbl_parts, sprintf("%.2f", bx$conf))
      }
      tag <- paste(lbl_parts, collapse = " ")

      if (nzchar(tag)) {
        tw <- graphics::strwidth(tag, cex = 0.7)
        th <- graphics::strheight(tag, cex = 0.7)
        y_top <- max(th + 6, bx$ymin)
        if (isTRUE(badge)) {
          bg_tag <- grDevices::adjustcolor(cb, alpha.f = 0.90)
          graphics::rect(bx$xmin, y_top - th - 6, bx$xmin + tw + 6, y_top, col = bg_tag, border = NA)
          graphics::text(bx$xmin + 3, y_top - 3, labels = tag, col = "#FFFFFF", cex = 0.7, font = 2, adj = c(0, 1))
        } else {
          graphics::text(bx$xmin + 3, y_top - 3, labels = tag, col = cb, cex = 0.7, font = 2, adj = c(0, 1))
        }
      }
    }
  }

  # 4. Draw exemplar points as small yellow dots
  if (!is.null(exemplar_points) && length(exemplar_points) > 0L) {
    if (is.matrix(exemplar_points) || is.data.frame(exemplar_points)) {
      graphics::points(exemplar_points[, 1], exemplar_points[, 2], col = "#FFFF00", pch = 19, cex = 1.3)
      graphics::points(exemplar_points[, 1], exemplar_points[, 2], col = "#000000", pch = 1, cex = 1.3, lwd = 1.5)
    }
  }

  # Title with metadata
  title_str <- sprintf("Image [%d/%d] (%s): %d %s [%s, conf >= %.2f]",
                       img_idx, total_imgs, img_name, n_inst, label, task, conf_threshold)
  graphics::title(main = title_str, font.main = 2, col.main = "#006622", cex.main = 1.05)
}

.interactive_fewshot_review_step <- function(curr_im,
                                             r_boxes,
                                             r_polys,
                                             curr_points,
                                             c_num_arr,
                                             cw,
                                             ch,
                                             pre_proto,
                                             latest_proto = NULL,
                                             task = "segment",
                                             label = "object",
                                             conf_threshold = 0.50,
                                             radius = NULL,
                                             min_dist = 12,
                                             max_objects = NULL,
                                             method = "sam",
                                             engine = "gpu",
                                             device_id = -1,
                                             rainbow = TRUE,
                                             show_class = TRUE,
                                             show_conf = TRUE,
                                             badge = TRUE,
                                             img_name = "",
                                             img_idx = 1L,
                                             total_imgs = 1L) {
  if (is.null(pre_proto) && !is.null(latest_proto)) {
    pre_proto <- latest_proto
  }

  .render_fewshot_review_plot(
    im = curr_im,
    boxes = r_boxes,
    polygons = r_polys,
    exemplar_points = curr_points,
    task = task,
    label = label,
    conf_threshold = conf_threshold,
    rainbow = rainbow,
    show_class = show_class,
    show_conf = show_conf,
    badge = badge,
    img_name = img_name,
    img_idx = img_idx,
    total_imgs = total_imgs
  )

  deleted_centers <- list()
  curr_conf <- conf_threshold

  repeat {
    n_inst <- if (identical(task, "segment") && length(r_polys) > 0L) length(r_polys) else nrow(r_boxes)
    cat("\n")
    cli::cli_rule(left = paste0("Few-Shot Review [", img_name, "]"))
    cli::cli_alert_success("{n_inst} instance(s) detected [task: {.val {task}}, conf >= {curr_conf}].")
    cli::cli_bullets(c(
      "*" = "{.bold [Enter] / c}: Accept detections and proceed to next image",
      "*" = "{.bold r}: Remove false positives (click on masks/boxes to delete, <Esc> when done)",
      "*" = "{.bold a}: Add missed objects (click on objects, <Esc> when done)",
      "*" = "{.bold t <val>}: Adjust confidence threshold (e.g. {.code t 0.70} or {.code t 0.50})",
      "*" = "{.bold s}: Sample again (clear exemplars and re-click points for this image)",
      "*" = "{.bold q}: Quit / finish model fitting now with completed images"
    ))

    ans <- tolower(trimws(readline("Select action [Enter to accept]: ")))

    if (ans == "" || ans %in% c("c", "continue", "y", "yes", "ok")) {
      return(list(
        action = "continue",
        boxes = r_boxes,
        polygons = r_polys,
        points = curr_points,
        conf = curr_conf,
        prototypes = pre_proto
      ))
    }

    if (ans %in% c("s", "sample", "resample", "retry")) {
      return(list(action = "resample"))
    }

    if (ans %in% c("q", "quit", "stop", "abort")) {
      return(list(
        action = "stop",
        boxes = r_boxes,
        polygons = r_polys,
        points = curr_points,
        conf = curr_conf,
        prototypes = pre_proto
      ))
    }

    if (ans %in% c("r", "remove", "delete", "d")) {
      cli::cli_alert_info("Click on masks or bounding boxes to remove them.")
      cli::cli_alert_info("Press <Esc> or right-click when finished.")

      del_pts <- tryCatch(
        graphics::locator(n = 200, type = "p", pch = 4, col = "#FF0033", cex = 1.6, lwd = 2.5),
        error = function(e) NULL
      )

      if (!is.null(del_pts) && length(del_pts$x) > 0L) {
        remove_idx <- integer()
        for (pt_i in seq_along(del_pts$x)) {
          cx <- del_pts$x[pt_i]; cy <- del_pts$y[pt_i]
          best_match <- NULL
          best_dist <- Inf

          if (identical(task, "segment") && length(r_polys) > 0L) {
            for (pi in seq_along(r_polys)) {
              if (pi %in% remove_idx) next
              poly <- r_polys[[pi]]
              if (is.null(poly) || length(poly) < 6L) next
              px <- poly[c(TRUE, FALSE)] * cw
              py <- poly[c(FALSE, TRUE)] * ch
              cent_x <- mean(px); cent_y <- mean(py)
              d_cent <- (cent_x - cx)^2 + (cent_y - cy)^2

              inside <- .point_in_polygon(cx, cy, px, py)
              if (inside) {
                if (d_cent < best_dist) {
                  best_dist <- d_cent
                  best_match <- pi
                }
              } else if (d_cent < (min(cw, ch) * 0.08)^2 && d_cent < best_dist) {
                best_dist <- d_cent
                best_match <- pi
              }
            }
          } else if (nrow(r_boxes) > 0L) {
            for (bi in seq_len(nrow(r_boxes))) {
              if (bi %in% remove_idx) next
              bx <- r_boxes[bi, ]
              cent_x <- (bx$xmin + bx$xmax) / 2
              cent_y <- (bx$ymin + bx$ymax) / 2
              d_cent <- (cent_x - cx)^2 + (cent_y - cy)^2

              inside <- (cx >= bx$xmin && cx <= bx$xmax && cy >= bx$ymin && cy <= bx$ymax)
              if (inside) {
                if (d_cent < best_dist) {
                  best_dist <- d_cent
                  best_match <- bi
                }
              } else if (d_cent < (min(cw, ch) * 0.08)^2 && d_cent < best_dist) {
                best_dist <- d_cent
                best_match <- bi
              }
            }
          }

          if (!is.null(best_match)) {
            remove_idx <- c(remove_idx, best_match)
          }
        }

        remove_idx <- unique(remove_idx)
        if (length(remove_idx) > 0L) {
          if (nrow(r_boxes) > 0L) {
            for (ri in remove_idx) {
              if (ri <= nrow(r_boxes)) {
                deleted_centers[[length(deleted_centers) + 1L]] <- c(
                  (r_boxes$xmin[ri] + r_boxes$xmax[ri]) / 2,
                  (r_boxes$ymin[ri] + r_boxes$ymax[ri]) / 2
                )
              }
            }
            r_boxes <- r_boxes[-remove_idx, , drop = FALSE]
          }
          if (length(r_polys) > 0L) {
            r_polys <- r_polys[-remove_idx]
          }

          cli::cli_alert_success("Removed {length(remove_idx)} instance(s). {if (identical(task, 'segment') && length(r_polys) > 0L) length(r_polys) else nrow(r_boxes)} remaining.")
          .render_fewshot_review_plot(
            im = curr_im,
            boxes = r_boxes,
            polygons = r_polys,
            exemplar_points = curr_points,
            task = task,
            label = label,
            conf_threshold = curr_conf,
            rainbow = rainbow,
            show_class = show_class,
            show_conf = show_conf,
            badge = badge,
            img_name = img_name,
            img_idx = img_idx,
            total_imgs = total_imgs
          )
        } else {
          cli::cli_alert_warning("No instance was close to the clicked point(s).")
        }
      }
      next
    }

    if (ans %in% c("a", "add")) {
      cli::cli_alert_info("Click on missed objects to add them.")
      cli::cli_alert_info("Press <Esc> or right-click when finished.")

      add_pts <- tryCatch(
        graphics::locator(n = 64, type = "p", pch = 3, col = "#00FF66", cex = 1.4, lwd = 2.5),
        error = function(e) NULL
      )

      if (!is.null(add_pts) && length(add_pts$x) > 0L) {
        new_coords <- cbind(x = add_pts$x, y = add_pts$y)
        curr_points <- if (!is.null(curr_points)) rbind(curr_points, new_coords) else new_coords

        cli::cli_progress_step(
          msg = "Re-evaluating SAM 2.1 with {nrow(new_coords)} added exemplar(s) [{toupper(engine)}]... ",
          msg_done = "SAM 2.1 re-evaluation complete"
        )
        new_pres <- tryCatch({
          .run_persam(
            mat = c_num_arr,
            exemplar_points = new_coords,
            precomputed_prototypes = pre_proto,
            sim_threshold = curr_conf,
            min_dist = min_dist,
            max_objects = max_objects,
            feat_res = 256L,
            engine = engine,
            device_id = device_id,
            fill_hull = TRUE,
            verbose = FALSE
          )
        }, error = function(e) NULL)

        cli::cli_progress_done()

        if (!is.null(new_pres) && !is.null(new_pres$prototypes)) {
          if (!is.null(pre_proto)) {
            pre_proto$semantic <- c(pre_proto$semantic, new_pres$prototypes$semantic)
            pre_proto$hr0 <- c(pre_proto$hr0, new_pres$prototypes$hr0)
            if (!is.null(new_pres$prototypes$avg_area) && new_pres$prototypes$avg_area > 0) {
              pre_proto$avg_area <- mean(c(pre_proto$avg_area, new_pres$prototypes$avg_area))
            }
          } else {
            pre_proto <- new_pres$prototypes
          }
        }

        nbx <- data.frame(
          xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
          conf = numeric(0), class_id = integer(0), class_name = character(0),
          stringsAsFactors = FALSE
        )
        npl <- list()

        if (!is.null(new_pres) && nrow(new_pres$boxes) > 0L) {
          nbx <- data.frame(
            xmin = as.numeric(new_pres$boxes$xmin),
            ymin = as.numeric(new_pres$boxes$ymin),
            xmax = as.numeric(new_pres$boxes$xmax),
            ymax = as.numeric(new_pres$boxes$ymax),
            conf = round(as.numeric(new_pres$boxes$score), 3),
            class_id = 0L,
            class_name = label,
            stringsAsFactors = FALSE
          )
          for (ci in seq_along(new_pres$contours)) {
            cnt <- new_pres$contours[[ci]]
            if (is.matrix(cnt) && nrow(cnt) >= 3L) {
              gx <- cnt[, 1] / cw; gy <- cnt[, 2] / ch
              if (length(gx) > 32L) {
                idx <- round(seq(1, length(gx), length.out = 32))
                gx <- gx[idx]; gy <- gy[idx]
              }
              npl[[ci]] <- as.vector(rbind(gx, gy))
            }
          }
        }

        # Seamlessly merge previous r_boxes with newly detected boxes
        if (nrow(r_boxes) > 0L) {
          keep_old <- rep(TRUE, nrow(r_boxes))
          for (oi in seq_len(nrow(r_boxes))) {
            ocx <- (r_boxes$xmin[oi] + r_boxes$xmax[oi]) / 2
            ocy <- (r_boxes$ymin[oi] + r_boxes$ymax[oi]) / 2
            for (ni in seq_len(nrow(nbx))) {
              ncx <- (nbx$xmin[ni] + nbx$xmax[ni]) / 2
              ncy <- (nbx$ymin[ni] + nbx$ymax[ni]) / 2
              if (((ocx - ncx)^2 + (ocy - ncy)^2) < (min(cw, ch) * 0.04)^2) {
                keep_old[oi] <- FALSE
                break
              }
            }
          }
          if (any(keep_old)) {
            nbx <- rbind(r_boxes[keep_old, , drop = FALSE], nbx)
            if (length(r_polys) > 0L) {
              npl <- c(r_polys[keep_old[seq_along(r_polys)]], npl)
            }
          }
        }

        # Fallback guarantee: ensure EVERY clicked point has a box on screen!
        def_w <- if (nrow(nbx) > 0L) stats::median(nbx$xmax - nbx$xmin) else min(cw, ch) * 0.04
        def_h <- if (nrow(nbx) > 0L) stats::median(nbx$ymax - nbx$ymin) else min(cw, ch) * 0.04
        if (is.na(def_w) || def_w <= 0) def_w <- min(cw, ch) * 0.04
        if (is.na(def_h) || def_h <= 0) def_h <- min(cw, ch) * 0.04

        for (pi in seq_len(nrow(new_coords))) {
          px <- new_coords[pi, 1]; py <- new_coords[pi, 2]
          has_box <- any(px >= nbx$xmin & px <= nbx$xmax & py >= nbx$ymin & py <= nbx$ymax)
          if (!has_box) {
            fallback_row <- data.frame(
              xmin = max(0, px - def_w / 2),
              ymin = max(0, py - def_h / 2),
              xmax = min(cw, px + def_w / 2),
              ymax = min(ch, py + def_h / 2),
              conf = 1.0,
              class_id = 0L,
              class_name = label,
              stringsAsFactors = FALSE
            )
            nbx <- rbind(nbx, fallback_row)
            if (identical(task, "segment")) {
              th <- seq(0, 2 * pi, length.out = 16)
              circ_x <- (px + (def_w / 2) * cos(th)) / cw
              circ_y <- (py + (def_h / 2) * sin(th)) / ch
              npl[[length(npl) + 1L]] <- as.vector(rbind(circ_x, circ_y))
            }
          }
        }

        # Filter any previously deleted centers
        if (length(deleted_centers) > 0L && nrow(nbx) > 0L) {
          keep_idx <- rep(TRUE, nrow(nbx))
          for (bi in seq_len(nrow(nbx))) {
            bc_x <- (nbx$xmin[bi] + nbx$xmax[bi]) / 2
            bc_y <- (nbx$ymin[bi] + nbx$ymax[bi]) / 2
            for (dc in deleted_centers) {
              if (((bc_x - dc[1])^2 + (bc_y - dc[2])^2) < (min(cw, ch) * 0.04)^2) {
                keep_idx[bi] <- FALSE
                break
              }
            }
          }
          nbx <- nbx[keep_idx, , drop = FALSE]
          if (length(npl) > 0L) {
            npl <- npl[keep_idx[seq_along(npl)]]
          }
        }

        r_boxes <- nbx
        r_polys <- npl

        cli::cli_alert_success("Model updated! {if (identical(task, 'segment') && length(r_polys) > 0L) length(r_polys) else nrow(r_boxes)} total instance(s) detected.")
        .render_fewshot_review_plot(
          im = curr_im,
          boxes = r_boxes,
          polygons = r_polys,
          exemplar_points = curr_points,
          task = task,
          label = label,
          conf_threshold = curr_conf,
          rainbow = rainbow,
          show_class = show_class,
          show_conf = show_conf,
          badge = badge,
          img_name = img_name,
          img_idx = img_idx,
          total_imgs = total_imgs
        )
      }
      next
    }

    if (startsWith(ans, "t")) {
      new_val_str <- trimws(sub("^t\\s*", "", ans))
      new_val <- if (nzchar(new_val_str)) as.numeric(new_val_str) else NA_real_
      if (is.na(new_val)) {
        ans_t <- readline("Enter new confidence threshold (e.g. 0.70): ")
        new_val <- as.numeric(trimws(ans_t))
      }

      if (!is.na(new_val) && new_val > 0 && new_val < 1) {
        curr_conf <- new_val
        cli::cli_progress_step(
          msg = "Re-evaluating detections with threshold = {curr_conf} [{toupper(engine)}]... ",
          msg_done = "Threshold update complete"
        )
        new_pres <- tryCatch({
          .run_persam(
            mat = c_num_arr,
            exemplar_points = curr_points,
            precomputed_prototypes = pre_proto,
            sim_threshold = curr_conf,
            min_dist = min_dist,
            max_objects = max_objects,
            feat_res = 256L,
            engine = engine,
            device_id = device_id,
            fill_hull = TRUE,
            verbose = FALSE
          )
        }, error = function(e) NULL)

        cli::cli_progress_done()

        if (!is.null(new_pres) && nrow(new_pres$boxes) > 0L) {
          nbx <- data.frame(
            xmin = as.numeric(new_pres$boxes$xmin),
            ymin = as.numeric(new_pres$boxes$ymin),
            xmax = as.numeric(new_pres$boxes$xmax),
            ymax = as.numeric(new_pres$boxes$ymax),
            conf = round(as.numeric(new_pres$boxes$score), 3),
            class_id = 0L,
            class_name = label,
            stringsAsFactors = FALSE
          )
          npl <- list()
          for (ci in seq_along(new_pres$contours)) {
            cnt <- new_pres$contours[[ci]]
            if (is.matrix(cnt) && nrow(cnt) >= 3L) {
              gx <- cnt[, 1] / cw; gy <- cnt[, 2] / ch
              if (length(gx) > 32L) {
                idx <- round(seq(1, length(gx), length.out = 32))
                gx <- gx[idx]; gy <- gy[idx]
              }
              npl[[ci]] <- as.vector(rbind(gx, gy))
            }
          }

          if (length(deleted_centers) > 0L) {
            keep_idx <- rep(TRUE, nrow(nbx))
            for (bi in seq_len(nrow(nbx))) {
              bc_x <- (nbx$xmin[bi] + nbx$xmax[bi]) / 2
              bc_y <- (nbx$ymin[bi] + nbx$ymax[bi]) / 2
              for (dc in deleted_centers) {
                if (((bc_x - dc[1])^2 + (bc_y - dc[2])^2) < (min(cw, ch) * 0.05)^2) {
                  keep_idx[bi] <- FALSE
                  break
                }
              }
            }
            nbx <- nbx[keep_idx, , drop = FALSE]
            if (length(npl) > 0L) {
              npl <- npl[keep_idx[seq_along(npl)]]
            }
          }

          r_boxes <- nbx
          r_polys <- npl
        } else {
          r_boxes <- data.frame()
          r_polys <- list()
        }

        cli::cli_alert_success("Updated threshold to {curr_conf}. {if (identical(task, 'segment') && length(r_polys) > 0L) length(r_polys) else nrow(r_boxes)} instance(s) detected.")
        .render_fewshot_review_plot(
          im = curr_im,
          boxes = r_boxes,
          polygons = r_polys,
          exemplar_points = curr_points,
          task = task,
          label = label,
          conf_threshold = curr_conf,
          rainbow = rainbow,
          show_class = show_class,
          show_conf = show_conf,
          badge = badge,
          img_name = img_name,
          img_idx = img_idx,
          total_imgs = total_imgs
        )
      } else {
        cli::cli_alert_warning("Invalid threshold value. Must be between 0 and 1 (e.g. 't 0.70').")
      }
      next
    }

    cli::cli_alert_warning("Unrecognized option: {.val {ans}}. Please select [Enter], 'r', 'a', 't', 's', or 'q'.")
  }
}

#' @title Interactive Few-Shot Exemplar Fitting for Object Detection and Segmentation
#' @name yolo_fewshot_fit
#' @description
#' `yolo_fewshot_fit()` implements an interactive few-shot visual learning engine powered by
#' Meta AI's Segment Anything 2.1 (SAM 2.1 / PerSAM). The user selects 3–5 sample instances
#' of an object (e.g. seeds, beans, logs, leaves, lesions, fruits) by clicking on an image.
#'
#' Using the PerSAM exemplar matching engine, the function computes visual prototypes from the
#' exemplars and segments all other matching instances across the entire image with high precision.
#'
#' The user can specify `task = "segment"` to obtain instance masks and polygonal boundaries
#' for YOLO segmentation datasets, or `task = "detect"` for oriented/axis-aligned bounding boxes.
#'
#' The resulting model object (`yolo_fewshot_model`) can be inspected, plotted, and directly
#' propagated to hundreds of unannotated images in a folder using [yolo_dataset_autolabel()].
#'
#' @param img An `image` object, file path, or matrix.
#' @param points,sample Optional $(x, y)$ coordinates of exemplar points. If `NULL` (default),
#'   an interactive plotting window opens for the user to click on 3–5 representative objects.
#' @param label Character vector naming the target class(es). Can be a single class (e.g., `"bean"`, `"leaf"`, `"seed"`, `"tora"`)
#'   or a character vector of multiple classes (e.g., `c("leaf", "lesion")`, `c("healthy", "damaged")`). Defaults to `"object"`.
#'   When multiple classes are specified, `yolo_fewshot_fit()` interactively collects exemplars for each class in sequence,
#'   or uses the respective elements from `points` if a list is provided.
#' @param task Target task: `"segment"` (instance segmentation polygons and masks) or
#'   `"detect"` (bounding boxes). Defaults to `"segment"`.
#' @param conf_threshold Numeric scalar between 0 and 1 specifying detection confidence/similarity threshold. Defaults to `0.50` (or `0.25` for adaptive).
#' @param radius Optional estimated radius in pixels for target objects. If `NULL` (default),
#'   automatically detected from exemplar edge gradients.
#' @param min_dist Minimum distance in pixels between detected instances. Defaults to `12`.
#' @param min_area Minimum object area in pixels. If `NULL` (default), automatically calibrated
#'   from the exemplar scale (`0.20 * median_exemplar_area`) to filter out background noise.
#' @param max_area Maximum object area in pixels. If `NULL` (default), automatically calibrated
#'   (`3.80 * median_exemplar_area`).
#' @param max_objects Maximum number of instances to detect. Defaults to `NULL` (no limit).
#' @param method Segmentation method: `"sam"` (or `"persam"`, Meta Segment Anything 2.1 via ONNX)
#'   or `"adaptive"` (fast color-texture morphology). Defaults to `"sam"`.
#' @param engine Execution engine for SAM: `"cpu"` or `"gpu"`. Defaults to `"cpu"`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param interactive Logical. Whether to enable interactive review and refinement after exemplar segmentation (allows adding missed objects, removing false positives, adjusting threshold, or resampling). Defaults to `interactive()`.
#' @param rainbow Logical. If `TRUE` (default), plots instances in distinct rainbow colors.
#' @param show_class Logical. Whether to display class name on labels/badges. Defaults to `TRUE`.
#' @param show_conf Logical. Whether to display confidence score on labels/badges. Defaults to `TRUE`.
#' @param badge Logical. Whether to draw a filled background badge behind labels. Defaults to `TRUE`.
#' @param plot Logical. Display the detection and segmentation overlay. Defaults to `TRUE`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#' @param ... Additional arguments passed to downstream methods or aliases (e.g., `show_labels`, `show_scores`, `label_box`).
#'
#' @return An object of class `yolo_fewshot_model` containing exemplar signatures,
#'   detected instances on the reference image, and parameters for batch propagation.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # 1. Fit exemplar model on reference photo (clicks or supplied coordinates)
#' modelo_tora <- yolo_fewshot_fit("toras/toras22.png", label = "tora", task = "segment", method = "sam")
#'
#' # 2. Inspect segmentation overlay on reference photo
#' plot(modelo_tora)
#'
#' # 3. Auto-label entire directory of 300 photos into a YOLO dataset
#' ds <- yolo_dataset_autolabel("toras_raw/", model = modelo_tora, output_dir = "dataset_toras_yolo")
#' }
yolo_fewshot_fit <- function(img,
                             points = NULL,
                             sample = NULL,
                             label = "object",
                             task = c("segment", "detect"),
                             conf_threshold = 0.50,
                             radius = NULL,
                             min_dist = 12,
                             min_area = NULL,
                             max_area = NULL,
                             max_objects = NULL,
                             method = c("sam", "persam", "adaptive"),
                             engine = c("gpu", "cpu"),
                             device_id = -1,
                             rainbow = TRUE,
                             show_class = TRUE,
                             show_conf = TRUE,
                             badge = TRUE,
                             interactive = base::interactive(),
                             plot = TRUE,
                             verbose = TRUE,
                             ...) {
  dots <- list(...)
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  task <- match.arg(task)
  method <- match.arg(method)
  engine <- match.arg(engine)

  # Support multi-class few-shot fitting
  if (length(label) > 1L) {
    if (isTRUE(verbose)) {
      cli::cli_h1("Multi-Class Few-Shot Training (pliman)")
      cli::cli_alert_info("Fitting {length(label)} visual classes: {.val {label}}")
    }
    models <- vector("list", length(label))
    names(models) <- label
    for (i in seq_along(label)) {
      lbl <- label[i]
      if (isTRUE(verbose)) {
        cli::cli_h2("Class {i}/{length(label)}: {.val {lbl}}")
      }
      pts_i <- if (is.list(points) && lbl %in% names(points)) {
        points[[lbl]]
      } else if (is.list(points) && length(points) >= i && !is.data.frame(points)) {
        points[[i]]
      } else if (i == 1L && (is.matrix(points) || is.data.frame(points))) {
        points
      } else {
        NULL
      }
      mod_i <- yolo_fewshot_fit(
        img = img,
        points = pts_i,
        sample = sample,
        label = lbl,
        task = task,
        conf_threshold = conf_threshold,
        radius = radius,
        min_dist = min_dist,
        min_area = min_area,
        max_area = max_area,
        max_objects = max_objects,
        method = method,
        engine = engine,
        device_id = device_id,
        rainbow = rainbow,
        show_class = show_class,
        show_conf = show_conf,
        badge = badge,
        interactive = interactive,
        plot = FALSE,
        verbose = verbose
      )
      mod_i$class_id <- as.integer(i - 1L)
      if (!is.null(mod_i$detected_boxes) && nrow(mod_i$detected_boxes) > 0L) {
        mod_i$detected_boxes$class_id <- as.integer(i - 1L)
      }
      models[[i]] <- mod_i
    }

    all_ref_b <- list()
    all_ref_p <- list()
    for (m_idx in seq_along(models)) {
      if (!is.null(models[[m_idx]]$detected_boxes) && nrow(models[[m_idx]]$detected_boxes) > 0L) {
        all_ref_b[[length(all_ref_b) + 1L]] <- models[[m_idx]]$detected_boxes
        all_ref_p <- c(all_ref_p, models[[m_idx]]$detected_polygons)
      }
    }
    comb_ref_boxes <- if (length(all_ref_b) > 0L) do.call(rbind, all_ref_b) else data.frame()
    if (nrow(comb_ref_boxes) > 0L) {
      nms_ref <- .cross_class_nms(comb_ref_boxes, all_ref_p, iou_thresh = 0.40)
      comb_ref_boxes <- nms_ref$boxes
      comb_ref_polys <- nms_ref$polygons
    } else {
      comb_ref_polys <- list()
    }

    comb_obj <- structure(
      list(
        label = label,
        classes = label,
        models = models,
        task = task,
        conf_threshold = conf_threshold,
        method = method,
        engine = engine,
        rainbow = rainbow,
        show_class = show_class,
        show_conf = show_conf,
        badge = badge,
        detected_boxes = comb_ref_boxes,
        detected_polygons = comb_ref_polys,
        reference_image = models[[1]]$reference_image,
        reference_images = models[[1]]$reference_images
      ),
      class = c("yolo_fewshot_multimodel", "yolo_fewshot_model")
    )

    if (isTRUE(verbose)) {
      cli::cli_alert_success("Multi-class few-shot model successfully fitted for classes: {.val {label}}!")
    }

    if (isTRUE(plot)) {
      plot(comb_obj, rainbow = rainbow, show_class = show_class, show_conf = show_conf, badge = badge)
    }

    return(invisible(comb_obj))
  }

  # Resolve input image(s): directory with sample, vector of file paths, or single image
  img_list <- list()
  if (is.character(img) && length(img) == 1L && dir.exists(img[1])) {
    img_files <- list.files(img[1], pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
    if (length(img_files) == 0L) {
      cli::cli_abort("No image files found in directory {.path {img[1]}}.")
    }
    n_s <- if (!is.null(sample) && sample > 0L) min(length(img_files), as.integer(sample)) else min(3L, length(img_files))
    img_list <- as.list(sample(img_files, n_s))
  } else if (is.character(img) && length(img) > 1L) {
    if (!is.null(sample) && sample > 0L && sample < length(img)) {
      img_list <- as.list(sample(img, as.integer(sample)))
    } else {
      img_list <- as.list(img)
    }
  } else if (is.list(img) && !inherits(img, c("image", "Image"))) {
    img_list <- img
  } else {
    img_list <- list(img)
  }

  # Multi-Image Training Mode
  if (length(img_list) > 1L) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Multi-Image Few-Shot Training (pliman)")
      cli::cli_alert_info("Training visual prototypes across {length(img_list)} reference images for robust generalization.")
    }

    all_prototypes_sem <- list()
    all_prototypes_hr0 <- list()
    all_areas <- numeric()
    all_min_areas <- numeric()
    all_max_areas <- numeric()
    all_ref_boxes <- list()
    all_ref_polys <- list()
    all_ref_paths <- character()
    all_ref_imgs <- list()
    all_ref_points <- list()
    all_exemplar_colors <- list()
    all_radii <- numeric()
    stop_early <- FALSE

    for (s in seq_along(img_list)) {
      curr_target <- img_list[[s]]
      ref_id <- if (is.character(curr_target)) as.character(curr_target[1]) else paste0("img_", s)
      ref_norm <- if (is.character(curr_target)) {
        if (file.exists(curr_target[1])) {
          normalizePath(curr_target[1], winslash = "/", mustWork = FALSE)
        } else {
          m_fls <- list.files(".", pattern = paste0("^", basename(curr_target[1]), "\\."), full.names = TRUE)
          if (length(m_fls) > 0L) normalizePath(m_fls[1], winslash = "/", mustWork = FALSE) else ref_id
        }
      } else {
        ref_id
      }

      curr_im <- if (inherits(curr_target, c("image", "Image"))) {
        curr_target
      } else if (is.character(curr_target)) {
        image_import(curr_target[1])
      } else {
        as_image(curr_target)
      }

      canon_id <- if (ref_norm != ref_id && file.exists(ref_norm)) ref_norm else ref_id
      all_ref_paths <- c(all_ref_paths, canon_id)
      all_ref_imgs[[canon_id]] <- curr_im

      c_dims <- dim(curr_im)
      cw <- c_dims[1]; ch <- c_dims[2]

      curr_points <- NULL
      if (is.list(points) && length(points) >= s) {
        curr_points <- points[[s]]
      } else if (is.matrix(points) && s == 1L) {
        curr_points <- points
      }

      img_proto_sem <- list()
      img_proto_hr0 <- list()
      img_areas <- numeric()
      img_min_areas <- numeric()
      img_max_areas <- numeric()

      repeat {
        if (is.null(curr_points)) {
          repeat {
            if (isTRUE(verbose)) {
              if (s == 1L || length(all_prototypes_sem) == 0L) {
                cli::cli_alert_info("Image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}): Click on 2 to 5 exemplar instances of {.val {label}} in the plot window.")
                cli::cli_alert_info("Press <Esc> or right-click when finished.")
              } else {
                cli::cli_alert_info("Image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}): Click on exemplar instances of {.val {label}}, or press <Esc> to finish or continue.")
              }
            }
            plot(curr_im)
            pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.3), error = function(e) NULL)

            if (is.null(pts) || length(pts$x) == 0L) {
              if (s == 1L && length(all_prototypes_sem) == 0L) {
                if (!interactive()) {
                  cli::cli_abort("No exemplar points selected on first reference image.")
                }
                cli::cli_alert_warning("No point clicked on first image ({basename(as.character(curr_target[1]))}) (<Esc> pressed).")
                cli::cli_bullets(c(
                  "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image",
                  "*" = "{.bold [s] / stop}: Abort model fitting"
                ))
                ans <- tolower(trimws(readline("Select action [r = retry / s = stop, default: r]: ")))
                if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar")) {
                  cli::cli_abort("Model fitting cancelled by user.")
                } else {
                  next
                }
              } else {
                if (interactive()) {
                  cli::cli_alert_warning("No point clicked on image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}) (<Esc> pressed).")
                  cli::cli_bullets(c(
                    "*" = "{.bold [s] / stop}: Finalize model now with the {length(all_ref_paths)} reference image(s) completed so far",
                    "*" = "{.bold [c] / continue}: Continue to remaining images (auto-detect on this image)",
                    "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image"
                  ))
                  ans <- tolower(trimws(readline("Select action [s = stop / c = continue / r = retry, default: s]: ")))
                  if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar", "")) {
                    cli::cli_alert_info("Finalizing model with the {length(all_ref_paths)} reference image(s) completed so far...")
                    stop_early <- TRUE
                    break
                  } else if (ans %in% c("r", "retry", "tentar", "t")) {
                    next
                  } else {
                    curr_points <- NULL
                    break
                  }
                } else {
                  curr_points <- NULL
                  break
                }
              }
            } else {
              curr_points <- cbind(x = pts$x, y = pts$y)
              break
            }
          }

          if (isTRUE(stop_early)) {
            break
          }
        } else if (is.vector(curr_points) && length(curr_points) >= 2L && length(curr_points) %% 2 == 0L) {
          half <- length(curr_points) / 2
          curr_points <- cbind(x = curr_points[seq_len(half)], y = curr_points[(half + 1L):length(curr_points)])
        } else if (is.data.frame(curr_points)) {
          curr_points <- as.matrix(curr_points[, 1:2])
        }

        c_num_arr <- image_data(curr_im, type = "numeric")
        pre_proto <- if (is.null(curr_points) && length(all_prototypes_sem) > 0L) {
          list(
            semantic = all_prototypes_sem,
            hr0 = all_prototypes_hr0,
            avg_area = if (length(all_areas) > 0) mean(all_areas) else 100.0,
            min_area = if (length(all_min_areas) > 0) min(all_min_areas) else 10.0,
            max_area = if (length(all_max_areas) > 0) max(all_max_areas) else 500.0
          )
        } else {
          NULL
        }

        r_boxes <- data.frame()
        r_polys <- list()

        if (identical(method, "adaptive")) {
          norm_arr <- image_data(curr_im, type = "normalized")
          if (isTRUE(verbose)) {
            im_name_lbl <- basename(as.character(curr_target[1]))
            if (!is.null(curr_points) && nrow(curr_points) > 0L) {
              cli::cli_alert_success("Collected {nrow(curr_points)} exemplar point{?s} on {.file {im_name_lbl}}.")
            } else {
              cli::cli_alert_info("Auto-detecting on {.file {im_name_lbl}} using learned adaptive color model.")
            }
          }

          # Extract exemplar colors and calibrate scale from curr_points
          if (!is.null(curr_points) && nrow(curr_points) > 0L) {
            angles <- seq(0, 2 * pi, length.out = 16)[-16]
            max_search <- min(140L, as.integer(round(min(cw, ch) * 0.20)))
            for (k in seq_len(nrow(curr_points))) {
              px <- as.integer(round(curr_points[k, 1]))
              py <- as.integer(round(curr_points[k, 2]))
              if (px < 1L || px > cw || py < 1L || py > ch) next
              c_rgb <- norm_arr[px, py, ]
              r_search <- 5:max_search
              diffs <- numeric(length(r_search))
              for (ri in seq_along(r_search)) {
                r_rad <- r_search[ri]
                xs <- pmax(1L, pmin(cw, as.integer(round(px + r_rad * cos(angles)))))
                ys <- pmax(1L, pmin(ch, as.integer(round(py + r_rad * sin(angles)))))
                d_rgb <- sqrt((norm_arr[cbind(xs, ys, 1)] - c_rgb[1])^2 +
                              (norm_arr[cbind(xs, ys, 2)] - c_rgb[2])^2 +
                              (norm_arr[cbind(xs, ys, 3)] - c_rgb[3])^2)
                diffs[ri] <- mean(d_rgb)
              }
              grad <- c(0, diff(diffs))
              start_idx <- min(8L, length(grad))
              edge_r <- r_search[which.max(grad[start_idx:length(grad)]) + (start_idx - 1L)]
              all_radii <- c(all_radii, edge_r)

              rk <- max(3L, as.integer(round(edge_r * 0.35)))
              sub_r <- norm_arr[max(1L, px - rk):min(cw, px + rk), max(1L, py - rk):min(ch, py + rk), 1]
              sub_g <- norm_arr[max(1L, px - rk):min(cw, px + rk), max(1L, py - rk):min(ch, py + rk), 2]
              sub_b <- norm_arr[max(1L, px - rk):min(cw, px + rk), max(1L, py - rk):min(ch, py + rk), 3]
              all_exemplar_colors[[length(all_exemplar_colors) + 1L]] <- c(mean(sub_r), mean(sub_g), mean(sub_b))
            }
          }

          cur_col_mat <- if (length(all_exemplar_colors) > 0L) do.call(rbind, all_exemplar_colors) else matrix(c(0.5, 0.5, 0.5), nrow = 1)
          cur_target_rgb <- colMeans(cur_col_mat)
          cur_rgb_sd <- if (nrow(cur_col_mat) > 1L) pmax(apply(cur_col_mat, 2, stats::sd), 0.05) else c(0.15, 0.15, 0.15)
          cur_median_rad <- if (length(all_radii) > 0L) median(all_radii) else 20.0

          # Fast Adaptive Gaussian Mahalanobis similarity map
          dr <- (norm_arr[, , 1] - cur_target_rgb[1]) / cur_rgb_sd[1]
          dg <- (norm_arr[, , 2] - cur_target_rgb[2]) / cur_rgb_sd[2]
          db <- (norm_arr[, , 3] - cur_target_rgb[3]) / cur_rgb_sd[3]
          dist_sq <- (dr^2 + dg^2 + db^2) / 3.0
          sim_raw <- exp(-0.5 * dist_sq)

          blur_sigma <- max(2, as.integer(round(cur_median_rad * 0.12)))
          sim_im <- as_image(matrix(sim_raw, cw, ch), colormode = "Grayscale")
          sim_blur <- image_blur(sim_im, sigma = blur_sigma)
          sim_mat <- image_data(sim_blur, type = "normalized")

          adapt_t <- if (conf_threshold > 0.4) 0.25 else conf_threshold
          bin_mask <- sim_mat >= adapt_t
          lbl <- bwlabel_cpp(bin_mask)
          counts <- table(lbl[lbl > 0])
          min_a <- if (!is.null(min_area) && min_area > 0) min_area else max(10, round(pi * (cur_median_rad * 0.4)^2))
          max_a <- if (!is.null(max_area) && max_area > 0) max_area else round(pi * (cur_median_rad * 2.5)^2)
          valid <- as.integer(names(counts[counts >= min_a & counts <= max_a]))

          cand_b <- list()
          cand_c <- numeric()
          conts_all <- if (length(valid) > 0L) extract_contours_cpp(lbl) else list()
          cand_p <- list()

          for (obj in valid) {
            coords <- which(lbl == obj, arr.ind = TRUE)
            bx1 <- min(coords[, 1]); bx2 <- max(coords[, 1])
            by1 <- min(coords[, 2]); by2 <- max(coords[, 2])
            bw <- bx2 - bx1 + 1L; bh <- by2 - by1 + 1L
            ar <- bw / pmax(1L, bh)
            if (ar < 0.38 || ar > 2.60) next

            if (bw >= 12L && bh >= 12L) {
              sub_p <- norm_arr[bx1:bx2, by1:by2, 1]
              nr <- nrow(sub_p); nc <- ncol(sub_p)
              mr <- max(2L, as.integer(round(nr * 0.15))); mc <- max(2L, as.integer(round(nc * 0.15)))
              c_sub <- sub_p[mr:(nr - mr), mc:(nc - mc)]
              b_sub <- c(sub_p[1:mr, ], sub_p[(nr - mr + 1L):nr, ],
                         sub_p[, 1:mc], sub_p[, (nc - mc + 1L):nc])
              if (mean(c_sub) - mean(b_sub) < 0.035) next
            }

            obj_conf <- mean(sim_mat[coords])
            if (obj_conf < adapt_t) next

            poly_norm <- NULL
            if (length(conts_all) >= obj && is.matrix(conts_all[[obj]]) && nrow(conts_all[[obj]]) >= 3L) {
              poly_pts <- conts_all[[obj]]
              gx <- poly_pts[, 1] / cw; gy <- poly_pts[, 2] / ch
              if (length(gx) > 32L) {
                idx <- round(seq(1, length(gx), length.out = 32))
                gx <- gx[idx]; gy <- gy[idx]
              }
              poly_norm <- as.vector(rbind(gx, gy))
            }

            cand_b[[length(cand_b) + 1L]] <- c(bx1, by1, bx2, by2)
            cand_c <- c(cand_c, obj_conf)
            cand_p[[length(cand_p) + 1L]] <- poly_norm
          }

          if (length(cand_b) > 0L) {
            box_mat <- do.call(rbind, cand_b)
            if (nrow(box_mat) > 1L) {
              keep_idx <- .nms_boxes(box_mat, cand_c, iou_thresh = 0.35)
              box_mat <- box_mat[keep_idx, , drop = FALSE]
              cand_c <- cand_c[keep_idx]
              cand_p <- cand_p[keep_idx]
            }
            r_boxes <- data.frame(
              xmin = as.numeric(box_mat[, 1]),
              ymin = as.numeric(box_mat[, 2]),
              xmax = as.numeric(box_mat[, 3]),
              ymax = as.numeric(box_mat[, 4]),
              conf = round(cand_c, 3),
              class_id = 0L,
              class_name = label,
              stringsAsFactors = FALSE
            )
            r_polys <- cand_p
          }
        } else {
          # SAM 2.1 / PerSAM branch
          if (isTRUE(verbose)) {
            im_name_lbl <- basename(as.character(curr_target[1]))
            if (!is.null(curr_points) && nrow(curr_points) > 0L) {
              cli::cli_alert_success("Collected {nrow(curr_points)} exemplar point{?s} on {.file {im_name_lbl}}.")
            } else {
              cli::cli_alert_info("Auto-detecting on {.file {im_name_lbl}} using {length(all_prototypes_sem)} learned prototype{?s}.")
            }
            cli::cli_progress_step(
              msg = "Running SAM 2.1 inference on {.file {im_name_lbl}} [{toupper(engine)}]... ",
              msg_done = "SAM 2.1 inference complete on {.file {im_name_lbl}}"
            )
          }

          p_res <- tryCatch({
            .run_persam(
              mat = c_num_arr,
              exemplar_points = curr_points,
              precomputed_prototypes = pre_proto,
              sim_threshold = conf_threshold,
              min_dist = min_dist,
              max_objects = max_objects,
              feat_res = 256L,
              engine = engine,
              device_id = device_id,
              fill_hull = TRUE,
              verbose = FALSE
            )
          }, error = function(e) {
            if (isTRUE(verbose)) {
              cli::cli_alert_warning("Exemplar matching error on image {s}: {e$message}")
            }
            NULL
          })

          if (isTRUE(verbose)) {
            cli::cli_progress_done()
          }

          if (!is.null(p_res) && !is.null(p_res$prototypes)) {
            img_proto_sem <- p_res$prototypes$semantic
            img_proto_hr0 <- p_res$prototypes$hr0
            if (!is.null(p_res$prototypes$avg_area)) {
              img_areas <- as.numeric(p_res$prototypes$avg_area)
            }
            if (!is.null(p_res$prototypes$min_area)) {
              img_min_areas <- as.numeric(p_res$prototypes$min_area)
            }
            if (!is.null(p_res$prototypes$max_area)) {
              img_max_areas <- as.numeric(p_res$prototypes$max_area)
            }
          }

          if (!is.null(p_res) && nrow(p_res$boxes) > 0L) {
            r_boxes <- data.frame(
              xmin = as.numeric(p_res$boxes$xmin),
              ymin = as.numeric(p_res$boxes$ymin),
              xmax = as.numeric(p_res$boxes$xmax),
              ymax = as.numeric(p_res$boxes$ymax),
              conf = round(as.numeric(p_res$boxes$score), 3),
              class_id = 0L,
              class_name = label,
              stringsAsFactors = FALSE
            )
            for (ci in seq_along(p_res$contours)) {
              cnt <- p_res$contours[[ci]]
              if (is.matrix(cnt) && nrow(cnt) >= 3L) {
                gx <- cnt[, 1] / cw; gy <- cnt[, 2] / ch
                if (length(gx) > 32L) {
                  idx <- round(seq(1, length(gx), length.out = 32))
                  gx <- gx[idx]; gy <- gy[idx]
                }
                r_polys[[ci]] <- as.vector(rbind(gx, gy))
              }
            }
          }
        }

        if (isTRUE(interactive)) {
          rev <- .interactive_fewshot_review_step(
            curr_im = curr_im,
            r_boxes = r_boxes,
            r_polys = r_polys,
            curr_points = curr_points,
            c_num_arr = c_num_arr,
            cw = cw,
            ch = ch,
            pre_proto = pre_proto,
            latest_proto = if (!is.null(p_res) && !is.null(p_res$prototypes)) p_res$prototypes else NULL,
            task = task,
            label = label,
            conf_threshold = conf_threshold,
            radius = radius,
            min_dist = min_dist,
            max_objects = max_objects,
            method = method,
            engine = engine,
            device_id = device_id,
            rainbow = rainbow,
            show_class = show_class,
            show_conf = show_conf,
            badge = badge,
            img_name = basename(as.character(curr_target[1])),
            img_idx = s,
            total_imgs = length(img_list)
          )

          if (identical(rev$action, "resample")) {
            curr_points <- NULL
            next
          }

          r_boxes <- rev$boxes
          r_polys <- rev$polygons
          curr_points <- rev$points
          if (!is.null(rev$conf)) conf_threshold <- rev$conf

          if (identical(rev$action, "stop")) {
            stop_early <- TRUE
          }
        }
        break
      }

      if (isTRUE(stop_early) && (is.null(r_boxes) || nrow(r_boxes) == 0L) && length(all_ref_boxes) > 0L) {
        break
      }

      # Commit learned prototypes from this image
      if (length(img_proto_sem) > 0L) {
        all_prototypes_sem <- c(all_prototypes_sem, img_proto_sem)
        all_prototypes_hr0 <- c(all_prototypes_hr0, img_proto_hr0)
        all_areas <- c(all_areas, img_areas)
        all_min_areas <- c(all_min_areas, img_min_areas)
        all_max_areas <- c(all_max_areas, img_max_areas)
      }

      all_ref_boxes[[canon_id]] <- r_boxes
      if (ref_id != canon_id) all_ref_boxes[[ref_id]] <- r_boxes
      all_ref_polys[[canon_id]] <- r_polys
      if (ref_id != canon_id) all_ref_polys[[ref_id]] <- r_polys
      if (!is.null(curr_points)) {
        all_ref_points[[canon_id]] <- curr_points
        if (ref_id != canon_id) all_ref_points[[ref_id]] <- curr_points
      }

      if (isTRUE(stop_early)) {
        break
      }
    }

    if (length(all_ref_paths) == 0L) {
      cli::cli_abort("Model fitting aborted. No reference images were annotated.")
    }

    all_areas <- as.numeric(unlist(all_areas, use.names = FALSE))
    all_min_areas <- as.numeric(unlist(all_min_areas, use.names = FALSE))
    all_max_areas <- as.numeric(unlist(all_max_areas, use.names = FALSE))
    all_areas <- all_areas[!is.na(all_areas)]
    all_min_areas <- all_min_areas[!is.na(all_min_areas)]
    all_max_areas <- all_max_areas[!is.na(all_max_areas)]

    merged_prototypes <- list(
      semantic = all_prototypes_sem,
      hr0 = all_prototypes_hr0,
      avg_area = if (length(all_areas) > 0) mean(all_areas) else 100.0,
      min_area = if (length(all_min_areas) > 0) min(all_min_areas) else 10.0,
      max_area = if (length(all_max_areas) > 0) max(all_max_areas) else 500.0
    )

    primary_id <- all_ref_paths[1]
    primary_im <- all_ref_imgs[[primary_id]]
    primary_boxes <- all_ref_boxes[[primary_id]]
    primary_polys <- all_ref_polys[[primary_id]]

    final_target_rgb <- if (length(all_exemplar_colors) > 0L) colMeans(do.call(rbind, all_exemplar_colors)) else c(0.5, 0.5, 0.5)
    final_rgb_sd <- if (length(all_exemplar_colors) > 1L) pmax(apply(do.call(rbind, all_exemplar_colors), 2, stats::sd), 0.05) else c(0.15, 0.15, 0.15)
    final_median_rad <- if (length(all_radii) > 0L) median(all_radii) else 20.0

    model_obj <- structure(
      list(
        label = label,
        task = task,
        target_rgb = final_target_rgb,
        rgb_sd = final_rgb_sd,
        median_radius = final_median_rad,
        median_area = pi * final_median_rad^2,
        conf_threshold = conf_threshold,
        min_area = if (!is.null(min_area)) min_area else 0,
        max_area = max_area,
        exemplar_points = if (length(all_ref_points) > 0L) all_ref_points[[1]] else NULL,
        exemplar_points_map = all_ref_points,
        detected_boxes = primary_boxes,
        detected_polygons = primary_polys,
        detected_labels = NULL,
        detected_mask = NULL,
        similarity_map = NULL,
        prototypes = if (identical(method, "sam")) merged_prototypes else NULL,
        reference_image_path = primary_id,
        reference_image_paths = all_ref_paths,
        reference_boxes_map = all_ref_boxes,
        reference_polys_map = all_ref_polys,
        reference_image = primary_im,
        reference_images = all_ref_imgs,
        method = method,
        engine = engine,
        rainbow = rainbow,
        show_class = show_class,
        show_conf = show_conf,
        badge = badge
      ),
      class = "yolo_fewshot_model"
    )

    if (isTRUE(verbose)) {
      unique_ref_keys <- intersect(all_ref_paths, names(all_ref_boxes))
      if (length(unique_ref_keys) == 0L) unique_ref_keys <- names(all_ref_imgs)
      total_det <- sum(vapply(unique_ref_keys, function(k) nrow(all_ref_boxes[[k]]), integer(1)))
      proto_info <- if (identical(method, "sam")) paste0("Learned ", length(all_prototypes_sem), " prototypes (") else "Adaptive color model ("
      cli::cli_alert_success("Model fitted across {length(unique_ref_keys)} image{?s}! {proto_info}{total_det} reference instances segmented).")
    }

    if (isTRUE(plot) && !is.null(primary_im)) {
      plot(model_obj, rainbow = rainbow, show_class = show_class, show_conf = show_conf, badge = badge)
    }

    return(invisible(model_obj))
  }

  # Single Image Training Mode
  img_single <- img_list[[1]]
  im <- if (inherits(img_single, c("image", "Image"))) {
    img_single
  } else if (is.character(img_single)) {
    image_import(img_single[1])
  } else {
    as_image(img_single)
  }

  dims <- dim(im)
  w <- dims[1]
  h <- dims[2]

  # Interactive clicking if points not provided
  if (is.null(points)) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Few-Shot Exemplar Selection (pliman)")
      cli::cli_alert_info("Click on 3 to 5 sample instances of {.val {label}} in the plot window.")
      cli::cli_alert_info("Press <Esc> or right-click when finished.")
    }
    repeat {
      plot(im)
      pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.3), error = function(e) NULL)
      if (is.null(pts) || length(pts$x) == 0L) {
        if (!interactive()) {
          cli::cli_abort("No exemplar points selected.")
        }
        cli::cli_alert_warning("No point clicked on image (<Esc> pressed).")
        cli::cli_bullets(c(
          "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image",
          "*" = "{.bold [s] / stop}: Abort model fitting"
        ))
        ans <- tolower(trimws(readline("Select action [r = retry / s = stop, default: r]: ")))
        if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar")) {
          cli::cli_abort("Model fitting cancelled by user.")
        } else {
          next
        }
      } else {
        points <- cbind(x = pts$x, y = pts$y)
        break
      }
    }
  } else if (is.vector(points) && length(points) >= 2L && length(points) %% 2 == 0L) {
    half <- length(points) / 2
    points <- cbind(x = points[seq_len(half)], y = points[(half + 1L):length(points)])
  } else if (is.data.frame(points)) {
    points <- as.matrix(points[, 1:2])
  }

  n_exemplars <- nrow(points)
  if (isTRUE(verbose)) {
    cli::cli_alert_success("Collected {n_exemplars} exemplar point{?s} for class {.val {label}}.")
  }

  norm_arr <- image_data(im, type = "normalized")
  num_arr <- image_data(im, type = "numeric")

  # 1. Scale auto-calibration from exemplars
  angles <- seq(0, 2 * pi, length.out = 16)[-16]
  max_search <- min(140L, as.integer(round(min(w, h) * 0.20)))
  radii <- numeric(n_exemplars)

  for (k in seq_len(n_exemplars)) {
    px <- as.integer(round(points[k, 1]))
    py <- as.integer(round(points[k, 2]))
    c_rgb <- norm_arr[px, py, ]
    r_search <- 5:max_search
    diffs <- numeric(length(r_search))
    for (i in seq_along(r_search)) {
      r <- r_search[i]
      xs <- pmax(1L, pmin(w, as.integer(round(px + r * cos(angles)))))
      ys <- pmax(1L, pmin(h, as.integer(round(py + r * sin(angles)))))
      d_rgb <- sqrt((norm_arr[cbind(xs, ys, 1)] - c_rgb[1])^2 +
                    (norm_arr[cbind(xs, ys, 2)] - c_rgb[2])^2 +
                    (norm_arr[cbind(xs, ys, 3)] - c_rgb[3])^2)
      diffs[i] <- mean(d_rgb)
    }
    grad <- c(0, diff(diffs))
    start_idx <- min(8L, length(grad))
    edge_r <- r_search[which.max(grad[start_idx:length(grad)]) + (start_idx - 1L)]
    radii[k] <- edge_r
  }

  median_rad <- if (!is.null(radius)) as.numeric(radius) else max(8, median(radii))
  median_area <- pi * median_rad^2

  if (is.null(min_area) || identical(min_area, 25)) {
    min_area <- 0
  }

  # Extract exemplar color signatures
  exemplar_colors <- list()
  for (k in seq_len(n_exemplars)) {
    px <- as.integer(round(points[k, 1]))
    py <- as.integer(round(points[k, 2]))
    rk <- max(3L, as.integer(round(radii[k] * 0.35)))
    sub_r <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 1]
    sub_g <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 2]
    sub_b <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 3]
    exemplar_colors[[k]] <- c(mean(sub_r), mean(sub_g), mean(sub_b))
  }
  col_mat <- do.call(rbind, exemplar_colors)
  target_rgb <- colMeans(col_mat)
  rgb_sd <- pmax(apply(col_mat, 2, stats::sd), 0.05)

  # 2. Execute segmentation engine
  used_method <- method
  df_boxes <- data.frame()
  cand_polys <- list()
  det_labels <- NULL
  det_mask <- NULL
  det_sim <- NULL

  if (method %in% c("sam", "persam")) {
    used_method <- "sam"
    persam_res <- tryCatch({
      .run_persam(
        mat = num_arr,
        exemplar_points = points,
        sim_threshold = conf_threshold,
        min_dist = if (!is.null(radius)) min(min_dist, radius * 0.4) else min_dist,
        max_objects = max_objects,
        feat_res = 256L,
        engine = engine,
        device_id = device_id,
        fill_hull = TRUE,
        verbose = verbose
      )
    }, error = function(e) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning("PerSAM execution encountered: {e$message}. Falling back to adaptive morphology.")
      }
      NULL
    })

    if (!is.null(persam_res) && nrow(persam_res$boxes) > 0L) {
      raw_boxes <- persam_res$boxes
      raw_conts <- persam_res$contours
      det_labels <- persam_res$labels
      det_mask <- persam_res$mask
      det_sim <- persam_res$similarity_map

      # Filter by scale if specified
      valid_idx <- seq_len(nrow(raw_boxes))
      if (!is.null(min_area) && min_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- which(b_areas >= min_area)
      }
      if (!is.null(max_area) && max_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- intersect(valid_idx, which(b_areas <= max_area))
      }

      if (length(valid_idx) > 0L) {
        raw_boxes <- raw_boxes[valid_idx, , drop = FALSE]
        raw_conts <- raw_conts[valid_idx]

        df_boxes <- data.frame(
          xmin = as.numeric(raw_boxes$xmin),
          ymin = as.numeric(raw_boxes$ymin),
          xmax = as.numeric(raw_boxes$xmax),
          ymax = as.numeric(raw_boxes$ymax),
          conf = round(as.numeric(raw_boxes$score), 3),
          class_id = 0L,
          class_name = label,
          stringsAsFactors = FALSE
        )

        for (i in seq_along(raw_conts)) {
          cnt <- raw_conts[[i]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            gx <- cnt[, 1] / w
            gy <- cnt[, 2] / h
            if (length(gx) > 32L) {
              idx <- round(seq(1, length(gx), length.out = 32))
              gx <- gx[idx]; gy <- gy[idx]
            }
            cand_polys[[i]] <- as.vector(rbind(gx, gy))
          } else {
            cand_polys[[i]] <- NULL
          }
        }
      }
    } else {
      used_method <- "adaptive"
    }
  }

  if (identical(used_method, "adaptive") || nrow(df_boxes) == 0L) {
    used_method <- "adaptive"
    # Adaptive Gaussian Mahalanobis similarity map
    dr <- (norm_arr[, , 1] - target_rgb[1]) / rgb_sd[1]
    dg <- (norm_arr[, , 2] - target_rgb[2]) / rgb_sd[2]
    db <- (norm_arr[, , 3] - target_rgb[3]) / rgb_sd[3]
    dist_sq <- (dr^2 + dg^2 + db^2) / 3.0
    sim_raw <- exp(-0.5 * dist_sq)

    blur_sigma <- max(2, as.integer(round(median_rad * 0.12)))
    sim_im <- as_image(matrix(sim_raw, w, h), colormode = "Grayscale")
    sim_blur <- image_blur(sim_im, sigma = blur_sigma)
    sim_mat <- image_data(sim_blur, type = "normalized")

    adapt_t <- if (conf_threshold > 0.4) 0.25 else conf_threshold
    bin_mask <- sim_mat >= adapt_t
    lbl <- bwlabel_cpp(bin_mask)
    counts <- table(lbl[lbl > 0])
    valid <- as.integer(names(counts[counts >= min_area & counts <= max_area]))

    cand_b <- list()
    cand_c <- numeric()
    conts_all <- if (length(valid) > 0L) extract_contours_cpp(lbl) else list()
    cand_p <- list()

    for (obj in valid) {
      coords <- which(lbl == obj, arr.ind = TRUE)
      bx1 <- min(coords[, 1]); bx2 <- max(coords[, 1])
      by1 <- min(coords[, 2]); by2 <- max(coords[, 2])
      bw <- bx2 - bx1 + 1L; bh <- by2 - by1 + 1L
      ar <- bw / pmax(1L, bh)
      if (ar < 0.38 || ar > 2.60) next

      if (bw >= 12L && bh >= 12L) {
        sub_p <- norm_arr[bx1:bx2, by1:by2, 1]
        nr <- nrow(sub_p); nc <- ncol(sub_p)
        mr <- max(2L, as.integer(round(nr * 0.15))); mc <- max(2L, as.integer(round(nc * 0.15)))
        c_sub <- sub_p[mr:(nr - mr), mc:(nc - mc)]
        b_sub <- c(sub_p[1:mr, ], sub_p[(nr - mr + 1L):nr, ],
                   sub_p[, 1:mc], sub_p[, (nc - mc + 1L):nc])
        if (mean(c_sub) - mean(b_sub) < 0.035) next
      }

      obj_conf <- mean(sim_mat[coords])
      if (obj_conf < adapt_t) next

      poly_norm <- NULL
      if (length(conts_all) >= obj && is.matrix(conts_all[[obj]]) && nrow(conts_all[[obj]]) >= 3L) {
        poly_pts <- conts_all[[obj]]
        gx <- poly_pts[, 1] / w; gy <- poly_pts[, 2] / h
        if (length(gx) > 32L) {
          idx <- round(seq(1, length(gx), length.out = 32))
          gx <- gx[idx]; gy <- gy[idx]
        }
        poly_norm <- as.vector(rbind(gx, gy))
      }

      cand_b[[length(cand_b) + 1L]] <- c(bx1, by1, bx2, by2)
      cand_c <- c(cand_c, obj_conf)
      cand_p[[length(cand_p) + 1L]] <- poly_norm
    }

    if (length(cand_b) > 0L) {
      box_mat <- do.call(rbind, cand_b)
      if (nrow(box_mat) > 1L) {
        keep_idx <- .nms_boxes(box_mat, cand_c, iou_thresh = 0.35)
        box_mat <- box_mat[keep_idx, , drop = FALSE]
        cand_c <- cand_c[keep_idx]
        cand_p <- cand_p[keep_idx]
      }
      df_boxes <- data.frame(
        xmin = as.numeric(box_mat[, 1]),
        ymin = as.numeric(box_mat[, 2]),
        xmax = as.numeric(box_mat[, 3]),
        ymax = as.numeric(box_mat[, 4]),
        conf = round(cand_c, 3),
        class_id = 0L,
        class_name = label,
        stringsAsFactors = FALSE
      )
      cand_polys <- cand_p
    }
  }

  if (isTRUE(interactive)) {
    rev <- .interactive_fewshot_review_step(
      curr_im = im,
      r_boxes = df_boxes,
      r_polys = cand_polys,
      curr_points = points,
      c_num_arr = num_arr,
      cw = w,
      ch = h,
      pre_proto = NULL,
      latest_proto = if (exists("persam_res") && !is.null(persam_res$prototypes)) persam_res$prototypes else NULL,
      task = task,
      label = label,
      conf_threshold = conf_threshold,
      radius = radius,
      min_dist = min_dist,
      max_objects = max_objects,
      method = used_method,
      engine = engine,
      device_id = device_id,
      rainbow = rainbow,
      show_class = show_class,
      show_conf = show_conf,
      badge = badge,
      img_name = if (is.character(img_single)) basename(img_single[1]) else "reference",
      img_idx = 1L,
      total_imgs = 1L
    )
    df_boxes <- rev$boxes
    cand_polys <- rev$polygons
    points <- rev$points
    conf_threshold <- rev$conf
  }

  ref_path <- if (is.character(img) && file.exists(img[1])) normalizePath(img[1], winslash = "/", mustWork = FALSE) else NULL
  prototypes <- if (exists("persam_res") && !is.null(persam_res)) persam_res$prototypes else NULL

  ref_boxes_map <- list()
  ref_polys_map <- list()
  if (!is.null(ref_path)) {
    ref_boxes_map[[ref_path]] <- df_boxes
    ref_polys_map[[ref_path]] <- cand_polys
  }

  model_obj <- structure(
    list(
      label = label,
      task = task,
      target_rgb = target_rgb,
      rgb_sd = rgb_sd,
      median_radius = median_rad,
      median_area = median_area,
      conf_threshold = conf_threshold,
      min_area = min_area,
      max_area = max_area,
      exemplar_points = points,
      exemplar_points_map = if (!is.null(ref_path)) stats::setNames(list(points), ref_path) else list(points),
      detected_boxes = df_boxes,
      detected_polygons = cand_polys,
      detected_labels = det_labels,
      detected_mask = det_mask,
      similarity_map = det_sim,
      prototypes = prototypes,
      reference_image_path = ref_path,
      reference_image_paths = if (!is.null(ref_path)) ref_path else character(0),
      reference_boxes_map = ref_boxes_map,
      reference_polys_map = ref_polys_map,
      reference_image = im,
      reference_images = if (!is.null(ref_path)) stats::setNames(list(im), ref_path) else list(im),
      method = used_method,
      engine = engine,
      rainbow = rainbow,
      show_class = show_class,
      show_conf = show_conf,
      badge = badge
    ),
    class = "yolo_fewshot_model"
  )

  if (isTRUE(verbose)) {
    n_det <- if (task == "segment" && length(cand_polys) > 0L) sum(!vapply(cand_polys, is.null, logical(1))) else nrow(df_boxes)
    cli::cli_alert_success("Model fitted! {n_det} instances {if (task == 'segment') 'segmented' else 'detected'} on reference image.")
  }

  if (isTRUE(plot) && !isTRUE(interactive)) {
    plot(model_obj, rainbow = rainbow, show_class = show_class, show_conf = show_conf, badge = badge)
  }

  invisible(model_obj)
}

# Cross-class Non-Maximum Suppression (competitive disambiguation across classes)
.cross_class_nms <- function(boxes, polygons = NULL, iou_thresh = 0.40) {
  if (is.null(boxes) || nrow(boxes) <= 1L) {
    return(list(boxes = boxes, polygons = polygons))
  }
  # Sort all candidate boxes across all classes by confidence descending
  ord <- order(boxes$conf, decreasing = TRUE)
  boxes <- boxes[ord, , drop = FALSE]
  if (!is.null(polygons) && length(polygons) == nrow(boxes)) {
    polygons <- polygons[ord]
  }

  x1 <- boxes$xmin
  y1 <- boxes$ymin
  x2 <- boxes$xmax
  y2 <- boxes$ymax
  areas <- pmax(0, x2 - x1) * pmax(0, y2 - y1)

  n <- nrow(boxes)
  suppressed <- rep(FALSE, n)

  for (i in seq_len(n)) {
    if (suppressed[i] || i >= n) next
    for (j in (i + 1L):n) {
      if (suppressed[j]) next
      xx1 <- max(x1[i], x1[j])
      yy1 <- max(y1[i], y1[j])
      xx2 <- min(x2[i], x2[j])
      yy2 <- min(y2[i], y2[j])
      w_int <- max(0, xx2 - xx1)
      h_int <- max(0, yy2 - yy1)
      inter <- w_int * h_int
      if (inter > 0) {
        union_area <- areas[i] + areas[j] - inter
        iou <- if (union_area > 0) inter / union_area else 0
        min_area <- min(areas[i], areas[j])
        iom <- if (min_area > 0) inter / min_area else 0
        if (iou >= iou_thresh || iom >= 0.60) {
          suppressed[j] <- TRUE
        }
      }
    }
  }

  keep <- which(!suppressed)
  res_boxes <- boxes[keep, , drop = FALSE]
  res_polys <- if (!is.null(polygons) && length(polygons) == n) polygons[keep] else polygons
  list(boxes = res_boxes, polygons = res_polys)
}

#' @title Combine Few-Shot Models
#' @description Combines multiple `yolo_fewshot_model` objects into a multi-class few-shot model.
#' @param ... `yolo_fewshot_model` objects to combine.
#' @return A `yolo_fewshot_multimodel` object.
#' @export
c.yolo_fewshot_model <- function(...) {
  dots <- list(...)
  dots <- dots[!vapply(dots, is.null, logical(1))]
  all_submodels <- list()
  for (item in dots) {
    if (inherits(item, "yolo_fewshot_multimodel")) {
      all_submodels <- c(all_submodels, item$models)
    } else if (inherits(item, "yolo_fewshot_model")) {
      all_submodels <- c(all_submodels, list(item))
    }
  }
  if (length(all_submodels) == 0L) {
    return(NULL)
  }
  cls_names <- character(length(all_submodels))
  for (i in seq_along(all_submodels)) {
    c_name <- all_submodels[[i]]$label[1]
    cls_names[i] <- c_name
    all_submodels[[i]]$class_id <- as.integer(i - 1L)
    if (!is.null(all_submodels[[i]]$detected_boxes) && nrow(all_submodels[[i]]$detected_boxes) > 0L) {
      all_submodels[[i]]$detected_boxes$class_id <- as.integer(i - 1L)
    }
  }
  names(all_submodels) <- cls_names
  structure(
    list(
      label = cls_names,
      classes = cls_names,
      models = all_submodels,
      task = all_submodels[[1]]$task,
      conf_threshold = all_submodels[[1]]$conf_threshold,
      method = all_submodels[[1]]$method,
      engine = all_submodels[[1]]$engine,
      reference_image = all_submodels[[1]]$reference_image,
      reference_images = all_submodels[[1]]$reference_images
    ),
    class = c("yolo_fewshot_multimodel", "yolo_fewshot_model")
  )
}

#' @title Predict Few-Shot Instances on a New Image
#' @name yolo_fewshot_predict
#' @description
#' Applies a trained `yolo_fewshot_model` to detect and segment instances in a new image.
#'
#' @param model A `yolo_fewshot_model` fitted via [yolo_fewshot_fit()].
#' @param img An `image` object or file path to evaluate.
#' @param task Task override: `"segment"` (polygons) or `"detect"` (bounding boxes). Defaults to `model$task`.
#' @param conf_threshold Numeric scalar override for detection confidence. Defaults to model threshold.
#' @param max_objects Maximum number of instances to detect.
#' @param engine Execution engine: `"cpu"` or `"gpu"`. If `NULL`, uses `model$engine`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param plot Logical. Display detection results. Defaults to `FALSE`.
#' @param rainbow Logical. If `TRUE`, plots instances in distinct rainbow colors. Defaults to `FALSE`.
#' @param show_class Logical. Whether to display class name on labels/badges. Defaults to `TRUE`.
#' @param show_conf Logical. Whether to display confidence score on labels/badges. Defaults to `TRUE`.
#' @param badge Logical. Whether to draw a filled background badge behind labels. Defaults to `TRUE`.
#' @param ... Additional arguments passed to plotting methods or aliases (e.g. `show_labels`, `show_scores`, `label_box`).
#'
#' @return A `list` containing `boxes` (bounding boxes data.frame) and `polygons` (list of normalized coordinates).
#' @export
yolo_fewshot_predict <- function(model,
                                 img,
                                 task = model$task,
                                 conf_threshold = NULL,
                                 max_objects = NULL,
                                 engine = "gpu",
                                 device_id = -1,
                                 plot = FALSE,
                                 rainbow = TRUE,
                                 show_class = if (!is.null(model$show_class)) model$show_class else TRUE,
                                 show_conf = if (!is.null(model$show_conf)) model$show_conf else TRUE,
                                 badge = if (!is.null(model$badge)) model$badge else TRUE,
                                 ...) {
  dots <- list(...)
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)
  if (is.list(model) && !inherits(model, "yolo_fewshot_model")) {
    if (all(vapply(model, function(m) inherits(m, "yolo_fewshot_model"), logical(1)))) {
      model <- do.call(c, model)
    }
  }

  if (!inherits(model, "yolo_fewshot_model")) {
    cli::cli_abort("{.arg model} must be a {.cls yolo_fewshot_model} from {.fn yolo_fewshot_fit}.")
  }

  im <- if (inherits(img, c("image", "Image"))) {
    img
  } else if (is.character(img)) {
    image_import(img[1])
  } else {
    as_image(img)
  }

  dims <- dim(im)
  w <- dims[1]; h <- dims[2]

  task <- if (!missing(task) && !is.null(task)) match.arg(task, c("segment", "detect")) else if (!is.null(model$task)) model$task else "segment"
  conf_t <- if (!is.null(conf_threshold)) as.numeric(conf_threshold) else model$conf_threshold

  # Handle multi-class models
  if (inherits(model, "yolo_fewshot_multimodel") || !is.null(model$models)) {
    all_b <- list()
    all_p <- list()
    for (m_idx in seq_along(model$models)) {
      sub_m <- model$models[[m_idx]]
      sub_res <- yolo_fewshot_predict(
        model = sub_m,
        img = im,
        task = task,
        conf_threshold = conf_threshold,
        max_objects = max_objects,
        engine = engine,
        device_id = device_id,
        plot = FALSE,
        rainbow = FALSE
      )
      if (!is.null(sub_res$boxes) && nrow(sub_res$boxes) > 0L) {
        all_b[[length(all_b) + 1L]] <- sub_res$boxes
        all_p <- c(all_p, sub_res$polygons)
      }
    }
    comb_boxes <- if (length(all_b) > 0L) do.call(rbind, all_b) else data.frame()
    if (nrow(comb_boxes) > 0L) {
      nms_res <- .cross_class_nms(comb_boxes, all_p, iou_thresh = 0.40)
      comb_boxes <- nms_res$boxes
      all_p <- nms_res$polygons
    }
    res <- list(boxes = comb_boxes, polygons = all_p, image = im)
    if (isTRUE(plot)) {
      plot(im)
      if (nrow(res$boxes) > 0L) {
        classes_vec <- if (!is.null(model$classes)) model$classes else unique(res$boxes$class_name)
        pal <- if (isTRUE(rainbow)) grDevices::rainbow(max(1L, length(classes_vec)), s = 0.85, v = 0.95) else rep("#00FFCC", length(classes_vec))
        names(pal) <- classes_vec
        if (task == "segment" && length(res$polygons) > 0L) {
          for (i in seq_along(res$polygons)) {
            poly <- res$polygons[[i]]
            c_name <- if (i <= nrow(res$boxes)) res$boxes$class_name[i] else classes_vec[1]
            c_base <- if (c_name %in% names(pal)) pal[[c_name]] else "#00FFCC"
            if (!is.null(poly) && length(poly) >= 6L) {
              px <- poly[c(TRUE, FALSE)] * w
              py <- poly[c(FALSE, TRUE)] * h
              graphics::polygon(px, py, col = grDevices::adjustcolor(c_base, alpha.f = 0.40), border = c_base, lwd = 1.5)
            }
          }
        }
        for (r in seq_len(nrow(res$boxes))) {
          c_name <- res$boxes$class_name[r]
          c_col <- if (c_name %in% names(pal)) pal[[c_name]] else "#00FFCC"
          if (task == "detect" || length(res$polygons) == 0L) {
            graphics::rect(res$boxes$xmin[r], res$boxes$ymin[r], res$boxes$xmax[r], res$boxes$ymax[r], border = c_col, lwd = 2)
          }
          graphics::text(res$boxes$xmin[r], res$boxes$ymin[r] - 5,
                         labels = sprintf("%s %.2f", res$boxes$class_name[r], res$boxes$conf[r]),
                         col = c_col, cex = 0.8, pos = 4)
        }
        cls_counts <- table(res$boxes$class_name)
        sum_str <- paste(sprintf("%d %s", as.integer(cls_counts), names(cls_counts)), collapse = ", ")
        graphics::title(sprintf("Few-Shot Multi-Class: %s [%s]", sum_str, task),
                        col.main = if (isTRUE(rainbow)) "#00CC66" else "#00FFCC",
                        cex.main = 1.0)
      }
    }
    return(invisible(res))
  }

  # 1. Check if this is one of the exact reference images used during fit
  is_ref <- FALSE
  ref_match_id <- NULL
  if (is.character(img)) {
    norm_in <- normalizePath(img[1], winslash = "/", mustWork = FALSE)
    raw_in <- as.character(img[1])
    raw_sans <- tools::file_path_sans_ext(basename(raw_in))
    ref_paths <- model$reference_image_paths
    if (!is.null(ref_paths) && length(ref_paths) > 0L) {
      ref_sans <- tools::file_path_sans_ext(basename(ref_paths))
      if (raw_in %in% names(model$reference_boxes_map)) {
        is_ref <- TRUE
        ref_match_id <- raw_in
      } else if (norm_in %in% names(model$reference_boxes_map)) {
        is_ref <- TRUE
        ref_match_id <- norm_in
      } else if (raw_sans %in% ref_sans) {
        is_ref <- TRUE
        ref_match_id <- ref_paths[which(raw_sans == ref_sans)[1]]
      } else if (norm_in %in% ref_paths) {
        is_ref <- TRUE
        ref_match_id <- norm_in
      } else if (raw_in %in% ref_paths) {
        is_ref <- TRUE
        ref_match_id <- raw_in
      }
    } else if (!is.null(model$reference_image_path) && identical(norm_in, model$reference_image_path)) {
      is_ref <- TRUE
      ref_match_id <- norm_in
    }
  }
  if (!is_ref && !is.null(model$reference_image)) {
    ref_d <- dim(model$reference_image)
    if (identical(dims, ref_d)) {
      s_idx <- seq(1, length(im), length.out = min(100L, length(im)))
      if (isTRUE(all.equal(as.numeric(im)[s_idx], as.numeric(model$reference_image)[s_idx]))) {
        is_ref <- TRUE
      }
    }
  }

  if (is_ref) {
    b_match <- if (!is.null(ref_match_id) && !is.null(model$reference_boxes_map[[ref_match_id]])) {
      model$reference_boxes_map[[ref_match_id]]
    } else {
      model$detected_boxes
    }
    p_match <- if (!is.null(ref_match_id) && !is.null(model$reference_polys_map[[ref_match_id]])) {
      model$reference_polys_map[[ref_match_id]]
    } else {
      model$detected_polygons
    }

    if (!is.null(b_match) && nrow(b_match) > 0L) {
      if (!is.null(conf_t) && !is.null(b_match$conf)) {
        keep_ref <- which(b_match$conf >= conf_t)
        b_match <- b_match[keep_ref, , drop = FALSE]
        if (length(p_match) > 0L) {
          p_match <- p_match[keep_ref]
        }
      }
      res <- list(boxes = b_match, polygons = p_match, image = im)
      if (isTRUE(plot)) {
        plot(im)
        n_inst <- if (task == "segment" && length(res$polygons) > 0L) length(res$polygons) else nrow(res$boxes)
        if (isTRUE(rainbow) && n_inst > 0L) {
          c_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
          c_poly <- grDevices::adjustcolor(c_border, alpha.f = 0.40)
        } else {
          c_border <- rep("#00FFCC", max(1L, n_inst))
          c_poly <- rep("#00FFCC55", max(1L, n_inst))
        }
        if (task == "segment" && length(res$polygons) > 0L) {
          for (i in seq_along(res$polygons)) {
            poly <- res$polygons[[i]]
            if (!is.null(poly) && length(poly) >= 6L) {
              px <- poly[c(TRUE, FALSE)] * w
              py <- poly[c(FALSE, TRUE)] * h
              graphics::polygon(px, py, col = c_poly[i], border = c_border[i], lwd = 1.5)
            }
          }
        } else if (nrow(res$boxes) > 0L) {
          df_b <- res$boxes
          for (r in seq_len(nrow(df_b))) {
            graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = c_border[r], lwd = 2)
            graphics::text(df_b$xmin[r], df_b$ymin[r] - 5, labels = sprintf("%s %.2f", df_b$class_name[r], df_b$conf[r]),
                           col = c_border[r], cex = 0.8, pos = 4)
          }
        }
      }
      return(invisible(res))
    }
  }

  df_b <- data.frame()
  cand_polys <- list()

  # 2. Run PerSAM visual prototype transfer if SAM method
  if (identical(model$method, "sam") && !is.null(model$prototypes)) {
    num_arr <- image_data(im, type = "numeric")
    persam_res <- tryCatch({
      .run_persam(
        mat = num_arr,
        exemplar_points = NULL,
        precomputed_prototypes = model$prototypes,
        sim_threshold = conf_t,
        min_dist = if (!is.null(model$median_radius)) min(12.0, max(6.0, model$median_radius * 0.25)) else 10.0,
        max_objects = if (!is.null(max_objects)) max_objects else model$max_objects,
        feat_res = 256L,
        engine = if (!is.null(engine)) engine else if (!is.null(model$engine)) model$engine else "cpu",
        device_id = if (!is.null(device_id) && device_id >= 0) device_id else if (!is.null(model$device_id)) model$device_id else -1L,
        fill_hull = TRUE,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(persam_res) && nrow(persam_res$boxes) > 0L) {
      raw_boxes <- persam_res$boxes
      raw_conts <- persam_res$contours

      valid_idx <- seq_len(nrow(raw_boxes))
      if (!is.null(model$min_area) && model$min_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- which(b_areas >= model$min_area)
      }
      if (!is.null(model$max_area) && model$max_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- intersect(valid_idx, which(b_areas <= model$max_area))
      }

      if (length(valid_idx) > 0L) {
        raw_boxes <- raw_boxes[valid_idx, , drop = FALSE]
        raw_conts <- raw_conts[valid_idx]

        df_b <- data.frame(
          xmin = as.numeric(raw_boxes$xmin),
          ymin = as.numeric(raw_boxes$ymin),
          xmax = as.numeric(raw_boxes$xmax),
          ymax = as.numeric(raw_boxes$ymax),
          conf = round(as.numeric(raw_boxes$score), 3),
          class_id = if (!is.null(model$class_id)) as.integer(model$class_id) else 0L,
          class_name = model$label[1],
          stringsAsFactors = FALSE
        )

        for (i in seq_along(raw_conts)) {
          cnt <- raw_conts[[i]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            gx <- cnt[, 1] / w
            gy <- cnt[, 2] / h
            if (length(gx) > 32L) {
              idx <- round(seq(1, length(gx), length.out = 32))
              gx <- gx[idx]; gy <- gy[idx]
            }
            cand_polys[[i]] <- as.vector(rbind(gx, gy))
          } else {
            cand_polys[[i]] <- NULL
          }
        }
      }
    }
  }

  # 3. Fallback to adaptive morphology if PerSAM did not find instances or method != "sam"
  if (nrow(df_b) == 0L) {
    norm_arr <- image_data(im, type = "normalized")
    dr <- (norm_arr[, , 1] - model$target_rgb[1]) / model$rgb_sd[1]
    dg <- (norm_arr[, , 2] - model$target_rgb[2]) / model$rgb_sd[2]
    db <- (norm_arr[, , 3] - model$target_rgb[3]) / model$rgb_sd[3]
    dist_sq <- (dr^2 + dg^2 + db^2) / 3.0
    sim_raw <- exp(-0.5 * dist_sq)

    blur_sigma <- max(2, as.integer(round(model$median_radius * 0.12)))
    sim_im <- as_image(matrix(sim_raw, w, h), colormode = "Grayscale")
    sim_blur <- image_blur(sim_im, sigma = blur_sigma)
    sim_mat <- image_data(sim_blur, type = "normalized")

    adapt_t <- if (conf_t > 0.4) 0.25 else conf_t
    bin_mask <- sim_mat >= adapt_t
    lbl <- bwlabel_cpp(bin_mask)
    counts <- table(lbl[lbl > 0])
    valid <- as.integer(names(counts[counts >= model$min_area & counts <= model$max_area]))

    cand_boxes <- list()
    cand_confs <- numeric()
    conts_all <- if (length(valid) > 0L) extract_contours_cpp(lbl) else list()
    cand_p <- list()

    for (obj in valid) {
      coords <- which(lbl == obj, arr.ind = TRUE)
      bx1 <- min(coords[, 1]); bx2 <- max(coords[, 1])
      by1 <- min(coords[, 2]); by2 <- max(coords[, 2])
      bw <- bx2 - bx1 + 1L; bh <- by2 - by1 + 1L
      ar <- bw / pmax(1L, bh)
      if (ar < 0.38 || ar > 2.60) next

      obj_conf <- mean(sim_mat[coords])
      if (obj_conf < adapt_t) next

      poly_norm <- NULL
      if (length(conts_all) >= obj && is.matrix(conts_all[[obj]]) && nrow(conts_all[[obj]]) >= 3L) {
        poly_pts <- conts_all[[obj]]
        gx <- poly_pts[, 1] / w; gy <- poly_pts[, 2] / h
        if (length(gx) > 32L) {
          idx <- round(seq(1, length(gx), length.out = 32))
          gx <- gx[idx]; gy <- gy[idx]
        }
        poly_norm <- as.vector(rbind(gx, gy))
      }

      cand_boxes[[length(cand_boxes) + 1L]] <- c(bx1, by1, bx2, by2)
      cand_confs <- c(cand_confs, obj_conf)
      cand_p[[length(cand_p) + 1L]] <- poly_norm
    }

    if (length(cand_boxes) > 0L) {
      box_mat <- do.call(rbind, cand_boxes)
      if (nrow(box_mat) > 1L) {
        keep_idx <- .nms_boxes(box_mat, cand_confs, iou_thresh = 0.35)
        box_mat <- box_mat[keep_idx, , drop = FALSE]
        cand_confs <- cand_confs[keep_idx]
        cand_p <- cand_p[keep_idx]
      }
      df_b <- data.frame(
        xmin = as.numeric(box_mat[, 1]),
        ymin = as.numeric(box_mat[, 2]),
        xmax = as.numeric(box_mat[, 3]),
        ymax = as.numeric(box_mat[, 4]),
        conf = round(cand_confs, 3),
        class_id = if (!is.null(model$class_id)) as.integer(model$class_id) else 0L,
        class_name = model$label[1],
        stringsAsFactors = FALSE
      )
      cand_polys <- cand_p
    }
  }

  if (nrow(df_b) == 0L) {
    df_b <- data.frame(xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
                       conf = numeric(0), class_id = integer(0), class_name = character(0), stringsAsFactors = FALSE)
  }

  res <- list(boxes = df_b, polygons = cand_polys, image = im)

  if (isTRUE(plot)) {
    plot(im)
    n_inst <- if (task == "segment" && length(cand_polys) > 0L) length(cand_polys) else nrow(df_b)
    if (isTRUE(rainbow) && n_inst > 0L) {
      c_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
      c_poly <- grDevices::adjustcolor(c_border, alpha.f = 0.40)
    } else {
      c_border <- rep("#00FFCC", max(1L, n_inst))
      c_poly <- rep("#00FFCC55", max(1L, n_inst))
    }
    if (task == "segment" && length(cand_polys) > 0L) {
      for (i in seq_along(cand_polys)) {
        poly <- cand_polys[[i]]
        if (!is.null(poly) && length(poly) >= 6L) {
          px <- poly[c(TRUE, FALSE)] * w
          py <- poly[c(FALSE, TRUE)] * h
          graphics::polygon(px, py, col = c_poly[i], border = c_border[i], lwd = 1.5)
        }
      }
    } else if (nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = c_border[r], lwd = 2)
      }
    }

    if (nrow(df_b) > 0L && (isTRUE(show_class) || isTRUE(show_conf))) {
      for (r in seq_len(nrow(df_b))) {
        cb <- c_border[((r - 1L) %% length(c_border)) + 1L]
        lbl_parts <- character(0)
        if (isTRUE(show_class)) {
          cur_lbl <- if ("class_name" %in% names(df_b) && nzchar(df_b$class_name[r])) df_b$class_name[r] else model$label[1]
          if (nzchar(cur_lbl)) lbl_parts <- c(lbl_parts, cur_lbl)
        }
        if (isTRUE(show_conf) && !is.null(df_b$conf) && length(df_b$conf) >= r && !is.na(df_b$conf[r])) {
          lbl_parts <- c(lbl_parts, sprintf("%.2f", df_b$conf[r]))
        }
        tag <- paste(lbl_parts, collapse = " ")

        if (nzchar(tag)) {
          tw <- graphics::strwidth(tag, cex = 0.8)
          th <- graphics::strheight(tag, cex = 0.8)
          y_top <- max(th + 6, df_b$ymin[r])
          if (isTRUE(badge)) {
            bg_tag <- grDevices::adjustcolor(cb, alpha.f = 0.90)
            graphics::rect(df_b$xmin[r], y_top - th - 6, df_b$xmin[r] + tw + 6, y_top, col = bg_tag, border = NA)
            graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = "#FFFFFF", cex = 0.8, font = 2, adj = c(0, 1))
          } else {
            graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = cb, cex = 0.8, font = 2, adj = c(0, 1))
          }
        }
      }
    }
  }

  invisible(res)
}

# Vectorized Non-Maximum Suppression (NMS) helper
.nms_boxes <- function(boxes, scores, iou_thresh = 0.35) {
  if (nrow(boxes) <= 1L) return(seq_len(nrow(boxes)))
  order_idx <- order(scores, decreasing = TRUE)
  keep <- integer()

  while (length(order_idx) > 0L) {
    curr <- order_idx[1]
    keep <- c(keep, curr)
    if (length(order_idx) == 1L) break

    rest <- order_idx[-1]
    xx1 <- pmax(boxes[curr, 1], boxes[rest, 1])
    yy1 <- pmax(boxes[curr, 2], boxes[rest, 2])
    xx2 <- pmin(boxes[curr, 3], boxes[rest, 3])
    yy2 <- pmin(boxes[curr, 4], boxes[rest, 4])

    w_inter <- pmax(0, xx2 - xx1 + 1)
    h_inter <- pmax(0, yy2 - yy1 + 1)
    inter <- w_inter * h_inter

    area_curr <- (boxes[curr, 3] - boxes[curr, 1] + 1) * (boxes[curr, 4] - boxes[curr, 2] + 1)
    area_rest <- (boxes[rest, 3] - boxes[rest, 1] + 1) * (boxes[rest, 4] - boxes[rest, 2] + 1)
    iou <- inter / (area_curr + area_rest - inter)

    order_idx <- rest[iou <= iou_thresh]
  }
  return(keep)
}

#' @title Auto-Label an Entire Directory into a Production YOLO Dataset
#' @name yolo_dataset_autolabel
#' @description
#' `yolo_dataset_autolabel()` batch-processes an entire folder of unannotated images using
#' a trained `yolo_fewshot_model` (from [yolo_fewshot_fit()]). It predicts bounding boxes
#' or instance segmentation polygons, writes standard YOLO `.txt` annotation files, splits
#' the data into `train`, `val`, and `test`, and generates a `data.yaml` descriptor.
#'
#' It also creates a `review/` visual inspection folder with rendered boxes/polygons so the user
#' can quickly audit the automated labels.
#'
#' @param input_dir Directory containing raw unannotated images.
#' @param model A `yolo_fewshot_model` object fitted with [yolo_fewshot_fit()].
#' @param output_dir Target directory to write the YOLO dataset. If `NULL` (default),
#'   the dataset structure (`images/`, `labels/`, `data.yaml`) is generated directly inside `input_dir`.
#'   Specify a new directory path only if you wish to create a separate dataset folder.
#' @param task Target task: `"segment"` (polygons) or `"detect"` (bounding boxes). Defaults to `model$task`.
#' @param classes Character vector of class names. If `NULL` (default), automatically inferred from the model object or metadata.
#' @param conf_threshold Minimum confidence / similarity threshold for accepting predicted instances (e.g. 0.50). If `NULL` (default), uses `model$conf_threshold`.
#' @param max_objects Maximum number of objects to detect per image. If `NULL` (default), no limit is applied.
#' @param split_ratio Numeric vector of length 3 for `c(train, val, test)` split. Defaults to `c(0.70, 0.20, 0.10)`.
#' @param crop_to_objects Logical. If `TRUE`, crops the image to the bounding box containing all detected objects plus a margin (default: `FALSE`).
#' @param crop_margin Numeric. Margin fraction added around the object bounding box when `crop_to_objects = TRUE` (default: `0.05`, i.e., 5%).
#' @param engine Execution engine override: `"cpu"` or `"gpu"`. If `NULL`, uses `model$engine`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param move Logical. If `TRUE` (default), moves the source images to the target split folders
#'   instead of copying them, saving disk space and avoiding duplicate files. An undo operation can be performed via [yolo_dataset_undo()]
#'   to restore files to their original location. Defaults to `TRUE`.
#' @param review Logical or integer. If `TRUE` (default), writes annotated review images to `review/`
#'   subfolder for visual inspection of all images. An integer can be provided (e.g. `review = 30`)
#'   to limit the number of review images generated. Setting `review = FALSE` avoids creating duplicate preview images.
#' @param rainbow Logical. If `TRUE`, review previews are rendered with distinct rainbow colors per instance. Defaults to `TRUE`.
#' @param show_class Logical. Whether to display class name on labels/badges in review images. Defaults to `TRUE`.
#' @param show_conf Logical. Whether to display confidence score on labels/badges in review images. Defaults to `TRUE`.
#' @param badge Logical. Whether to draw a filled background badge behind labels in review images. Defaults to `TRUE`.
#' @param seed Random seed for reproducible splitting. Defaults to `42`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#' @param ... Additional arguments passed to methods or aliases (e.g. `show_labels`, `show_scores`, `label_box`).
#'
#' @return A `yolo_dataset` S3 object pointing to the created dataset.
#' @export
yolo_dataset_autolabel <- function(input_dir,
                                   model,
                                   output_dir = NULL,
                                   task = c("detect", "segment"),
                                   classes = NULL,
                                   conf_threshold = NULL,
                                   max_objects = NULL,
                                   split_ratio = c(0.70, 0.20, 0.10),
                                   crop_to_objects = FALSE,
                                   crop_margin = 0.05,
                                   engine = "gpu",
                                   device_id = -1,
                                   move = TRUE,
                                   review = TRUE,
                                   rainbow = TRUE,
                                   show_class = TRUE,
                                   show_conf = TRUE,
                                   badge = TRUE,
                                   seed = 42,
                                   verbose = TRUE,
                                   ...) {
  dots <- list(...)
  if ("move" %in% names(dots)) move <- isTRUE(dots$move)
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)
  if ("crop_to_objects" %in% names(dots)) crop_to_objects <- isTRUE(dots$crop_to_objects)
  if ("crop_margin" %in% names(dots)) crop_margin <- as.numeric(dots$crop_margin[1])
  if (is.numeric(crop_to_objects) && length(crop_to_objects) == 1L && crop_to_objects > 0) {
    crop_margin <- as.numeric(crop_to_objects)
    crop_to_objects <- TRUE
  }

  # Support swapping input_dir and model if called as yolo_dataset_autolabel(model, input_dir, ...)
  if ((inherits(input_dir, c("yolo_workflow", "yolo_train", "yolo_fewshot_model")) ||
       (is.character(input_dir) && length(input_dir) == 1L && grepl("\\.(onnx|pt)$", input_dir))) &&
      (missing(model) || (is.character(model) && length(model) == 1L && dir.exists(model)))) {
    tmp_m <- input_dir
    input_dir <- if (!missing(model)) model else if (!is.null(tmp_m$workflow$data)) tmp_m$workflow$data else "."
    model <- tmp_m
  }

  if (missing(task) || is.null(task)) {
    if (inherits(model, c("yolo_workflow", "yolo_train")) && !is.null(model$task)) {
      task <- model$task
    } else if (inherits(model, "yolo_fewshot_model") && !is.null(model$task)) {
      task <- model$task
    } else {
      task <- "detect"
    }
  } else {
    task <- match.arg(task, c("detect", "segment"))
  }

  is_yolo_model <- FALSE
  yolo_model_path <- NULL
  if (inherits(model, c("yolo_workflow", "yolo_train"))) {
    is_yolo_model <- TRUE
    yolo_model_path <- model$model_file
    if ((is.null(yolo_model_path) || !file.exists(yolo_model_path)) && !is.null(model$pt_file) && file.exists(model$pt_file)) {
      yolo_model_path <- model$pt_file
    }
    if (!is.null(yolo_model_path) && file.exists(yolo_model_path)) {
      yolo_model_path <- normalizePath(yolo_model_path, winslash = "/", mustWork = FALSE)
    }
    if (missing(task) || is.null(task)) {
      task <- if (!is.null(model$task)) model$task else "detect"
    }
    if (is.null(classes)) {
      classes <- if (!is.null(model$classes)) as.character(model$classes) else "object"
    }
  } else if (is.character(model) && length(model) == 1L) {
    # Check if it's an RDS file containing a model object
    if (grepl("\\.rds$", model, ignore.case = TRUE) && file.exists(model)) {
      loaded_m <- try(readRDS(model), silent = TRUE)
      if (!inherits(loaded_m, "try-error")) {
        model <- loaded_m
      }
    }
    if (inherits(model, c("yolo_workflow", "yolo_train"))) {
      is_yolo_model <- TRUE
      yolo_model_path <- model$model_file
      if ((is.null(yolo_model_path) || !file.exists(yolo_model_path)) && !is.null(model$pt_file) && file.exists(model$pt_file)) {
        yolo_model_path <- model$pt_file
      }
      if (missing(task) || is.null(task)) {
        task <- if (!is.null(model$task)) model$task else "detect"
      }
      if (is.null(classes)) {
        classes <- if (!is.null(model$classes)) as.character(model$classes) else "object"
      }
    } else if (inherits(model, "yolo_fewshot_model")) {
      # Handled in fewshot block below
    } else {
      stem <- sub("\\.(onnx|pt)$", "", basename(model), ignore.case = TRUE)
      cand_paths <- c(
        model,
        file.path(getwd(), model),
        file.path(pliman_model_dir(), model),
        file.path(pliman_model_dir(), paste0(model, ".onnx")),
        file.path(pliman_model_dir(), paste0(model, ".pt")),
        file.path(pliman_model_dir(), paste0(stem, ".onnx")),
        file.path(pliman_model_dir(), paste0(stem, ".pt")),
        file.path(pliman_model_dir(), paste0(stem, "_results"), "weights", "best.onnx"),
        file.path(pliman_model_dir(), paste0(stem, "_results"), "weights", "best.pt"),
        file.path(getwd(), paste0(stem, "_results"), "weights", "best.onnx"),
        file.path(getwd(), paste0(stem, "_results"), "weights", "best.pt"),
        file.path(dirname(model), paste0(stem, "_results"), "weights", "best.onnx"),
        file.path(dirname(model), paste0(stem, "_results"), "weights", "best.pt")
      )
      cand_paths <- unique(cand_paths[nzchar(cand_paths)])
      ex_cand <- cand_paths[file.exists(cand_paths) & !dir.exists(cand_paths)]
      if (length(ex_cand) > 0L) {
        is_yolo_model <- TRUE
        yolo_model_path <- normalizePath(ex_cand[1], winslash = "/")
      } else {
        loc_m <- .find_local_model(model, dir = pliman_model_dir())
        if (!is.null(loc_m) && file.exists(loc_m)) {
          is_yolo_model <- TRUE
          yolo_model_path <- normalizePath(loc_m, winslash = "/")
        } else if (grepl("\\.(onnx|pt)$", model, ignore.case = TRUE) || grepl("yolo", model, ignore.case = TRUE)) {
          is_yolo_model <- TRUE
          yolo_model_path <- model
        } else {
          cli::cli_abort(c(
            "Model {.val {model}} was not found.",
            "i" = "Searched in working directory: {.path {getwd()}}",
            "i" = "Searched in pliman models directory: {.path {pliman_model_dir()}}"
          ))
        }
      }
      if (is_yolo_model) {
        if (is.null(classes)) {
          classes <- "object"
          model_stem <- sub("\\.(onnx|pt)$", "", basename(yolo_model_path), ignore.case = TRUE)
          res_dirs <- c(
            file.path(dirname(yolo_model_path), paste0(model_stem, "_results")),
            file.path(pliman_model_dir(), paste0(model_stem, "_results")),
            sub("/weights$", "", dirname(yolo_model_path))
          )
          res_dirs <- unique(res_dirs[dir.exists(res_dirs)])
          if (length(res_dirs) > 0L) {
            pr_f <- file.path(res_dirs[1], "pr_curve.csv")
            cm_f <- file.path(res_dirs[1], "confusion_matrix.csv")
            if (file.exists(pr_f)) {
              pr_hdr <- names(utils::read.csv(pr_f, nrows = 1, check.names = FALSE))
              if (length(pr_hdr) > 1L) classes <- pr_hdr[-1]
            } else if (file.exists(cm_f)) {
              cm_hdr <- names(utils::read.csv(cm_f, nrows = 1, check.names = FALSE))
              cm_classes <- cm_hdr[-1]
              cm_classes <- cm_classes[!tolower(cm_classes) %in% c("background", "fundo")]
              if (length(cm_classes) > 0L) classes <- cm_classes
            }
          }
        }
      }
    }
  }

  if (!is_yolo_model) {
    if (is.list(model) && !inherits(model, "yolo_fewshot_model")) {
      if (all(vapply(model, function(m) inherits(m, "yolo_fewshot_model"), logical(1)))) {
        model <- do.call(c, model)
      }
    }
    if (!inherits(model, "yolo_fewshot_model")) {
      cli::cli_abort("{.arg model} must be a trained YOLO model, file path/name, or a {.cls yolo_fewshot_model} from {.fn yolo_fewshot_fit}.")
    }
    if (is.null(classes)) {
      classes <- if (!is.null(model$classes)) {
        as.character(model$classes)
      } else if (!is.null(model$label)) {
        as.character(model$label)
      } else {
        "object"
      }
    }
  }

  if (!dir.exists(input_dir)) {
    cli::cli_abort("Directory {.path {input_dir}} does not exist.")
  }

  fls <- list.files(input_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE, recursive = FALSE)
  fls <- fls[!grepl("^review_", basename(fls))]
  if (length(fls) == 0L) {
    cand_fls <- list.files(input_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
    cand_fls <- cand_fls[!grepl("(_aug\\d+|_tile_\\d+_\\d+|^review_)", basename(cand_fls))]
    fls <- cand_fls
  }
  if (length(fls) == 0L) {
    cli::cli_abort("No valid images found in {.path {input_dir}}.")
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- input_dir
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  splits <- c("train", "val", "test")
  for (sp in splits) {
    dir.create(file.path(output_dir, "images", sp), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "labels", sp), recursive = TRUE, showWarnings = FALSE)
  }

  do_review <- isTRUE(review) || (is.numeric(review) && length(review) == 1L && review[1] > 0)
  max_review <- if (is.numeric(review) && length(review) == 1L && review[1] > 0) as.integer(review[1]) else if (isTRUE(review)) Inf else 0L

  review_dir <- file.path(output_dir, "review")
  if (do_review) dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  # Stratify images
  n_total <- length(fls)
  shuffled <- sample(fls)
  n_tr <- max(1L, round(n_total * split_ratio[1]))
  n_val <- max(1L, round(n_total * split_ratio[2]))

  split_assign <- character(n_total)
  split_assign[seq_len(n_tr)] <- "train"
  split_assign[(n_tr + 1L):min(n_total, n_tr + n_val)] <- "val"
  if (n_total > (n_tr + n_val)) {
    split_assign[(n_tr + n_val + 1L):n_total] <- "test"
  } else {
    split_assign[split_assign == ""] <- "val"
  }

  used_eng <- if (!is.null(engine)) engine else if (is.list(model) && !is.null(model$engine)) model$engine else "cpu"
  conf_t <- if (!is.null(conf_threshold)) as.numeric(conf_threshold) else if (is.list(model) && !is.null(model$conf_threshold)) as.numeric(model$conf_threshold) else 0.40

  cli_pb <- NULL
  if (isTRUE(verbose)) {
    cli::cli_h2("Automated YOLO Dataset Labeling (pliman)")
    cls_disp <- paste(classes, collapse = ", ")
    cls_label <- if (length(classes) > 1L) "Classes" else "Class"
    cli::cli_alert_info("{cls_label}: {.val {cls_disp}} | Task: {.val {task}} | Conf: {.val {conf_t}} | Images: {.val {n_total}} | Engine: {.val {used_eng}}")
    cli_pb <- cli::cli_progress_bar(
      format = "{cli::pb_spin} Auto-labeling [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      total = n_total,
      type = "iterator",
      clear = FALSE
    )
  }

  # Ensure unique destination stems across images (prevents label overwrite if same stem with different extensions)
  stems_raw <- tools::file_path_sans_ext(basename(fls))
  exts_raw  <- tools::file_ext(basename(fls))

  unique_stems <- character(length(fls))
  names(unique_stems) <- fls
  unique_fns <- character(length(fls))
  names(unique_fns) <- fls

  if (any(duplicated(stems_raw))) {
    dup_stems <- unique(stems_raw[duplicated(stems_raw)])
    for (idx in seq_along(fls)) {
      f <- fls[idx]
      s <- stems_raw[idx]
      e <- exts_raw[idx]
      if (s %in% dup_stems) {
        stem_clean <- paste0(s, "_", tolower(e))
        unique_stems[f] <- stem_clean
        unique_fns[f]   <- paste0(stem_clean, ".", e)
      } else {
        unique_stems[f] <- s
        unique_fns[f]   <- basename(f)
      }
    }
  } else {
    for (idx in seq_along(fls)) {
      unique_stems[fls[idx]] <- stems_raw[idx]
      unique_fns[fls[idx]]   <- basename(fls[idx])
    }
  }

  total_objects <- 0L
  logged_orig_paths <- character()
  logged_dest_imgs <- character()
  logged_dest_lbls <- character()
  logged_rev_imgs <- character()

  for (i in seq_along(shuffled)) {
    im_path <- shuffled[i]
    sp <- split_assign[i]
    if (is_yolo_model) {
      model_dir_ctx <- if (!is.null(yolo_model_path) && file.exists(yolo_model_path)) dirname(yolo_model_path) else pliman_model_dir()
      if (identical(task, "detect")) {
        det <- try(image_detect_dl(im_path, model = yolo_model_path,
                                   labels = classes,
                                   conf_threshold = if (!is.null(conf_t)) conf_t else 0.40,
                                   plot = FALSE, verbose = FALSE,
                                   dir = model_dir_ctx,
                                   engine = used_eng, device_id = device_id), silent = TRUE)
        if (inherits(det, "try-error")) {
          if (isTRUE(verbose) && !is.null(cli_pb)) cli::cli_progress_update(id = cli_pb)
          next
        }
        df_b <- as.data.frame(det)
        polys <- list()
        im_obj <- image_import(im_path)
      } else {
        seg <- try(image_segment_dl(im_path, model = yolo_model_path,
                                    labels = classes,
                                    box_threshold = if (!is.null(conf_t)) conf_t else 0.40,
                                    plot = FALSE, verbose = FALSE,
                                    dir = model_dir_ctx,
                                    engine = used_eng, device_id = device_id), silent = TRUE)
        if (inherits(seg, "try-error")) {
          if (isTRUE(verbose) && !is.null(cli_pb)) cli::cli_progress_update(id = cli_pb)
          next
        }
        df_b <- if (!is.null(seg$boxes)) as.data.frame(seg$boxes) else if (!is.null(attr(seg, "boxes"))) as.data.frame(attr(seg, "boxes")) else data.frame()
        polys <- list()
        cnts <- if (!is.null(seg$contours)) seg$contours else attr(seg, "contours")
        if (!is.null(cnts)) {
          im_tmp <- image_import(im_path)
          iw <- dim(im_tmp)[1]; ih <- dim(im_tmp)[2]
          for (ci in seq_along(cnts)) {
            cnt <- cnts[[ci]]
            if (is.matrix(cnt) && nrow(cnt) >= 3L) {
              gx <- cnt[, 1] / iw; gy <- cnt[, 2] / ih
              if (length(gx) > 32L) {
                idx <- round(seq(1, length(gx), length.out = 32))
                gx <- gx[idx]; gy <- gy[idx]
              }
              polys[[ci]] <- as.vector(rbind(gx, gy))
            }
          }
        }
        im_obj <- image_import(im_path)
      }
      if (nrow(df_b) > 0L) {
        if (!"conf" %in% names(df_b)) {
          if ("score" %in% names(df_b)) df_b$conf <- df_b$score
          else if ("confidence" %in% names(df_b)) df_b$conf <- df_b$confidence
        }
        if (!"class_id" %in% names(df_b)) {
          if ("label" %in% names(df_b) && length(classes) > 0L) {
            df_b$class_id <- match(df_b$label, classes) - 1L
            df_b$class_id[is.na(df_b$class_id)] <- 0L
          } else {
            df_b$class_id <- 0L
          }
        }
        if (!"class_name" %in% names(df_b)) {
          if ("label" %in% names(df_b)) {
            df_b$class_name <- df_b$label
          } else if ("class_id" %in% names(df_b) && length(classes) > 0L) {
            c_safe <- df_b$class_id + 1L
            c_safe[c_safe < 1L | c_safe > length(classes)] <- 1L
            df_b$class_name <- classes[c_safe]
          }
        }
      }
    } else {
      pred <- try(yolo_fewshot_predict(model, im_path, task = task, conf_threshold = conf_t, max_objects = max_objects, engine = engine, device_id = device_id, plot = FALSE), silent = TRUE)
      if (inherits(pred, "try-error")) {
        if (isTRUE(verbose) && !is.null(cli_pb)) cli::cli_progress_update(id = cli_pb)
        next
      }

      df_b <- pred$boxes
      polys <- pred$polygons

      # Filter strictly by conf_t
      if (!is.null(conf_t) && nrow(df_b) > 0L && !is.null(df_b$conf)) {
        pass_idx <- which(df_b$conf >= conf_t)
        df_b <- df_b[pass_idx, , drop = FALSE]
        if (length(polys) > 0L) {
          polys <- polys[pass_idx]
        }
      }

      im_obj <- pred$image
    }
    dims <- dim(im_obj)
    w <- dims[1]; h <- dims[2]

    # Optional crop to objects
    was_cropped <- FALSE
    if (isTRUE(crop_to_objects)) {
      cmarg <- if (is.numeric(crop_margin) && length(crop_margin) >= 1L && !is.na(crop_margin[1])) as.numeric(crop_margin[1]) else 0.05
      c_res <- .yolo_crop_to_objects(im_obj, boxes = df_b, polygons = polys, margin = cmarg)
      if (isTRUE(c_res$cropped)) {
        im_obj <- c_res$img
        df_b <- c_res$boxes
        polys <- c_res$polygons
        w <- dim(im_obj)[1]
        h <- dim(im_obj)[2]
        was_cropped <- TRUE
      }
    }

    # Destination paths with unique stems
    fn_dest <- unique_fns[[im_path]]
    stem_dest <- unique_stems[[im_path]]
    dest_img <- file.path(output_dir, "images", sp, fn_dest)
    dest_lbl <- file.path(output_dir, "labels", sp, paste0(stem_dest, ".txt"))

    logged_orig_paths <- c(logged_orig_paths, im_path)
    logged_dest_imgs <- c(logged_dest_imgs, dest_img)
    logged_dest_lbls <- c(logged_dest_lbls, dest_lbl)

    if (was_cropped) {
      image_export(im_obj, name = fn_dest, subfolder = file.path(output_dir, "images", sp), verbose = FALSE)
      if (isTRUE(move) && normalizePath(im_path, winslash = "/", mustWork = FALSE) != normalizePath(dest_img, winslash = "/", mustWork = FALSE)) {
        unlink(im_path, force = TRUE)
      }
    } else if (isTRUE(move)) {
      if (normalizePath(im_path, winslash = "/", mustWork = FALSE) != normalizePath(dest_img, winslash = "/", mustWork = FALSE)) {
        mv_ok <- file.rename(im_path, dest_img)
        if (!mv_ok) {
          file.copy(im_path, dest_img, overwrite = TRUE)
          unlink(im_path, force = TRUE)
        }
      }
    } else {
      if (normalizePath(im_path, winslash = "/", mustWork = FALSE) != normalizePath(dest_img, winslash = "/", mustWork = FALSE)) {
        file.copy(im_path, dest_img, overwrite = TRUE)
      }
    }

    # Build YOLO label lines
    lines <- character()
    if (nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        cid <- if ("class_id" %in% names(df_b)) as.integer(df_b$class_id[r]) else 0L
        if (task == "detect") {
          xc <- ((df_b$xmin[r] + df_b$xmax[r]) / 2) / w
          yc <- ((df_b$ymin[r] + df_b$ymax[r]) / 2) / h
          bw <- (df_b$xmax[r] - df_b$xmin[r]) / w
          bh <- (df_b$ymax[r] - df_b$ymin[r]) / h
          lines <- c(lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, bw, bh))
        } else {
          poly <- polys[[r]]
          if (!is.null(poly) && length(poly) >= 6L) {
            lines <- c(lines, paste(c(as.character(cid), sprintf("%.6f", poly)), collapse = " "))
          } else {
            # Bounding box fallback as 4-corner polygon
            xc1 <- df_b$xmin[r] / w; yc1 <- df_b$ymin[r] / h
            xc2 <- df_b$xmax[r] / w; yc2 <- df_b$ymax[r] / h
            lines <- c(lines, sprintf("%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f",
                                      cid, xc1, yc1, xc2, yc1, xc2, yc2, xc1, yc2))
          }
        }
      }
      total_objects <- total_objects + nrow(df_b)
    }

    writeLines(lines, dest_lbl)

    # Optional review image rendering (wrapped in tryCatch to prevent file locks or device errors from crashing autolabel)
    cur_rev <- NA_character_
    if (do_review && i <= max_review) {
      tryCatch({
        rev_out <- file.path(review_dir, paste0("review_", fn_dest))
        out_w <- 1200L
        out_h <- max(100L, as.integer(round(1200 * h / w)))
        dev_opened <- FALSE
        tryCatch({
          grDevices::png(rev_out, width = out_w, height = out_h, res = 120)
          dev_opened <- TRUE
          cur_rev <- rev_out
        }, error = function(e) {
          # Try alternative filename if primary is locked by image viewer
          alt_out <- file.path(review_dir, paste0("review_", i, "_", fn_dest))
          tryCatch({
            grDevices::png(alt_out, width = out_w, height = out_h, res = 120)
            dev_opened <<- TRUE
            cur_rev <<- alt_out
          }, error = function(e2) NULL)
        })

        if (dev_opened) {
          plot(im_obj)
          n_b <- nrow(df_b)
          if (isTRUE(rainbow) && n_b > 0L) {
            rev_border <- grDevices::rainbow(n_b, s = 0.85, v = 0.95)
            rev_poly <- grDevices::adjustcolor(rev_border, alpha.f = 0.40)
          } else {
            rev_border <- rep("#00FFCC", max(1L, n_b))
            rev_poly <- rep("#00FFCC55", max(1L, n_b))
          }
          if (task == "segment" && length(polys) > 0L) {
            for (r in seq_len(n_b)) {
              poly <- polys[[r]]
              if (!is.null(poly) && length(poly) >= 6L) {
                px <- poly[c(TRUE, FALSE)] * w
                py <- poly[c(FALSE, TRUE)] * h
                graphics::polygon(px, py, col = rev_poly[r], border = rev_border[r], lwd = 1.5)
              } else {
                graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = rev_border[r], lwd = 2)
              }
            }
          } else if (n_b > 0L) {
            for (r in seq_len(n_b)) {
              graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = rev_border[r], lwd = 2)
            }
          }

          if (n_b > 0L && (isTRUE(show_class) || isTRUE(show_conf))) {
            for (r in seq_len(n_b)) {
              lbl_parts <- character(0)
              if (isTRUE(show_class)) {
                cur_lbl <- if ("class_name" %in% names(df_b) && nzchar(df_b$class_name[r])) df_b$class_name[r] else classes[1]
                if (nzchar(cur_lbl)) lbl_parts <- c(lbl_parts, cur_lbl)
              }
              if (isTRUE(show_conf) && !is.null(df_b$conf) && !is.na(df_b$conf[r])) {
                lbl_parts <- c(lbl_parts, sprintf("%.2f", df_b$conf[r]))
              }
              tag <- paste(lbl_parts, collapse = " ")

              if (nzchar(tag)) {
                tw <- graphics::strwidth(tag, cex = 0.8)
                th <- graphics::strheight(tag, cex = 0.8)
                y_top <- max(th + 6, df_b$ymin[r])
                if (isTRUE(badge)) {
                  bg_tag <- grDevices::adjustcolor(rev_border[r], alpha.f = 0.90)
                  graphics::rect(df_b$xmin[r], y_top - th - 6, df_b$xmin[r] + tw + 6, y_top, col = bg_tag, border = NA)
                  graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = "#FFFFFF", cex = 0.8, font = 2, adj = c(0, 1))
                } else {
                  graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = rev_border[r], cex = 0.8, font = 2, adj = c(0, 1))
                }
              }
            }
          }
          graphics::title(sprintf("%s: %d %s", fn_dest, n_b, if (task == "segment") "segmented" else "detected"),
                          col.main = if (isTRUE(rainbow)) "#222222" else "#00FFCC",
                          cex.main = 1.0)
          grDevices::dev.off()
        }
      }, error = function(e) {
        # Ensure graphics device is closed if an unexpected error occurs during plotting
        if (grDevices::dev.cur() > 1L) {
          tryCatch(grDevices::dev.off(), error = function(e) NULL)
        }
      })
    }
    logged_rev_imgs <- c(logged_rev_imgs, cur_rev)

    if (isTRUE(verbose) && !is.null(cli_pb)) {
      cli::cli_progress_update(id = cli_pb)
    }
  }

  if (isTRUE(verbose) && !is.null(cli_pb)) {
    cli::cli_progress_done(id = cli_pb)
  }

  # Write data.yaml
  yaml_content <- c(
    paste0("path: ", output_dir),
    "train: images/train",
    "val: images/val",
    "test: images/test",
    "",
    "names:",
    paste0("  ", seq_along(classes) - 1L, ": ", classes)
  )
  writeLines(yaml_content, file.path(output_dir, "data.yaml"))

  # Save undo history manifest (.yolo_history.rds)
  hist_file <- file.path(output_dir, ".yolo_history.rds")
  existing_hist <- if (file.exists(hist_file)) tryCatch(readRDS(hist_file), error = function(e) list()) else list()
  if (!is.list(existing_hist)) existing_hist <- list()

  history_entry <- list(
    id = paste0("autolabel_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    timestamp = Sys.time(),
    operation = "autolabel",
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = FALSE),
    output_dir = output_dir,
    moved = isTRUE(move),
    images = data.frame(
      orig = logged_orig_paths,
      dest = logged_dest_imgs,
      lbl = logged_dest_lbls,
      rev = logged_rev_imgs,
      stringsAsFactors = FALSE
    ),
    yaml_file = file.path(output_dir, "data.yaml"),
    review_dir = if (do_review) review_dir else NULL
  )
  existing_hist[[length(existing_hist) + 1L]] <- history_entry
  saveRDS(existing_hist, hist_file)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Dataset auto-labeled successfully!")
    cli::cli_alert_info("Total instances generated: {.val {total_objects}}")
    if (do_review) {
      cli::cli_alert_info("Review previews written to: {.path {review_dir}}")
    }
  }

  # Ensure no unwanted 'runs' directory remains from any model initialization
  wd_runs <- file.path(getwd(), "runs")
  if (dir.exists(wd_runs)) {
    unlink(wd_runs, recursive = TRUE, force = TRUE)
  }
  if (!is.null(output_dir)) {
    out_runs <- file.path(output_dir, "runs")
    if (dir.exists(out_runs)) {
      unlink(out_runs, recursive = TRUE, force = TRUE)
    }
  }
  if (!is.null(input_dir) && is.character(input_dir) && length(input_dir) == 1L && dir.exists(input_dir)) {
    in_runs <- file.path(input_dir, "runs")
    if (dir.exists(in_runs)) {
      unlink(in_runs, recursive = TRUE, force = TRUE)
    }
  }

  .yolo_wrap_dataset(output_dir, task = task, classes = classes)
}

#' @title Plot a Few-Shot Exemplar Model
#' @name plot.yolo_fewshot_model
#' @description
#' Plots reference image(s) with bounding boxes, segmentation polygons, and clicked exemplar points.
#' When the model contains multiple reference images, they are arranged automatically in a multi-panel grid.
#' Setting \code{rainbow = TRUE} renders each detected instance in distinct vibrant colors.
#'
#' @param x A \code{yolo_fewshot_model} object.
#' @param task Target task: \code{"segment"} or \code{"detect"}. Defaults to \code{x$task}.
#' @param rainbow Logical. If \code{TRUE}, colors each detected instance with distinct vibrant hues. Defaults to \code{FALSE}.
#' @param show_boxes Logical. Show bounding boxes. Defaults to \code{TRUE} for detection tasks.
#' @param show_polygons Logical. Show segmentation polygons. Defaults to \code{TRUE} for segmentation tasks.
#' @param show_points Logical. Show clicked exemplar points. Defaults to \code{TRUE}.
#' @param show_class Logical. Whether to display class name on labels/badges. Defaults to \code{TRUE}.
#' @param show_conf Logical. Whether to display confidence score on labels/badges. Defaults to \code{TRUE}.
#' @param badge Logical. Whether to draw a filled background badge behind labels. Defaults to \code{TRUE}.
#' @param col_border Border color when \code{rainbow = FALSE}. Defaults to \code{"#00FFCC"}.
#' @param col_poly Fill color for polygons when \code{rainbow = FALSE}. Defaults to \code{"#00FFCC55"}.
#' @param col_point Point color for exemplars. Defaults to \code{"#FFFF00"}.
#' @param lwd Line width for borders. Defaults to \code{2}.
#' @param mfrow Numeric vector of length 2 for grid layout (e.g. \code{c(2, 2)}). If \code{NULL}, calculated automatically.
#' @param which Indices or names of reference images to plot. If \code{NULL} (default), plots all reference images.
#' @param ... Additional arguments passed to graphics methods.
#'
#' @return Invisible \code{x}.
#' @export
plot.yolo_fewshot_model <- function(x,
                                    task = x$task,
                                    rainbow = TRUE,
                                    show_boxes = NULL,
                                    show_polygons = NULL,
                                    show_points = TRUE,
                                    show_class = if (!is.null(x$show_class)) x$show_class else TRUE,
                                    show_conf = if (!is.null(x$show_conf)) x$show_conf else TRUE,
                                    badge = if (!is.null(x$badge)) x$badge else TRUE,
                                    col_border = "#00FFCC",
                                    col_poly = "#00FFCC55",
                                    col_point = "#FFFF00",
                                    lwd = 2,
                                    mfrow = NULL,
                                    which = NULL,
                                    ...) {
  dots <- list(...)
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  if (inherits(x, "yolo_fewshot_multimodel") && !is.null(x$models)) {
    im_ref <- if (!is.null(x$reference_image)) x$reference_image else x$models[[1]]$reference_image
    if (!is.null(im_ref)) {
      yolo_fewshot_predict(
        model = x,
        img = im_ref,
        task = task,
        plot = TRUE,
        rainbow = rainbow,
        show_class = show_class,
        show_conf = show_conf,
        badge = badge,
        conf_threshold = x$conf_threshold,
        engine = x$engine
      )
    }
    return(invisible(x))
  }
  imgs <- if (!is.null(x$reference_images) && length(x$reference_images) > 0L) {
    x$reference_images
  } else if (!is.null(x$reference_image)) {
    list(x$reference_image)
  } else {
    cli::cli_abort("No reference images found in the model to plot.")
  }

  if (!is.null(names(imgs))) {
    # If duplicate entries exist (e.g. raw stem and full path), deduplicate by stem
    stems <- tolower(tools::file_path_sans_ext(basename(names(imgs))))
    if (anyDuplicated(stems)) {
      has_ext <- grepl("\\.[a-zA-Z0-9]+$", names(imgs))
      ord <- order(!has_ext)
      stems_ord <- stems[ord]
      keep_ord <- !duplicated(stems_ord)
      keep_idx <- ord[keep_ord]
      imgs <- imgs[sort(keep_idx)]
    }
  }

  if (!is.null(which)) {
    imgs <- imgs[which]
  }

  n_imgs <- length(imgs)
  if (n_imgs == 0L) {
    cli::cli_abort("Selected index in {.arg which} does not match any reference image.")
  }

  task <- if (!is.null(task)) match.arg(task, c("segment", "detect")) else if (!is.null(x$task)) x$task else "segment"
  if (is.null(show_polygons)) show_polygons <- (task == "segment")
  if (is.null(show_boxes)) show_boxes <- (task == "detect" || !show_polygons)

  if (n_imgs > 1L) {
    if (is.null(mfrow)) {
      nc <- ceiling(sqrt(n_imgs))
      nr <- ceiling(n_imgs / nc)
      mfrow <- c(nr, nc)
    }
    op <- graphics::par(mfrow = mfrow, mar = c(1, 1, 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  img_keys <- names(imgs)

  for (k in seq_along(imgs)) {
    im_k <- imgs[[k]]
    k_id <- if (!is.null(img_keys) && length(img_keys) >= k) img_keys[k] else NULL

    df_b <- if (!is.null(k_id) && !is.null(x$reference_boxes_map[[k_id]])) {
      x$reference_boxes_map[[k_id]]
    } else if (!is.null(k_id) && !is.null(x$reference_boxes_map[[basename(k_id)]])) {
      x$reference_boxes_map[[basename(k_id)]]
    } else if (!is.null(k_id) && !is.null(x$reference_boxes_map[[tools::file_path_sans_ext(basename(k_id))]])) {
      x$reference_boxes_map[[tools::file_path_sans_ext(basename(k_id))]]
    } else if (k == 1L && !is.null(x$detected_boxes)) {
      x$detected_boxes
    } else {
      data.frame()
    }

    polys <- if (!is.null(k_id) && !is.null(x$reference_polys_map[[k_id]])) {
      x$reference_polys_map[[k_id]]
    } else if (!is.null(k_id) && !is.null(x$reference_polys_map[[basename(k_id)]])) {
      x$reference_polys_map[[basename(k_id)]]
    } else if (!is.null(k_id) && !is.null(x$reference_polys_map[[tools::file_path_sans_ext(basename(k_id))]])) {
      x$reference_polys_map[[tools::file_path_sans_ext(basename(k_id))]]
    } else if (k == 1L && !is.null(x$detected_polygons)) {
      x$detected_polygons
    } else {
      list()
    }

    pts_k <- if (!is.null(k_id) && !is.null(x$exemplar_points_map[[k_id]])) {
      x$exemplar_points_map[[k_id]]
    } else if (!is.null(k_id) && !is.null(x$exemplar_points_map[[basename(k_id)]])) {
      x$exemplar_points_map[[basename(k_id)]]
    } else if (!is.null(k_id) && !is.null(x$exemplar_points_map[[tools::file_path_sans_ext(basename(k_id))]])) {
      x$exemplar_points_map[[tools::file_path_sans_ext(basename(k_id))]]
    } else if (k == 1L && !is.null(x$exemplar_points)) {
      x$exemplar_points
    } else {
      NULL
    }

    plot(im_k)
    w <- dim(im_k)[1]
    h <- dim(im_k)[2]

    n_inst <- if (task == "segment" && length(polys) > 0L) {
      length(polys)
    } else {
      nrow(df_b)
    }

    if (isTRUE(rainbow) && n_inst > 0L) {
      cols_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
      cols_poly <- grDevices::adjustcolor(cols_border, alpha.f = 0.40)
    } else {
      cols_border <- rep(col_border, max(1L, n_inst))
      cols_poly <- rep(col_poly, max(1L, n_inst))
    }

    # Draw segmentation polygons if available
    if (isTRUE(show_polygons) && length(polys) > 0L) {
      for (i in seq_along(polys)) {
        poly <- polys[[i]]
        if (!is.null(poly) && length(poly) >= 6L) {
          px <- poly[c(TRUE, FALSE)] * w
          py <- poly[c(FALSE, TRUE)] * h
          cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
          cp <- cols_poly[((i - 1L) %% length(cols_poly)) + 1L]
          graphics::polygon(px, py, col = cp, border = cb, lwd = max(1, lwd - 1))
        }
      }
    }

    # Draw exemplar points
    if (isTRUE(show_points) && !is.null(pts_k) && length(pts_k) > 0L) {
      if (is.matrix(pts_k) || is.data.frame(pts_k)) {
        graphics::points(pts_k[, 1], pts_k[, 2], col = col_point, pch = 19, cex = 1.3)
      }
    }

    # Draw detected bounding boxes
    if (isTRUE(show_boxes) && nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        cb <- cols_border[((r - 1L) %% length(cols_border)) + 1L]
        graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = cb, lwd = lwd)
      }
    }

    # Draw labels and badges
    if (nrow(df_b) > 0L && (isTRUE(show_class) || isTRUE(show_conf))) {
      for (r in seq_len(nrow(df_b))) {
        cb <- cols_border[((r - 1L) %% length(cols_border)) + 1L]
        lbl_parts <- character(0)
        if (isTRUE(show_class)) {
          cname <- if (!is.null(df_b$class_name) && length(df_b$class_name) >= r && nzchar(df_b$class_name[r])) {
            df_b$class_name[r]
          } else if (!is.null(df_b$label) && length(df_b$label) >= r && nzchar(df_b$label[r])) {
            df_b$label[r]
          } else {
            x$label[1]
          }
          if (nzchar(cname)) lbl_parts <- c(lbl_parts, cname)
        }
        if (isTRUE(show_conf) && !is.null(df_b$conf) && length(df_b$conf) >= r && !is.na(df_b$conf[r])) {
          lbl_parts <- c(lbl_parts, sprintf("%.2f", df_b$conf[r]))
        }
        tag <- paste(lbl_parts, collapse = " ")

        if (nzchar(tag)) {
          tw <- graphics::strwidth(tag, cex = 0.8)
          th <- graphics::strheight(tag, cex = 0.8)
          y_top <- max(th + 6, df_b$ymin[r])
          if (isTRUE(badge)) {
            bg_tag <- grDevices::adjustcolor(cb, alpha.f = 0.90)
            graphics::rect(df_b$xmin[r], y_top - th - 6, df_b$xmin[r] + tw + 6, y_top, col = bg_tag, border = NA)
            graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = "#FFFFFF", cex = 0.8, font = 2, adj = c(0, 1))
          } else {
            graphics::text(df_b$xmin[r] + 3, y_top - 3, labels = tag, col = cb, cex = 0.8, font = 2, adj = c(0, 1))
          }
        }
      }
    }

    n_det_k <- if (task == "segment" && length(polys) > 0L) {
      sum(!vapply(polys, is.null, logical(1)))
    } else {
      nrow(df_b)
    }

    if (n_imgs > 1L) {
      fn_label <- if (!is.null(k_id) && nzchar(k_id)) basename(k_id) else sprintf("Reference %d", k)
      graphics::title(sprintf("%s: %d %s", fn_label, n_det_k, if (task == "segment") "segmented" else "detected"),
                      col.main = if (isTRUE(rainbow)) "#222222" else col_border,
                      cex.main = 1.0)
    } else {
      graphics::title(sprintf("Few-Shot Model: '%s' (%d %s) [%s]",
                              x$label, n_det_k,
                              if (task == "segment") "segmented" else "detected",
                              task),
                      col.main = col_border)
    }
  }

  invisible(x)
}

#' @export
print.yolo_fewshot_model <- function(x, ...) {
  cli::cli_h2("pliman Few-Shot Exemplar Model")
  cls_val <- if (!is.null(x$classes)) x$classes else x$label
  cls_len <- length(cls_val)
  cls_lbl <- if (cls_len > 1L) "Target Classes   " else "Target Class     "
  cls_txt <- paste(cls_val, collapse = ", ")
  cli::cli_text("{cls_lbl}: {.val {cls_txt}}")
  n_refs <- if (!is.null(x$reference_image_paths)) length(unique(x$reference_image_paths)) else if (!is.null(x$reference_images)) length(x$reference_images) else 1L
  cli::cli_text("Reference Images : {.val {n_refs}}")
  n_det <- if (!is.null(x$detected_boxes)) nrow(x$detected_boxes) else 0L
  if (!is.null(x$reference_boxes_map) && length(x$reference_boxes_map) > 0L) {
    ref_keys <- if (!is.null(x$reference_image_paths)) intersect(x$reference_image_paths, names(x$reference_boxes_map)) else names(x$reference_boxes_map)
    if (length(ref_keys) > 0L) {
      n_det <- sum(vapply(ref_keys, function(k) nrow(x$reference_boxes_map[[k]]), integer(1)))
    }
  }
  cli::cli_text("Detected on Refs : {.val {n_det}} instances")
  cli::cli_text("Confidence Thresh: {.val {x$conf_threshold}}")
  cli::cli_text("Task             : {.val {x$task}}")
  cli::cli_text("Engine           : {.val {x$engine}}")
  invisible(x)
}
