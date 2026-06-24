# ============================================================
# ACCC Petrol Price Chart Extractor  (hardened)
# ============================================================
library(png)
library(tesseract)
library(tidyr)
library(dplyr)
library(lubridate)

Sys.setenv(TZ = "Australia/Sydney")
today <- Sys.Date()

# Domain knowledge: ULP gridlines are 5-cent steps; real metro
# prices sit well inside this band and are always positive.
VALUE_STEP <- 5
PLAUSIBLE  <- c(120, 260)


# ---- FUNCTION: calibrate_axis --------------------------------
# Takes the gridline rows + raw (possibly garbage) OCR reads and
# returns a clean, evenly-spaced axis. The geometry of the
# gridlines fixes the SPACING; OCR is only used to pin the
# absolute LEVEL, and only plausible reads are allowed to vote.
calibrate_axis <- function(axis_map, value_step = VALUE_STEP,
                           plausible = PLAUSIBLE) {

  axis_map <- axis_map[order(axis_map$image_row), ]   # top -> bottom
  rows <- axis_map$image_row

  if (nrow(axis_map) < 2)
    stop("Fewer than 2 gridlines detected — cannot calibrate axis.")

  # Integer slot index for each gridline from geometry alone.
  # round(gap / median gap) tolerates an occasional missed line.
  gaps <- diff(rows)
  slot <- c(0, cumsum(round(gaps / median(gaps))))

  # Filter raw OCR to physically possible label values BEFORE
  # any of it is trusted. value decreases as image_row increases,
  # so the top gridline (slot 0) holds the largest value.
  ocr <- axis_map$axis_value
  ok  <- !is.na(ocr) &
         ocr >= plausible[1] & ocr <= plausible[2] &
         (ocr %% value_step == 0)

  if (!any(ok))
    stop("No plausible y-axis labels survived filtering — aborting ",
         "instead of fabricating an axis.")

  # Each good read implies a value for the TOP gridline. Vote.
  v_top_votes <- ocr[ok] + value_step * slot[ok]
  tbl   <- sort(table(v_top_votes), decreasing = TRUE)
  v_top <- as.numeric(names(tbl)[1])
  agree <- as.integer(tbl[1])

  cat(sprintf("  Axis calibration: top gridline = %g  (%d/%d reads agree)\n",
              v_top, agree, sum(ok)))
  if (agree < 2)
    cat("  WARNING: weak axis consensus (only one supporting read) — ",
        "verify the check plot.\n")

  # Rebuild a perfectly even axis from the agreed level.
  axis_map$axis_value <- v_top - value_step * slot
  axis_map
}


# ---- FUNCTION: extract_chart_data ----------------------------

extract_chart_data <- function(image_path, city_name, target_dots = 45,
                               start_radius = 3.0, step = 0.02,
                               max_iter = 200, dot_tol = 1,
                               value_step = VALUE_STEP,
                               plausible = PLAUSIBLE) {

  cat("\n========================================\n")
  cat("Processing:", city_name, "\n")
  cat("========================================\n\n")

  img <- readPNG(image_path)
  cat("Image:", nrow(img), "x", ncol(img), "x", dim(img)[3], "\n")

  r_chan <- img[,,1]
  g_chan <- img[,,2]
  b_chan <- img[,,3]

  # -- Purple mask --
  purple_mask <- (r_chan > 0.15 & r_chan < 0.65) &
    (g_chan < 0.25) &
    (b_chan > 0.15 & b_chan < 0.65)

  # -- Erosion function --
  erode_mask <- function(mask, radius) {
    nr <- nrow(mask); nc <- ncol(mask)
    full_r <- floor(radius); frac_r <- radius - full_r
    eroded <- mask
    for (pass in seq_len(full_r)) {
      new_mask <- eroded
      for (dr in -1:1) for (dc in -1:1) {
        r_from <- max(1,1+dr):min(nr,nr+dr); r_to <- max(1,1-dr):min(nr,nr-dr)
        c_from <- max(1,1+dc):min(nc,nc+dc); c_to <- max(1,1-dc):min(nc,nc-dc)
        new_mask[r_to, c_to] <- new_mask[r_to, c_to] & eroded[r_from, c_from]
      }
      eroded <- new_mask
    }
    if (frac_r > 0) {
      cardinal <- list(c(-1,0), c(1,0), c(0,-1), c(0,1))
      new_mask <- eroded
      for (nb in cardinal) {
        dr <- nb[1]; dc <- nb[2]
        r_from <- max(1,1+dr):min(nr,nr+dr); r_to <- max(1,1-dr):min(nr,nr-dr)
        c_from <- max(1,1+dc):min(nc,nc+dc); c_to <- max(1,1-dc):min(nc,nc-dc)
        new_mask[r_to, c_to] <- new_mask[r_to, c_to] & eroded[r_from, c_from]
      }
      eroded <- new_mask
    }
    return(eroded)
  }

  # -- Blob labelling --
  label_blobs <- function(mask) {
    nr <- nrow(mask); nc <- ncol(mask)
    labels <- matrix(0L, nr, nc); current_label <- 0L
    flood_fill <- function(r0, c0, lab) {
      queue <- list(c(r0, c0)); head <- 1
      while (head <= length(queue)) {
        rc <- queue[[head]]; head <- head + 1
        r <- rc[1]; c <- rc[2]
        if (r<1||r>nr||c<1||c>nc) next
        if (!mask[r,c]||labels[r,c]!=0) next
        labels[r,c] <<- lab
        for (dr in -1:1) for (dc in -1:1) {
          if (dr==0&&dc==0) next
          queue[[length(queue)+1]] <- c(r+dr, c+dc)
        }
      }
    }
    for (r in 1:nr) for (c in 1:nc) {
      if (mask[r,c] && labels[r,c]==0) {
        current_label <- current_label + 1L
        flood_fill(r, c, current_label)
      }
    }
    return(labels)
  }

  count_dots <- function(radius, min_blob_px = 5) {
    eroded <- erode_mask(purple_mask, radius)
    bl <- label_blobs(eroded)
    n <- max(bl); if (n==0) return(0)
    sum(tabulate(bl) >= min_blob_px)
  }

  # -- Auto-tune erosion radius (softened) ---------------------
  # Track the radius that gets CLOSEST to target; accept within
  # dot_tol; break early if the search starts oscillating.
  radius <- start_radius
  best_radius <- radius; best_diff <- Inf; best_n <- NA_integer_
  prev_dir <- 0L; converged <- FALSE
  cat("Auto-tuning for", target_dots, "dots (tolerance +/-", dot_tol, ")...\n")

  for (iter in 1:max_iter) {
    n_found <- count_dots(radius)
    d <- abs(n_found - target_dots)
    if (d < best_diff) { best_diff <- d; best_radius <- radius; best_n <- n_found }

    if (n_found == target_dots) {
      cat(sprintf("  Converged exactly: iter %d, radius %.2f\n", iter, radius))
      converged <- TRUE; break
    }

    dir <- if (n_found > target_dots) 1L else -1L  # +radius -> fewer dots
    # Oscillating around the target and already close enough -> stop.
    if (prev_dir != 0L && dir != prev_dir && best_diff <= dot_tol) {
      cat(sprintf("  Oscillating near target; accepting best: radius %.2f (%d dots)\n",
                  best_radius, best_n))
      converged <- TRUE; break
    }

    radius <- radius + dir * step
    prev_dir <- dir
    if (radius < 0.1) { radius <- 0.1; break }   # guard against runaway
  }

  radius <- best_radius
  if (!converged) {
    if (best_diff <= dot_tol) {
      cat(sprintf("  Accepted within tolerance: radius %.2f (%d dots)\n",
                  best_radius, best_n))
    } else {
      cat(sprintf("  WARNING: best was %d dots (target %d) at radius %.2f\n",
                  best_n, target_dots, best_radius))
    }
  }

  # -- Extract dot positions --
  eroded <- erode_mask(purple_mask, radius)
  bl <- label_blobs(eroded); n_blobs <- max(bl)
  pts <- data.frame(dot_id=integer(), median_col=numeric(),
                    median_row=numeric(), n_pixels=integer())
  for (i in 1:n_blobs) {
    idx <- which(bl==i, arr.ind=TRUE)
    if (nrow(idx) >= 5)
      pts <- rbind(pts, data.frame(dot_id=nrow(pts)+1, median_col=median(idx[,2]),
                                   median_row=median(idx[,1]), n_pixels=nrow(idx)))
  }
  pts <- pts[order(pts$median_col), ]; pts$dot_id <- 1:nrow(pts)

  # Refine positions from original mask
  for (i in 1:nrow(pts)) {
    cr <- round(pts$median_row[i]); cc <- round(pts$median_col[i])
    rr <- max(1,cr-10):min(nrow(img),cr+10)
    cr2 <- max(1,cc-10):min(ncol(img),cc+10)
    lm <- purple_mask[rr, cr2]; lc <- which(lm, arr.ind=TRUE)
    if (nrow(lc) > 0) {
      pts$median_row[i] <- median(rr[lc[,1]])
      pts$median_col[i] <- median(cr2[lc[,2]])
    }
  }
  cat("Dots detected:", nrow(pts), "\n")

  # -- Gridlines --
  grey_mask <- (abs(r_chan-g_chan)<0.05)&(abs(g_chan-b_chan)<0.05)&(r_chan>0.55)&(r_chan<0.95)
  row_gc <- rowSums(grey_mask)
  gl_rows <- which(row_gc > max(row_gc)*0.4)
  rd <- diff(gl_rows); cid <- cumsum(c(1, rd>3))
  grid_centres <- sort(tapply(gl_rows, cid, median))

  dark_mask <- (r_chan<0.3)&(g_chan<0.3)&(b_chan<0.3)
  rdc <- rowSums(dark_mask)
  dlr <- which(rdc > ncol(img)*0.3)
  if (length(dlr) > 0) {
    dld <- diff(dlr); dlc <- cumsum(c(1, dld>3))
    dlcen <- sort(tapply(dlr, dlc, median))
    bar <- max(dlcen)
    if (all(abs(grid_centres-bar)>5)) grid_centres <- sort(c(grid_centres, bar))
  }
  cat("Gridlines:", length(grid_centres), "\n")

  # -- Adaptive y-axis label crop -------------------------------
  # Find the y-axis border (tallest dark vertical line in the left
  # region); the price digits are the widest dark column-cluster
  # sitting just LEFT of it. Everything is re-derived each run, so
  # a few-pixel shift in axis position is handled automatically.
  left_region <- 1:max(2, round(ncol(img) * 0.30))
  col_dark    <- colSums(dark_mask[, left_region, drop = FALSE])
  axis_col    <- left_region[which.max(col_dark)]

  # Collect dark columns in a tight band around each gridline.
  label_cols <- c()
  for (gc in grid_centres) {
    rr <- max(1, round(gc)-15):min(nrow(img), round(gc)+15)
    cc <- 1:min(ncol(img), axis_col)
    ld <- which(dark_mask[rr, cc], arr.ind = TRUE)
    if (nrow(ld) > 0) label_cols <- c(label_cols, ld[,2])
  }

  # Default window in case clustering is degenerate.
  sl <- max(1, round(quantile(label_cols, 0.05, names = FALSE)))
  sr <- min(axis_col - 2, round(quantile(label_cols, 0.95, names = FALSE)))

  dig <- sort(unique(label_cols[label_cols < axis_col - 2]))
  if (length(dig) > 2) {
    grp   <- cumsum(c(1, diff(dig) > 10))          # split clusters on >10px gaps
    spans <- tapply(dig, grp, function(x) max(x) - min(x))
    best  <- as.integer(names(which.max(spans)))   # widest cluster = the number
    cluster <- dig[grp == best]
    sl <- max(1, min(cluster) - 4)
    sr <- min(axis_col - 2, max(cluster) + 4)
  }
  cat(sprintf("Label crop: cols %d..%d (axis border at %d)\n", sl, sr, axis_col))

  # -- OCR ensemble (multiple page-seg modes) -------------------
  ocr_engines <- lapply(c(7, 8, 10), function(psm)
    tesseract(options = list(tessedit_char_whitelist = "0123456789",
                             tessedit_pageseg_mode = psm)))

  axis_map <- data.frame(gridline = integer(), image_row = numeric(),
                         axis_value = numeric(), stringsAsFactors = FALSE)

  for (i in seq_along(grid_centres)) {
    rc <- round(grid_centres[i])
    rt <- max(1, rc-22); rb <- min(nrow(img), rc+22)
    strip <- img[rt:rb, sl:sr, , drop = FALSE]
    sg  <- 0.299*strip[,,1] + 0.587*strip[,,2] + 0.114*strip[,,3]
    sbw <- ifelse(sg < 0.4, 0, 1)
    so  <- array(sbw, dim = c(dim(sbw), 3))
    nr  <- nrow(so); nc <- ncol(so)
    sb  <- so[rep(1:nr, each = 4), rep(1:nc, each = 4), , drop = FALSE]
    pd  <- array(1, dim = c(nrow(sb)+40, ncol(sb)+40, 3))
    pd[21:(20+nrow(sb)), 21:(20+ncol(sb)), ] <- sb
    tf  <- tempfile(fileext = ".png"); writePNG(pd, tf)

    reads <- c()
    for (eng in ocr_engines) {
      ct <- gsub("[^0-9]", "", ocr(tf, engine = eng))
      nv <- suppressWarnings(as.numeric(ct))
      if (!is.na(nv)) reads <- c(reads, nv)
    }
    file.remove(tf)

    # Prefer ensemble reads that are physically possible; if none,
    # keep the most common raw read so calibrate_axis can still
    # judge it (and likely discard it).
    good   <- reads[reads >= plausible[1] & reads <= plausible[2] &
                    reads %% value_step == 0]
    pool   <- if (length(good) > 0) good else reads
    chosen <- if (length(pool) > 0)
                as.numeric(names(sort(table(pool), decreasing = TRUE))[1]) else NA

    cat(sprintf("  Grid %2d | row %4d | reads: %-18s | chosen: %s\n",
                i, rc,
                if (length(reads)) paste(reads, collapse = ",") else "-",
                ifelse(is.na(chosen), "NONE", chosen)))

    axis_map <- rbind(axis_map, data.frame(gridline = i,
                                           image_row = grid_centres[i],
                                           axis_value = chosen,
                                           stringsAsFactors = FALSE))
  }

  # -- Calibrate the axis (consensus vote) ----------------------
  axis_map <- calibrate_axis(axis_map, value_step = value_step,
                             plausible = plausible)

  ac <- axis_map[!is.na(axis_map$axis_value), ]
  tr <- ac$image_row[1];          tv <- ac$axis_value[1]
  br <- ac$image_row[nrow(ac)];   bv <- ac$axis_value[nrow(ac)]

  # -- Scale to real values --
  ppu <- (br-tr)/(bv-tv)
  pts$value <- round(tv + (pts$median_row-tr)/ppu, 1)
  cat("Values:", min(pts$value), "to", max(pts$value), "\n")

  # -- Output sanity check (last line of defence) --------------
  if (any(pts$value < plausible[1] | pts$value > plausible[2])) {
    n_bad <- sum(pts$value < plausible[1] | pts$value > plausible[2])
    cat(sprintf("  WARNING: %d extracted price(s) outside %g-%g c/L — ",
                n_bad, plausible[1], plausible[2]))
    cat("flagging", city_name, "for manual review.\n")
  }

  # -- Save overlay check plot --
  if (!dir.exists("data")) dir.create("data", recursive = TRUE)
  plot_file <- paste0("data/check_", gsub(" ", "_", tolower(city_name)), "_", today, ".png")
  png(plot_file, width = 1200, height = 800)

  plot(1, type = "n",
       xlim = c(0, ncol(img)), ylim = c(0, nrow(img)), asp = 1,
       xlab = "x (pixels)", ylab = "y (pixels)",
       main = paste0(city_name, " — ", today, " — gridlines (red) + dots (green)"))
  rasterImage(img, 0, 0, ncol(img), nrow(img))

  for (j in 1:nrow(ac)) {
    y_plot <- nrow(img) - ac$image_row[j]
    abline(h = y_plot, col = "red", lty = 2, lwd = 1.5)
    text(ncol(img) - 50, y_plot + 10,
         labels = ac$axis_value[j], col = "red", cex = 1, font = 2)
  }
  points(pts$median_col, nrow(img) - pts$median_row,
         col = "green", pch = 1, cex = 2, lwd = 2)
  text(pts$median_col, nrow(img) - pts$median_row + 15,
       labels = pts$dot_id, col = "green", cex = 0.6)

  dev.off()
  cat("Saved check plot:", plot_file, "\n")

  data.frame(point=pts$dot_id, city=city_name, value=pts$value)
}


# ============================================================
# MAIN: Download images, extract data, append to CSV
# ============================================================

base_url <- "https://www.accc.gov.au/sites/www.accc.gov.au/files/fuelwatch/"

cities <- data.frame(
  file = c("sydney-metro-ulp.png", "melbourne-metro-ulp.png",
           "brisbane-metro-ulp.png", "adelaide-metro-ulp.png",
           "perth-metro-ulp.png"),
  city = c("Sydney","Melbourne","Brisbane","Adelaide","Perth"),
  stringsAsFactors = FALSE
)

all_results <- data.frame()

for (i in 1:nrow(cities)) {
  img_url  <- paste0(base_url, cities$file[i])
  tmp_path <- tempfile(fileext = ".png")
  Sys.sleep(10)
  tryCatch({
    download.file(img_url, destfile = tmp_path, mode = "wb", quiet = TRUE)
    cat("Downloaded:", cities$city[i], "\n")
    result <- extract_chart_data(image_path = tmp_path, city_name = cities$city[i])
    all_results <- rbind(all_results, result)
  }, error = function(e) {
    cat("ERROR processing", cities$city[i], ":", conditionMessage(e), "\n")
  })

  if (file.exists(tmp_path)) file.remove(tmp_path)
}

# -- Pivot wide and add dates --
all_wide <- all_results %>%
  pivot_wider(names_from = city, values_from = value) %>%
  mutate(forecast_date = today,
         date = today -1 - days(max(point) - point)) %>%
  select(forecast_date, date, point, everything())

cat("\nToday's extraction:\n")
print(all_wide)

# -- Append to historical CSV --
output_file <- "data/petrol_prices.csv"

# Ensure the data directory exists
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

if (file.exists(output_file)) {
  # Read existing data, append new rows, remove duplicates
  existing <- read.csv(output_file, stringsAsFactors = FALSE)
  existing$date <- as.Date(existing$date)
  existing$forecast_date <- as.Date(existing$forecast_date)

  combined <- bind_rows(existing, all_wide) %>%
    group_by(forecast_date, date) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    arrange(date)

  write.csv(combined, output_file, row.names = FALSE)
  cat("\nAppended to", output_file, "- total rows:", nrow(combined), "\n")
} else {
  write.csv(all_wide, output_file, row.names = FALSE)
  cat("\nCreated", output_file, "with", nrow(all_wide), "rows\n")
}

cat("\nDone!\n")
