# =============================================================================
# benchmark/panelio/to_panelio.R
# -----------------------------------------------------------------------------
# Converts the A1..A12 benchmark datasets (stored as .rds) into panelio native
# archives. One archive = one balanced panel, so each of the 30 replicates
# becomes its own archive:
#     panelio/A1/A1_rep01_v1.tar.gz ... A1_rep30_v1.tar.gz
#     ...
#     panelio/A12/A12_rep01_v1.tar.gz ... A12_rep30_v1.tar.gz
#
# The abundance value archived is relab_true: the relative-abundance panel,
# which is the panelio design contract (every observation x time layer sums to
# 1). Raw counts Y, depth and effect truth cannot be represented by panelio's
# schema and are retained in the source .rds files.
#
# Usage: Rscript to_panelio.R [rds_dir] [out_dir]
# =============================================================================

rds_dir <- if (length(args <- commandArgs(trailingOnly = TRUE)) >= 1) args[1] else "."
out_dir <- if (length(args) >= 2) args[2] else file.path(rds_dir, "panelio")

pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE)

# Build a canonical panelio long table from one replicate's relab_true.
build_long_from_relab <- \(sim, r) {
  N <- dim(sim$relab_true)[1L]; D <- dim(sim$relab_true)[2L]; T <- dim(sim$relab_true)[3L]
  obs  <- sprintf("subject_%03d", seq_len(N))
  feat <- sprintf("taxa_%03d", seq_len(D))
  tm   <- sprintf("time_%03d", seq_len(T))
  long <- expand.grid(feature = feat, observation = obs, time_identifier = tm,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  oidx <- match(long$observation, obs)
  fidx <- match(long$feature, feat)
  tidx <- match(long$time_identifier, tm)
  long$abundance_value <- sim$relab_true[, , , r][cbind(oidx, fidx, tidx)]
  long$elapsed_time    <- sim$time[tidx]
  long$label           <- factor(ifelse(sim$treatment[oidx] == 1L, "treat", "ctrl"),
                                 levels = c("ctrl", "treat"))
  long$is_intervention <- tidx >= as.integer(sim$t0)
  long[, c("observation", "time_identifier", "elapsed_time", "feature",
           "abundance_value", "label", "is_intervention"), drop = FALSE]
}

cat("== convert A1..A12 (.rds) -> panelio archives ==\n")
for (i in 1:12) {
  nm <- sprintf("A%d", i)
  sim <- readRDS(file.path(rds_dir, sprintf("%s.rds", nm)))
  ds_dir <- file.path(out_dir, nm)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  truth_comment <- if (length(sim$effect_taxa) > 0L) {
    sprintf("effect_taxa=%s; effect_signs=%s; t0=%d",
            paste(sim$effect_taxa, collapse = ","),
            paste(sim$effect_signs, collapse = ","),
            sim$t0)
  } else {
    sprintf("null scenario; effect_size=0; t0=%d", sim$t0)
  }
  for (r in seq_len(30L)) {
    long_df <- build_long_from_relab(sim, r)
    archive <- panelio::write_panel_archive(
      long_df,
      dataset_name = sprintf("%s_rep%02d", nm, r),
      version = "v1",
      out_dir = ds_dir,
      manifest_comments = c(
        sprintf("%s replicate %02d of %d (relab_true panel)", nm, r, dim(sim$Y)[4L]),
        truth_comment
      )
    )
    if (r %% 10L == 0L) cat(sprintf("  %s rep %02d done -> %s\n", nm, r, basename(archive)))
  }
  cat(sprintf("  %-3s: %d archives written\n", nm, 30L))
}
cat("== done ==\n")
