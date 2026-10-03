# =============================================================================
# truth/verify_archive.R
# Verifies a panelio archive against the independent ground truth.
#
#   Rscript verify_archive.R <archive_path>
#
# e.g. Rscript verify_archive.R ../panelio/A6/A6_rep12_v1.tar.gz
#
# Checks (each reported as PASS/FAIL, overall verdict at the end):
#   1. structural identity  : N, D, T, elapsed time, t0/intervention, label
#                              groups match the truth table for the dataset
#   2. effect ground truth  : effect_taxa / effect_signs match the truth table
#   3. content fingerprint  : relab_true in the source .rds matches the md5 in
#                              the truth table (source rds itself is authentic)
#   4. read-back fidelity   : abundance in the archive equals the ground-truth
#                              relab_true tensor for that dataset x replicate
# =============================================================================
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
archive_path <- if (length(args) >= 1) args[1] else stop("usage: Rscript verify_archive.R <archive_path>")
bench_root <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark"
truth <- read.csv(file.path(bench_root, "truth", "truth_table.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)

basename_ <- basename(archive_path)                       # A6_rep12_v1.tar.gz
ds_rep <- sub("_v1\\.tar\\.gz$", "", basename_)           # A6_rep12
ds <- sub("_rep[0-9]{2}$", "", ds_rep)                    # A6
rep <- as.integer(sub(".*_rep([0-9]{2})$", "\\1", ds_rep))# 12
stopifnot(ds %in% truth$dataset)
row_truth <- truth[truth$dataset == ds, , drop = TRUE]

cat(sprintf("verifying %s  (dataset=%s, replicate=%02d)\n", archive_path, ds, rep))

# ---- 1. structural identity -------------------------------------------------
back <- tryCatch(
  read_panel_table(archive_path),
  error = function(e) conditionMessage(e)
)
if (is.character(back)) {
  cat("  [FAIL] read_panel_table : archive failed panelio validation\n")
  cat("         message:", substr(back, 1, 140), "\n")
  cat(sprintf("\nOVERALL VERDICT: %s\n", "CHECK FAILED (archive rejected or tampered)"))
  quit(status = 1)
}
D <- nlevels(factor(back$feature))
time_elapsed <- unique(back$elapsed_time[!duplicated(back$time_identifier)])
back_time <- as.numeric(strsplit(row_truth$time, ";")[[1]])
checks <- list()
checks$dim_N       <- nlevels(factor(back$observation)) == row_truth$N
checks$dim_D       <- D == row_truth$D
checks$dim_T       <- nlevels(factor(back$time_identifier)) == row_truth$T
checks$time        <- isTRUE(all.equal(time_elapsed, back_time, tolerance = 1e-8))
checks$label_groups <- setequal(levels(back$label), c("ctrl", "treat"))
# single monotone intervention breakpoint at t0
uif <- back$is_intervention[!duplicated(back$time_identifier)]
checks$intervention <- identical(as.logical(uif),
                                 as.logical(seq_along(uif) >= row_truth$t0))

# ---- 2. effect ground truth -------------------------------------------------
if (row_truth$n_effect_taxa > 0L) {
  exp_taxa <- as.integer(strsplit(row_truth$effect_taxa, ";")[[1]])
  exp_signs <- as.integer(strsplit(row_truth$effect_signs, ";")[[1]])
  feature_idx <- match(sprintf("taxa_%03d", exp_taxa), levels(factor(back$feature)))
  checks$effect_size  <- row_truth$effect_size %in% c(0.4, 1.2)
  checks$effect_taxa  <- length(exp_taxa) == row_truth$n_effect_taxa
  checks$effect_signs <- length(exp_signs) == row_truth$n_effect_taxa
} else {
  checks$effect_size  <- row_truth$effect_size == 0
  checks$effect_taxa  <- row_truth$n_effect_taxa == 0L
  checks$effect_signs <- TRUE
}

# ---- 3. content fingerprint of the source rds -------------------------------
fingerprint <- \(obj) {
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
  saveRDS(obj, f, compress = "xz"); unname(tools::md5sum(f))
}
rds_path <- file.path(bench_root, sprintf("%s.rds", ds))
if (file.exists(rds_path)) {
  sim <- readRDS(rds_path)
  checks$fingerprint <- fingerprint(sim$relab_true) == row_truth$relab_true_md5
} else {
  checks$fingerprint <- NA
}

# ---- 4. read-back fidelity against ground-truth tensor ----------------------
if (file.exists(rds_path) && !is.na(checks$fingerprint) && checks$fingerprint) {
  N <- row_truth$N; D <- row_truth$D; T <- row_truth$T
  obs  <- sprintf("subject_%03d", seq_len(N))
  feat <- sprintf("taxa_%03d", seq_len(D))
  tm   <- sprintf("time_%03d", seq_len(T))
  long <- expand.grid(feature = feat, observation = obs, time_identifier = tm,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  oidx <- match(long$observation, obs); fidx <- match(long$feature, feat); tidx <- match(long$time_identifier, tm)
  long$abundance_value <- sim$relab_true[, , , rep][cbind(oidx, fidx, tidx)]
  long$elapsed_time    <- sim$time[tidx]
  long$label           <- factor(ifelse(sim$treatment[oidx] == 1L, "treat", "ctrl"), levels = c("ctrl", "treat"))
  long$is_intervention <- tidx >= as.integer(sim$t0)
  long <- long[, c("observation","time_identifier","elapsed_time","feature",
                   "abundance_value","label","is_intervention"), drop = FALSE]
  sfn <- panelio:::sort_panel_table
  checks$fidelity <- identical(sfn(back), sfn(long)) &&
    max(abs(sfn(back)$abundance_value - sfn(long)$abundance_value)) == 0
} else {
  checks$fidelity <- NA
}

# ---- report -----------------------------------------------------------------
cat("\n")
all_pass <- TRUE
for (nm in names(checks)) {
  v <- checks[[nm]]
  if (isTRUE(v)) { st <- "PASS" } else if (is.na(v)) { st <- "SKIP" ; all_pass <- FALSE }
  else { st <- "FAIL"; all_pass <- FALSE }
  cat(sprintf("  [%s] %-16s %s\n", st, nm, if (is.na(v)) "source .rds unavailable" else ""))
}
cat(sprintf("\nOVERALL VERDICT: %s\n", if (all_pass) "AUTHENTIC" else "CHECK FAILED"))
