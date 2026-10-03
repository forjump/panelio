# 14_summarize.R
# Consolidate perfeature CSVs into metric tables for the paper.
# Usage: Rscript 14_summarize.R sim|abl|sen_filter|sen_normalize|sen_zero|sen_band|sen_bootstrap
# Reads results/sim_s{1..4}_perfeature.csv, results/abl_*.csv, results/sen_*.csv.
# Metrics per (dataset x method | dataset x variant | dataset x config):
#   sign_hit, power_fdr, power_placebo, signed_power_fdr, typeI_fdr, typeI_placebo,
#   bias, mae, recon_eff (|est-truth| on effect taxa), recon_null (mean|est| on null taxa),
#   pre_rmse, pre_snr, weight_sum (scm only), n_kept
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1L]
res_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results"
truth <- read.csv(file.path(res_dir, "../truth/truth_table.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)

read_all <- function(pattern) {
  files <- list.files(res_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) stop("no files match ", pattern)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
}

truth_vec <- function(ds) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  if (nchar(tr$effect_taxa) == 0L) {
    data.frame(feature = character(0L), truth = numeric(0L), sign = integer(0L),
               stringsAsFactors = FALSE)
  } else {
    idx <- as.integer(strsplit(tr$effect_taxa, ";")[[1L]])
    sgn <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])
    data.frame(feature = sprintf("taxa_%03d", idx),
               truth = sgn * as.numeric(tr$effect_size), sign = sgn,
               stringsAsFactors = FALSE)
  }
}

summarize_frame <- function(pf, group_col) {
  ds_list <- unique(pf$dataset)
  out <- lapply(ds_list, function(ds) {
    tv <- truth_vec(ds)
    eff <- tv$feature
    grp <- unique(pf[[group_col]])
    lapply(grp, function(g) {
      sub <- pf[pf$dataset == ds & pf[[group_col]] == g, , drop = FALSE]
      is_scm <- unique(sub$method)[1L] == "scm" || "variant" %in% names(sub)
      est <- sub$estimate
      fdr <- sub$fdr
      pfd <- if ("placebo_fdr" %in% names(sub)) sub$placebo_fdr else NA_real_
      if (length(eff) > 0L) {
        eff_mask <- sub$feature %in% eff
        tv_map <- tv$truth[match(sub$feature[eff_mask], tv$feature)]
        tv_sign <- tv$sign[match(sub$feature[eff_mask], tv$feature)]
        null_mask <- !eff_mask
        sign_hit <- mean(sign(est[eff_mask]) == tv_sign, na.rm = TRUE)
        power_fdr <- mean(fdr[eff_mask] < 0.05, na.rm = TRUE)
        signed_power_fdr <- mean(fdr[eff_mask] < 0.05 & sign(est[eff_mask]) == tv_sign, na.rm = TRUE)
        power_placebo <- if (any(!is.na(pfd))) mean(pfd[eff_mask] < 0.05, na.rm = TRUE) else NA_real_
        bias <- mean(est[eff_mask] - tv_map, na.rm = TRUE)
        mae <- mean(abs(est[eff_mask] - tv_map), na.rm = TRUE)
        recon_null <- mean(abs(est[null_mask]), na.rm = TRUE)
        typeI_fdr <- mean(fdr[null_mask] < 0.05, na.rm = TRUE)
        typeI_placebo <- if (any(!is.na(pfd))) mean(pfd[null_mask] < 0.05, na.rm = TRUE) else NA_real_
      } else {
        sign_hit <- power_fdr <- signed_power_fdr <- NA_real_
        bias <- mae <- recon_null <- NA_real_
        power_placebo <- NA_real_
        typeI_fdr <- mean(fdr < 0.05, na.rm = TRUE)
        typeI_placebo <- if (any(!is.na(pfd))) mean(pfd < 0.05, na.rm = TRUE) else NA_real_
      }
      data.frame(
        dataset = ds, group = g,
        n_est = nrow(sub), n_kept = if ("n_kept" %in% names(sub)) mean(sub$n_kept, na.rm = TRUE) else NA_real_,
        sign_hit = sign_hit, power_fdr = power_fdr, signed_power_fdr = signed_power_fdr,
        power_placebo = power_placebo, typeI_fdr = typeI_fdr, typeI_placebo = typeI_placebo,
        bias = bias, mae = mae, recon_null = recon_null,
        pre_rmse = if ("pre_rmse" %in% names(sub)) mean(sub$pre_rmse, na.rm = TRUE) else NA_real_,
        pre_snr = if ("pre_snr" %in% names(sub)) mean(sub$pre_snr, na.rm = TRUE) else NA_real_,
        weight_sum = if ("weight_sum" %in% names(sub)) mean(sub$weight_sum, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    })
  })
  do.call(rbind, unlist(out, recursive = FALSE))
}

if (mode == "sim") {
  pf <- read_all("sim_s._perfeature\\.csv")
  cat("sim rows:", nrow(pf), " expected 12*30*4*50 = 72000\n")
  # reconciliation
  n_check <- table(pf$dataset, pf$method)
  cat("dataset x method counts (rows/1500):\n"); print(n_check / 1500)
  sm <- summarize_frame(pf, "method")
  sm$method <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_simulation.csv"), row.names = FALSE)
  cat("WROTE summary_simulation.csv\n")
  print(sm)
} else if (mode == "abl") {
  pf <- read_all("abl_.*_perfeature\\.csv")
  cat("abl rows:", nrow(pf), "\n")
  cat("variants:", paste(unique(pf$variant), collapse = ", "), "\n")
  sm <- summarize_frame(pf, "variant")
  sm$variant <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_ablation.csv"), row.names = FALSE)
  cat("WROTE summary_ablation.csv\n")
  print(sm)
} else if (mode == "sen_filter") {
  pf <- read_all("sen_.*_filter_perfeature\\.csv")
  sm <- summarize_frame(pf, "config")
  sm$config <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_sensitivity_filter.csv"), row.names = FALSE)
  cat("WROTE summary_sensitivity_filter.csv\n"); print(sm)
} else if (mode == "sen_normalize") {
  pf <- read_all("sen_.*_normalize_perfeature\\.csv")
  sm <- summarize_frame(pf, "config")
  sm$config <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_sensitivity_normalize.csv"), row.names = FALSE)
  cat("WROTE summary_sensitivity_normalize.csv\n"); print(sm)
} else if (mode == "sen_zero") {
  pf <- read_all("sen_.*_zero_perfeature\\.csv")
  sm <- summarize_frame(pf, "config")
  sm$config <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_sensitivity_zero.csv"), row.names = FALSE)
  cat("WROTE summary_sensitivity_zero.csv\n"); print(sm)
} else if (mode == "sen_band") {
  pf <- read_all("sen_.*_band_perfeature\\.csv")
  sm <- summarize_frame(pf, "config")
  sm$config <- sm$group
  sm$group <- NULL
  write.csv(sm, file.path(res_dir, "summary_sensitivity_band.csv"), row.names = FALSE)
  cat("WROTE summary_sensitivity_band.csv\n"); print(sm)
} else if (mode == "sen_bootstrap") {
  pf <- read_all("sen_.*_bootstrap_perfeature\\.csv")
  tv <- lapply(unique(pf$dataset), truth_vec)
  names(tv) <- unique(pf$dataset)
  boot_sum <- lapply(split(seq_len(nrow(pf)), interaction(pf$dataset, pf$config)), function(ix) {
    sub <- pf[ix, , drop = FALSE]
    ds <- sub$dataset[1L]
    cfg <- sub$config[1L]
    t <- tv[[ds]]
    eff <- t$feature
    eff_mask <- sub$feature %in% eff
    width <- mean(sub$bootstrap_upper[eff_mask] - sub$bootstrap_lower[eff_mask], na.rm = TRUE)
    se_boot <- mean((sub$bootstrap_upper[eff_mask] - sub$bootstrap_lower[eff_mask]) /
                    (2 * qnorm(0.975)), na.rm = TRUE)
    cover <- if (nrow(t) > 0L) {
      cov_flag <- sub$bootstrap_lower[eff_mask] <= t$truth[match(sub$feature[eff_mask], t$feature)] &
                  sub$bootstrap_upper[eff_mask] >= t$truth[match(sub$feature[eff_mask], t$feature)]
      mean(cov_flag, na.rm = TRUE)
    } else NA_real_
    data.frame(dataset = ds, config = cfg, boot_width = width, boot_se = se_boot,
               coverage = cover, stringsAsFactors = FALSE)
  })
  bm <- do.call(rbind, boot_sum)
  write.csv(bm, file.path(res_dir, "summary_sensitivity_bootstrap.csv"), row.names = FALSE)
  cat("WROTE summary_sensitivity_bootstrap.csv\n"); print(bm)
} else {
  stop("unknown mode")
}
