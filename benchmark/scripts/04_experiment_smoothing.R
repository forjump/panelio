# 04_experiment_smoothing.R
# Choose the paper-default smoothing depth empirically on A1 reps 01-05 (scm).
# Variants (all CLR-preprocessed unless noted):
#   P0 raw CLR, no smoothing
#   P1 ekf gaussian qratio=0.1     P2 rts gaussian qratio=0.1
#   P3 rts gaussian qratio=1.0     P4 loess span=0.75 (brk=t0)
#   P5 savgol hw=2 (brk=t0)        P6 gaussian kernel sigma=2 (brk=t0)
#   P7 raw abundance + rts NB hierarchical (paper-depth NB, expected rel-abundance scale)
# Metric per variant (over 5 reps): bias, MAE on effect taxa, sign hit, mean |est| on null taxa.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
tr <- truth[truth$dataset == "A1", , drop = FALSE]
t0 <- as.integer(tr$t0)
eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
signs <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])
truth_est <- signs * as.numeric(tr$effect_size)
reps <- sprintf("%02d", 1:5)

run_variant <- function(rep_str, variant) {
  arch <- file.path(bench_dir, "A1", sprintf("A1_rep%s_v1.tar.gz", rep_str))
  pa <- as_panel_array(read_panel_table(arch))
  ab <- pa$abundance_array
  pp <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = list(method = "half_min"), normalize = "clr")
  data_in <- switch(variant,
    P0 = pp$abundance_array,
    P7 = {
      pp_raw <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                                 zero = NULL, normalize = "none")
      smooth_array(pp_raw$abundance_array, method = "rts", observation_model = "nb",
                   hierarchical = TRUE, em_iterations = 80, time_scaling = TRUE,
                   time_meta = pa$time_meta)
    },
    smooth_array(pp$abundance_array, method = switch(variant,
        P1 = "ekf", P2 = "rts", P3 = "rts", P4 = "loess", P5 = "savgol", P6 = "gaussian"),
      observation_model = "gaussian",
      qratio = if (variant == "P3") 1.0 else NULL,
      span = if (variant == "P4") 0.75 else NULL,
      halfwindow = if (variant == "P5") 2 else NULL,
      sigma = if (variant == "P6") 2 else NULL,
      breakpoints = if (variant %in% c("P4", "P5", "P6")) t0 else NULL,
      time_meta = pa$time_meta)
  )
  est <- panel_causal_estimate(feature ~ label | observation + time,
                               data = data_in, method = "scm",
                               scm_constraint = "simplex", scm_offset = TRUE,
                               intervention_index = t0,
                               observation_meta = pa$observation_meta,
                               time_meta = pa$time_meta)
  est$effects
}

variants <- c("P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7")
for (v in variants) {
  res <- lapply(reps, function(r) run_variant(r, v))
  est_mat <- sapply(res, function(e) e$estimate[match(eff, e$feature)])
  null_mat <- sapply(res, function(e) e$estimate[!(e$feature %in% eff)])
  bias <- mean(est_mat - truth_est)
  mae <- mean(abs(est_mat - truth_est))
  sign_hit <- mean(sign(est_mat) == matrix(signs, nrow = length(signs), ncol = length(reps)))
  null_abs <- mean(abs(null_mat))
  cat(sprintf("%s: bias=%+.3f MAE=%.3f sign_hit=%.2f null|est|=%.4f\n",
              v, bias, mae, sign_hit, null_abs))
  cat("   mean est per effect taxa:",
      paste(sprintf("%+.3f", rowMeans(est_mat)), collapse = " "), "\n")
  cat("   truth:                    ",
      paste(sprintf("%+.3f", truth_est), collapse = " "), "\n")
}
