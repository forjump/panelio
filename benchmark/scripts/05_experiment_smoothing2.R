# 05_experiment_smoothing2.R
# Second round: finer smoothing grid for SCM (simplex + offset) on A1 reps 01-05.
# Variants on CLR-preprocessed trajectories:
#   S1 rts qr=0.1 no-ts   S2 rts qr=1 no-ts   S3 rts qr=10 no-ts
#   S4 ekf qr=1 no-ts     S5 loess span=.75   S6 savgol hw=2
#   S7 gauss sigma=2      S8 MA hw=2          S9 none (no smoothing)
# Also run DID on CLR as reference (should work well) and A2 rep01 S2/S9 sanity.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

clr_prep <- function(pa) {
  preprocess_array(pa$abundance_array,
                   filter = list(min_prevalence = 0.1, min_abundance = 0),
                   zero = list(method = "half_min"), normalize = "clr")$abundance_array
}

smooth_in <- function(pp, variant, t0, pa) {
  if (variant == "S9") return(pp)
  smooth_array(pp, method = switch(variant,
      S1 = "rts", S2 = "rts", S3 = "rts", S4 = "ekf",
      S6 = "savgol", S7 = "gaussian", S8 = "moving_average", S10 = "ewma", S11 = "spline"),
    observation_model = "gaussian",
    qratio = switch(variant, S1 = 0.1, S2 = 1.0, S3 = 10.0, S4 = 1.0, NULL),
    time_scaling = FALSE,
    halfwindow = if (variant == "S6") 2 else if (variant == "S8") 2 else NULL,
    sigma = if (variant == "S7") 2 else NULL,
    alpha = if (variant == "S10") 0.3 else NULL,
    lambda = if (variant == "S11") 1.0 else NULL,
    breakpoints = if (variant %in% c("S6", "S7", "S8", "S10", "S11")) t0 else NULL,
    time_meta = pa$time_meta)
}

run_reps <- function(ds, reps, variant, method = "scm") {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
  signs <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])
  truth_est <- signs * as.numeric(tr$effect_size)
  one <- function(r) {
    arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, r))
    pa <- as_panel_array(read_panel_table(arch))
    data_in <- smooth_in(clr_prep(pa), variant, t0, pa)
    est <- panel_causal_estimate(feature ~ label | observation + time,
                                 data = data_in, method = method,
                                 scm_constraint = "simplex", scm_offset = TRUE,
                                 intervention_index = t0,
                                 observation_meta = pa$observation_meta,
                                 time_meta = pa$time_meta)
    est$effects
  }
  res <- lapply(reps, one)
  list(ds = ds, eff = eff, truth_est = truth_est, res = res)
}

summarize <- function(run, label) {
  eff <- run$eff
  truth_est <- run$truth_est
  est_mat <- sapply(run$res, function(e) e$estimate[match(eff, e$feature)])
  null_mat <- sapply(run$res, function(e) e$estimate[!(e$feature %in% eff)])
  bias <- mean(est_mat - truth_est)
  mae <- mean(abs(est_mat - truth_est))
  sign_hit <- mean(sign(est_mat) == matrix(signs, nrow = length(signs), ncol = ncol(est_mat)))
  cat(sprintf("%-4s: bias=%+.3f MAE=%.3f sign_hit=%.2f null|est|=%.4f  mean_est=%s\n",
              label, bias, mae, sign_hit, mean(abs(null_mat)),
              paste(sprintf("%+.2f", rowMeans(est_mat)), collapse = " ")))
}

reps <- sprintf("%02d", 1:5)
signs <- as.integer(strsplit(truth$effect_signs[truth$dataset == "A1"], ";")[[1L]])
cat("=== A1 reps 01-05, SCM simplex+offset, CLR ===\n")
for (v in c("S1", "S2", "S3", "S4", "S6", "S7", "S8", "S10", "S11", "S9")) {
  run <- run_reps("A1", reps, v, "scm")
  summarize(run, v)
}
cat("\n=== A1 reps 01-05, DID on CLR (reference) ===\n")
run_did <- run_reps("A1", reps, "S9", "did")
summarize(run_did, "DID")

cat("\n=== A2 (effect 1.2) reps 01-05, SCM S2 vs S9 vs DID ===\n")
tr2 <- truth[truth$dataset == "A2", , drop = FALSE]
signs2 <- as.integer(strsplit(tr2$effect_signs, ";")[[1L]])
for (v in c("S2", "S9")) {
  run <- run_reps("A2", reps, v, "scm")
  summarize(run, paste0("A2-", v))
}
run_did2 <- run_reps("A2", reps, "S9", "did")
summarize(run_did2, "A2-DID")
