# 07_choose_pipeline.R
# Decide the paper-default pipeline on full A1 (30 reps) and A11 (null, Type I).
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
all_reps <- sprintf("%02d", 1:30)

prep_clr <- function(pa) {
  preprocess_array(pa$abundance_array,
                   filter = list(min_prevalence = 0.1, min_abundance = 0),
                   zero = list(method = "half_min"), normalize = "clr")$abundance_array
}
smooth_variant <- function(pp, variant, t0, pa) {
  if (variant == "C1") return(pp)
  smooth_array(pp, method = switch(variant,
      C2 = "ewma", C3 = "spline", C4 = "gaussian", C5 = "moving_average"),
    alpha = if (variant == "C2") 0.3 else NULL,
    lambda = if (variant == "C3") 1.0 else NULL,
    sigma = if (variant == "C4") 2 else NULL,
    halfwindow = if (variant == "C5") 2 else NULL,
    breakpoints = t0, time_meta = pa$time_meta)
}

run_dataset <- function(ds, variant, method) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  eff <- if (nchar(tr$effect_taxa) > 0L) {
    sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
  } else character(0L)
  signs <- if (nchar(tr$effect_signs) > 0L) {
    as.integer(strsplit(tr$effect_signs, ";")[[1L]])
  } else integer(0L)
  truth_est <- signs * as.numeric(tr$effect_size)
  one <- function(r) {
    arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, r))
    pa <- as_panel_array(read_panel_table(arch))
    data_in <- smooth_variant(prep_clr(pa), variant, t0, pa)
    est <- panel_causal_estimate(feature ~ label | observation + time,
                                 data = data_in, method = method,
                                 scm_constraint = "simplex", scm_offset = TRUE,
                                 intervention_index = t0,
                                 observation_meta = pa$observation_meta,
                                 time_meta = pa$time_meta)
    est$effects
  }
  res <- lapply(all_reps, one)
  list(ds = ds, eff = eff, signs = signs, truth_est = truth_est, res = res)
}

summarize <- function(run) {
  ds <- run$ds
  eff <- run$eff
  signs <- run$signs
  truth_est <- run$truth_est
  n_rep <- length(run$res)
  if (length(eff) > 0L) {
    est_mat <- sapply(run$res, function(e) e$estimate[match(eff, e$feature)])
    fdr_mat <- sapply(run$res, function(e) e$fdr[match(eff, e$feature)])
    null_fdr <- sapply(run$res, function(e) e$fdr[!(e$feature %in% eff)])
    sign_hit <- mean(sign(est_mat) == matrix(signs, nrow = length(signs), ncol = n_rep))
    power <- mean(fdr_mat < 0.05, na.rm = TRUE)
    bias <- mean(est_mat - truth_est)
    mae <- mean(abs(est_mat - truth_est))
    typeI <- mean(null_fdr < 0.05, na.rm = TRUE)
    cat(sprintf("%s: sign_hit=%.2f power(fdr)=%.2f bias=%+.3f MAE=%.3f typeI=%.3f\n",
                ds, sign_hit, power, bias, mae, typeI))
  } else {
    fdr_all <- sapply(run$res, function(e) e$fdr)
    typeI <- mean(fdr_all < 0.05, na.rm = TRUE)
    cat(sprintf("%s (null): typeI=%.3f  mean_est_null=%.4f\n", ds, typeI,
                mean(abs(sapply(run$res, function(e) e$estimate)))))
  }
}

cat("=== A1 (effect 0.4), scm candidates, 30 reps ===\n")
for (v in c("C1", "C2", "C3", "C4", "C5")) {
  run <- run_dataset("A1", v, "scm")
  cat(sprintf("%-4s ", v)); summarize(run)
}
cat("D1   "); summarize(run_dataset("A1", "C1", "did"))

cat("\n=== A11 (null), type I for leading candidates ===\n")
for (v in c("C1", "C2", "C3")) {
  run <- run_dataset("A11", v, "scm")
  cat(sprintf("%-4s ", v)); summarize(run)
}
cat("D1   "); summarize(run_dataset("A11", "C1", "did"))
