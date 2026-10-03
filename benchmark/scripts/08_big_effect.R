# 08_big_effect.R
# How do the candidates behave on large-effect datasets A2 (T=5) and A6 (T=9)?
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
  eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
  signs <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])
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
  est_mat <- sapply(res, function(e) e$estimate[match(eff, e$feature)])
  fdr_mat <- sapply(res, function(e) e$fdr[match(eff, e$feature)])
  null_fdr <- sapply(res, function(e) e$fdr[!(e$feature %in% eff)])
  null_est <- sapply(res, function(e) e$estimate[!(e$feature %in% eff)])
  sign_hit <- mean(sign(est_mat) == matrix(signs, nrow = length(signs), ncol = 30L))
  power <- mean(fdr_mat < 0.05, na.rm = TRUE)
  bias <- mean(est_mat - truth_est)
  mae <- mean(abs(est_mat - truth_est))
  typeI <- mean(null_fdr < 0.05, na.rm = TRUE)
  cat(sprintf("%s %-4s: sign_hit=%.2f power=%.2f bias=%+.3f MAE=%.3f typeI=%.3f null|est|=%.3f\n",
              ds, variant, sign_hit, power, bias, mae, typeI, mean(abs(null_est))))
}

for (ds in c("A2", "A6")) {
  cat(sprintf("=== %s (effect 1.2) 30 reps ===\n", ds))
  for (v in c("C1", "C2", "C3", "C5")) {
    run_dataset(ds, v, "scm")
  }
  run_dataset(ds, "C1", "did")
  run_dataset(ds, "C1", "pooled_ols")
  run_dataset(ds, "C1", "event_study")
}
