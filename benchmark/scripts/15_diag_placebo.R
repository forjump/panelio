# 15_diag_placebo.R
# Reproduce the scm + placebo + pre_fit failure with full traceback.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))
bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
arch <- file.path(bench_dir, "A1", "A1_rep01_v1.tar.gz")
pa <- as_panel_array(read_panel_table(arch))
pp <- preprocess_array(pa$abundance_array,
                       filter = list(min_prevalence = 0.1, min_abundance = 0),
                       zero = list(method = "half_min"), normalize = "clr")
cat("prep ok, dims:", dim(pp$abundance_array), "\n")
cat("=== scm no placebo/pre_fit ===\n")
e0 <- panel_causal_estimate(feature ~ label | observation + time, data = pp$abundance_array,
                            method = "scm", scm_constraint = "simplex", scm_offset = TRUE,
                            intervention_index = 4L,
                            observation_meta = pa$observation_meta, time_meta = pa$time_meta)
cat("ok, effects rows:", nrow(e0$effects), " weights class:", class(e0$weights),
    " dim:", if (is.array(e0$weights)) dim(e0$weights) else length(e0$weights), "\n")
cat("=== scm + placebo only ===\n")
e1 <- tryCatch({
  panel_causal_estimate(feature ~ label | observation + time, data = pp$abundance_array,
                        method = "scm", scm_constraint = "simplex", scm_offset = TRUE,
                        placebo = TRUE, intervention_index = 4L,
                        observation_meta = pa$observation_meta, time_meta = pa$time_meta)
}, error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n")
  traceback()
  NULL
})
cat("=== scm + pre_fit only ===\n")
e2 <- tryCatch({
  panel_causal_estimate(feature ~ label | observation + time, data = pp$abundance_array,
                        method = "scm", scm_constraint = "simplex", scm_offset = TRUE,
                        pre_fit = TRUE, intervention_index = 4L,
                        observation_meta = pa$observation_meta, time_meta = pa$time_meta)
}, error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n")
  traceback()
  NULL
})
