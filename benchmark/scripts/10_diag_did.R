# 10_diag_did.R
# Cross-check did$estimate against hand-computed CLR DiD on A2; inspect sign
# hit distribution across 30 reps.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
tr <- truth[truth$dataset == "A2", , drop = FALSE]
t0 <- as.integer(tr$t0)
eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
signs <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])

# rep01: did estimate vs hand DiD
arch <- file.path(bench_dir, "A2", "A2_rep01_v1.tar.gz")
pa <- as_panel_array(read_panel_table(arch))
pp <- preprocess_array(pa$abundance_array,
                       filter = list(min_prevalence = 0.1, min_abundance = 0),
                       zero = list(method = "half_min"), normalize = "clr")$abundance_array
est <- panel_causal_estimate(feature ~ label | observation + time,
                             data = pp, method = "did",
                             intervention_index = t0,
                             observation_meta = pa$observation_meta,
                             time_meta = pa$time_meta)
e <- est$effects
cat("=== rep01 did effects head ===\n")
print(head(e[, c("feature", "estimate", "p_value", "fdr")], 8L))
cat("=== rep01 did on effect taxa ===\n")
er <- e[e$feature %in% eff, , drop = FALSE]
print(er[, c("feature", "estimate")])
cat("signs order (feature):", paste(er$feature, collapse = ","), "\n")
cat("truth signs:", paste(signs, collapse = ","), "\n")

# all reps: per-rep sign hit
one <- function(r) {
  arch <- file.path(bench_dir, "A2", sprintf("A2_rep%s_v1.tar.gz", r))
  pa <- as_panel_array(read_panel_table(arch))
  pp <- preprocess_array(pa$abundance_array,
                         filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = list(method = "half_min"), normalize = "clr")$abundance_array
  e <- panel_causal_estimate(feature ~ label | observation + time,
                             data = pp, method = "did",
                             intervention_index = t0,
                             observation_meta = pa$observation_meta,
                             time_meta = pa$time_meta)$effects
  e
}
res <- lapply(sprintf("%02d", 1:30), one)
hit_per_rep <- sapply(res, function(e) {
  est_v <- e$estimate[match(eff, e$feature)]
  mean(sign(est_v) == signs)
})
cat("\n=== per-rep sign hit distribution (did, A2) ===\n")
print(table(hit_per_rep))
cat("mean sign_hit:", mean(hit_per_rep), "\n")
# sanity: print rep02 effect taxa estimates
e2 <- res[[2L]]
print(e2[e2$feature %in% eff, c("feature", "estimate")])
