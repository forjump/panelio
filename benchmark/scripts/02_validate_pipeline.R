# 02_validate_pipeline.R
# Validate the paper-default pipeline on 1 rep (A1, A11):
#   read -> as_panel_array -> preprocess(filter only) -> smooth(rts, nb, hierarchical)
#   -> panel_causal_estimate(method, simplex, offset) -> check sign / magnitude / Type I.
# Also times each method to size the full run.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

effect_features <- function(ds) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  if (nchar(tr$effect_taxa) == 0L) return(character(0L))
  idx <- as.integer(strsplit(tr$effect_taxa, ";")[[1L]])
  sprintf("taxa_%03d", idx)
}

run_estimate <- function(ds, rep_str, method, t0) {
  arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, rep_str))
  lt <- read_panel_table(arch)
  pa <- as_panel_array(lt)
  pp <- preprocess_array(pa$abundance_array,
                         filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = NULL, normalize = "none")
  sp <- smooth_array(pp$abundance_array, method = "rts",
                     observation_model = "nb", hierarchical = TRUE,
                     em_iterations = 80, time_scaling = TRUE,
                     time_meta = pa$time_meta)
  est <- panel_causal_estimate(feature ~ label | observation + time,
                               data = sp, method = method,
                               scm_constraint = "simplex", scm_offset = TRUE,
                               intervention_index = t0,
                               observation_meta = pa$observation_meta,
                               time_meta = pa$time_meta)
  list(kept = pp$kept_features, est = est)
}

cat("=== A1 rep01 (t0=4), full pipeline, all methods ===\n")
tr1 <- truth[truth$dataset == "A1", , drop = FALSE]
eff1 <- effect_features("A1")
for (m in c("scm", "did", "pooled_ols", "event_study")) {
  tt <- system.time(res <- run_estimate("A1", "01", m, as.integer(tr1$t0)))
  e <- res$est$effects
  cat(sprintf("method %-11s kept=%d  wall=%.1fs  n_signif=%d\n", m,
              length(res$kept), tt[["elapsed"]], sum(e$fdr < 0.05, na.rm = TRUE)))
  eff_rows <- e[e$feature %in% eff1, , drop = FALSE]
  print(eff_rows[, c("feature", "estimate", "standard_error", "statistic", "p_value", "fdr")])
}

cat("\n=== A11 rep01 (null), scm ===\n")
tr11 <- truth[truth$dataset == "A11", , drop = FALSE]
tt <- system.time(res11 <- run_estimate("A11", "01", "scm", as.integer(tr11$t0)))
e11 <- res11$est$effects
cat(sprintf("kept=%d  wall=%.1fs  n_signif_fdr005=%d  typeI_rate=%.3f\n",
            length(res11$kept), tt[["elapsed"]],
            sum(e11$fdr < 0.05, na.rm = TRUE),
            sum(e11$fdr < 0.05, na.rm = TRUE) / nrow(e11)))
cat("top5 smallest fdr:\n")
print(head(e11[order(e11$fdr), c("feature", "estimate", "fdr")], 5L))

cat("\n=== A12 rep01 (large effect, t0=7), scm direction check ===\n")
tr12 <- truth[truth$dataset == "A12", , drop = FALSE]
eff12 <- effect_features("A12")
signs12 <- as.integer(strsplit(tr12$effect_signs, ";")[[1L]])
tt <- system.time(res12 <- run_estimate("A12", "01", "scm", as.integer(tr12$t0)))
e12 <- res12$est$effects
eff_rows <- e12[e12$feature %in% eff12, , drop = FALSE]
print(eff_rows[, c("feature", "estimate", "p_value", "fdr")])
cat("true signs:", signs12, " est signs:", sign(eff_rows$estimate), "\n")
cat("kept effect taxa:", sum(eff12 %in% res12$kept), "/ 6\n")
