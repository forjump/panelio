# 11_sim_shard.R
# Shared simulation shard for the panelio benchmark evaluation.
# Usage: Rscript 11_sim_shard.R <ds1> <ds2> ... <out_tag>
# For each dataset x rep x method runs the default pipeline
#   preprocess(filter=0.1, zero=half_min, normalize=clr)  [no smoothing]
#   did / pooled_ols / event_study: normal-theory p + BH-FDR
#   scm: simplex + offset, placebo permutation + BH-FDR, pre-fit diagnostics
# Writes results/sim_<out_tag>_perfeature.csv and results/sim_<out_tag>_diag.csv.
# Vectorized only (lapply / vapply / apply); no for/while loops.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
datasets <- args[seq_len(length(args) - 1L)]
out_tag <- args[length(args)]
res_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results"
bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
all_reps <- sprintf("%02d", 1:30)
methods <- c("did", "pooled_ols", "event_study", "scm")

prep_clr <- function(pa) {
  preprocess_array(pa$abundance_array,
                   filter = list(min_prevalence = 0.1, min_abundance = 0),
                   zero = list(method = "half_min"), normalize = "clr")$abundance_array
}

run_rep_method <- function(ds, rep_str, method) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, rep_str))
  tryCatch({
    pa <- as_panel_array(read_panel_table(arch))
    pp <- prep_clr(pa)
    est <- panel_causal_estimate(feature ~ label | observation + time,
                                 data = pp, method = method,
                                 scm_constraint = "simplex", scm_offset = TRUE,
                                 placebo = identical(method, "scm"),
                                 pre_fit = identical(method, "scm"),
                                 intervention_index = t0,
                                 observation_meta = pa$observation_meta,
                                 time_meta = pa$time_meta)
    e <- est$effects
    pfit <- est$pre_fit
    pre_rmse <- pre_snr <- NA_real_
    if (!is.null(pfit)) {
      pre_rmse <- vapply(e$feature, function(f) {
        vals <- pfit$pre_rmse[pfit$feature == f]
        if (length(vals) > 0L) mean(vals) else NA_real_
      }, numeric(1L))
      pre_snr <- vapply(e$feature, function(f) {
        vals <- pfit$pre_snr[pfit$feature == f]
        if (length(vals) > 0L) mean(vals) else NA_real_
      }, numeric(1L))
    }
    wsum <- if (is.null(est$weight_sum)) NA_real_ else mean(est$weight_sum)
    data.frame(
      dataset = ds, rep = rep_str, method = method, feature = e$feature,
      estimate = e$estimate, standard_error = e$standard_error,
      statistic = e$statistic, p_value = e$p_value, fdr = e$fdr,
      placebo_p_value = if ("placebo_p_value" %in% names(e)) e$placebo_p_value else NA_real_,
      placebo_fdr = if ("placebo_fdr" %in% names(e)) e$placebo_fdr else NA_real_,
      pre_rmse = pre_rmse, pre_snr = pre_snr, weight_sum = wsum,
      n_kept = dim(pp)[2L], error = "", stringsAsFactors = FALSE)
  }, error = function(err) {
    data.frame(dataset = ds, rep = rep_str, method = method,
               feature = sprintf("taxa_%03d", 1:50),
               estimate = NA_real_, standard_error = NA_real_, statistic = NA_real_,
               p_value = NA_real_, fdr = NA_real_, placebo_p_value = NA_real_,
               placebo_fdr = NA_real_, pre_rmse = NA_real_, pre_snr = NA_real_,
               weight_sum = NA_real_, n_kept = NA_integer_, error = conditionMessage(err),
               stringsAsFactors = FALSE)
  })
}

flatten_df <- function(x) {
  if (is.data.frame(x)) return(list(x))
  unlist(lapply(x, flatten_df), recursive = FALSE)
}

rows <- lapply(datasets, function(ds) {
  lapply(all_reps, function(r) {
    lapply(methods, function(m) run_rep_method(ds, r, m))
  })
})
perfeature <- do.call(rbind, flatten_df(rows))
diag <- unique(perfeature[, c("dataset", "rep", "method", "weight_sum", "n_kept", "error")])

perfeature_path <- file.path(res_dir, sprintf("sim_%s_perfeature.csv", out_tag))
diag_path <- file.path(res_dir, sprintf("sim_%s_diag.csv", out_tag))
write.csv(perfeature, perfeature_path, row.names = FALSE)
write.csv(diag, diag_path, row.names = FALSE)
cat("WROTE", perfeature_path, "rows:", nrow(perfeature), "\n")
cat("WROTE", diag_path, "rows:", nrow(diag), "\n")
cat("errors:\n")
print(table(diag$error[diag$error != ""]))
