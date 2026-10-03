# 12_ablation.R
# Ablation analysis: toggle one paper-depth feature at a time, scm estimator.
# Usage: Rscript 12_ablation.R <ds1> <ds2> ... <out_tag>
# Variants (each changes exactly one thing vs `base` unless coupled by design):
#   base            CLR(half_min) + scm(simplex, offset=TRUE)          [= main pipeline]
#   ab_nnls         scm constraint nnls
#   ab_no_offset    scm offset FALSE
#   ab_smooth_spline  spline(lambda=1, breakpoint=t0) on CLR + scm
#   ab_smooth_ewma    ewma(alpha=0.3, breakpoint=t0) on CLR + scm
#   ab_kalman_gauss   rts gaussian observation (qratio=0.1) on CLR + scm
#   ab_no_filter    filter off (all 50 taxa)
#   ab_rclr_nofill  rclr (zero-robust, no zero fill) instead of clr+half_min
#   ab_no_normalize raw relative-abundance scale (normalize="none")
#   ab_nb_hier      rts NB hierarchical (paper-depth count model) on raw + scm
#   ab_nb_nohier    rts NB non-hierarchical (EM off) on raw + scm
# Writes results/abl_<out_tag>_perfeature.csv and results/abl_<out_tag>_diag.csv.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
datasets <- args[seq_len(length(args) - 1L)]
out_tag <- args[length(args)]
res_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results"
bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
all_reps <- sprintf("%02d", 1:30)

# ---- variant spec table ------------------------------------------------------
filter_on <- list(min_prevalence = 0.1, min_abundance = 0)
variant_specs <- list(
  base = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
              smooth = NULL, constraint = "simplex", offset = TRUE),
  ab_nnls = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
                 smooth = NULL, constraint = "nnls", offset = TRUE),
  ab_no_offset = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
                      smooth = NULL, constraint = "simplex", offset = FALSE),
  ab_smooth_spline = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
                          smooth = list(method = "spline", lambda = 1.0, breakpoints = "t0"),
                          constraint = "simplex", offset = TRUE),
  ab_smooth_ewma = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
                        smooth = list(method = "ewma", alpha = 0.3, breakpoints = "t0"),
                        constraint = "simplex", offset = TRUE),
  ab_kalman_gauss = list(filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
                         smooth = list(method = "rts", observation_model = "gaussian", qratio = 0.1,
                                       time_scaling = FALSE),
                         constraint = "simplex", offset = TRUE),
  ab_no_filter = list(filter = NULL, zero = list(method = "half_min"), normalize = "clr",
                      smooth = NULL, constraint = "simplex", offset = TRUE),
  ab_rclr_nofill = list(filter = filter_on, zero = NULL, normalize = "rclr",
                        smooth = NULL, constraint = "simplex", offset = TRUE),
  ab_no_normalize = list(filter = filter_on, zero = NULL, normalize = "none",
                         smooth = NULL, constraint = "simplex", offset = TRUE),
  ab_nb_hier = list(filter = filter_on, zero = NULL, normalize = "none",
                    smooth = list(method = "rts", observation_model = "nb", hierarchical = TRUE,
                                  em_iterations = 80, time_scaling = TRUE),
                    constraint = "simplex", offset = TRUE),
  ab_nb_nohier = list(filter = filter_on, zero = NULL, normalize = "none",
                      smooth = list(method = "rts", observation_model = "nb", hierarchical = FALSE,
                                    time_scaling = TRUE),
                      constraint = "simplex", offset = TRUE)
)
variants <- names(variant_specs)

run_one <- function(ds, rep_str, vname) {
  spec <- variant_specs[[vname]]
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, rep_str))
  tryCatch({
    pa <- as_panel_array(read_panel_table(arch))
    pp <- preprocess_array(pa$abundance_array, filter = spec$filter,
                           zero = spec$zero, normalize = spec$normalize)
    data_in <- pp$abundance_array
    kept <- pp$kept_features
    if (!is.null(spec$smooth)) {
      sargs <- spec$smooth
      sargs$abundance_array <- data_in
      sargs$time_meta <- pa$time_meta
      if ("breakpoints" %in% names(sargs) && identical(sargs$breakpoints, "t0")) {
        sargs$breakpoints <- t0
      }
      data_in <- do.call(smooth_array, sargs)
    }
    est <- panel_causal_estimate(feature ~ label | observation + time,
                                 data = data_in, method = "scm",
                                 scm_constraint = spec$constraint,
                                 scm_offset = spec$offset,
                                 intervention_index = t0,
                                 observation_meta = pa$observation_meta,
                                 time_meta = pa$time_meta)
    e <- est$effects
    wsum <- if (is.null(est$weight_sum)) NA_real_ else mean(est$weight_sum)
    data.frame(dataset = ds, rep = rep_str, variant = vname, feature = e$feature,
               estimate = e$estimate, standard_error = e$standard_error,
               statistic = e$statistic, p_value = e$p_value, fdr = e$fdr,
               weight_sum = wsum, n_kept = length(kept), error = "",
               stringsAsFactors = FALSE)
  }, error = function(err) {
    data.frame(dataset = ds, rep = rep_str, variant = vname,
               feature = sprintf("taxa_%03d", 1:50),
               estimate = NA_real_, standard_error = NA_real_, statistic = NA_real_,
               p_value = NA_real_, fdr = NA_real_, weight_sum = NA_real_,
               n_kept = NA_integer_, error = conditionMessage(err), stringsAsFactors = FALSE)
  })
}

flatten_df <- function(x) {
  if (is.data.frame(x)) return(list(x))
  unlist(lapply(x, flatten_df), recursive = FALSE)
}

rows <- lapply(datasets, function(ds) {
  lapply(all_reps, function(r) {
    lapply(variants, function(v) run_one(ds, r, v))
  })
})
perfeature <- do.call(rbind, flatten_df(rows))
diag <- unique(perfeature[, c("dataset", "rep", "variant", "weight_sum", "n_kept", "error")])
pf_path <- file.path(res_dir, sprintf("abl_%s_perfeature.csv", out_tag))
dg_path <- file.path(res_dir, sprintf("abl_%s_diag.csv", out_tag))
write.csv(perfeature, pf_path, row.names = FALSE)
write.csv(diag, dg_path, row.names = FALSE)
cat("WROTE", pf_path, "rows:", nrow(perfeature), "\n")
cat("errors:\n")
print(table(diag$error[diag$error != ""]))
