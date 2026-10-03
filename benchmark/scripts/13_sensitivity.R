# 13_sensitivity.R
# Sensitivity analysis: one-variable scans of key hyperparameters, scm (simplex+offset).
# Usage: Rscript 13_sensitivity.R <config_type> <ds1> <ds2> ... <out_tag>
# config_type one of: filter | normalize | zero | band | bootstrap
#   filter    min_prevalence in {0.05, 0.1, 0.2}
#   normalize zero_fill=half_min + method in {clr, ilr, rclr, z_score, min_max, css, rle, tss, log1p}
#             (rclr uses zero=NULL; others auto-fill half_min via normalize when zeros present)
#   zero      normalize=clr + zero method in {constant, half_min, sqrt_min, multiplicative, bayesian}
#             plus rclr (no fill) as the "none" arm
#   band      spline lambda {0.5, 1, 2} and ewma alpha {0.2, 0.3, 0.5} (breakpoint=t0)
#   bootstrap samples {200, 500, 1000} x level {0.90, 0.95, 0.99} (dedup 500/0.95), reps 01-05
# Writes results/sen_<out_tag>_<config>_perfeature.csv (+ diag for band/normalize/zero).
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_type <- args[1L]
datasets <- args[seq.int(2L, length(args) - 1L)]
out_tag <- args[length(args)]
res_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results"
bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
all_reps <- sprintf("%02d", 1:30)
boot_reps <- sprintf("%02d", 1:5)
filter_on <- list(min_prevalence = 0.1, min_abundance = 0)

# ---- config tables -----------------------------------------------------------
configs <- switch(config_type,
  filter = lapply(c(0.05, 0.1, 0.2), function(p) {
    list(name = sprintf("prev_%g", p),
         filter = list(min_prevalence = p, min_abundance = 0),
         zero = list(method = "half_min"), normalize = "clr", smooth = NULL)
  }),
  normalize = lapply(c("clr", "ilr", "rclr", "z_score", "min_max", "css", "rle", "tss", "log1p"),
    function(m) {
      zero <- if (m == "rclr") NULL else list(method = "half_min")
      list(name = sprintf("norm_%s", m),
           filter = filter_on, zero = zero, normalize = m, smooth = NULL)
    }),
  zero = {
    zc <- lapply(c("constant", "half_min", "sqrt_min", "multiplicative", "bayesian"),
      function(m) list(name = sprintf("zero_%s", m),
                       filter = filter_on, zero = list(method = m),
                       normalize = "clr", smooth = NULL))
    rclr_arm <- list(name = "zero_none_rclr", filter = filter_on, zero = NULL,
                     normalize = "rclr", smooth = NULL)
    c(zc, list(rclr_arm))
  },
  band = {
    spl <- lapply(c(0.5, 1.0, 2.0), function(l) {
      list(name = sprintf("spline_l%.1f", l), filter = filter_on,
           zero = list(method = "half_min"), normalize = "clr",
           smooth = list(method = "spline", lambda = l, breakpoints = "t0"))
    })
    ewm <- lapply(c(0.2, 0.3, 0.5), function(a) {
      list(name = sprintf("ewma_a%.1f", a), filter = filter_on,
           zero = list(method = "half_min"), normalize = "clr",
           smooth = list(method = "ewma", alpha = a, breakpoints = "t0"))
    })
    c(spl, ewm)
  },
  bootstrap = {
    combos <- list(c(200, 0.95), c(500, 0.90), c(500, 0.95), c(500, 0.99), c(1000, 0.95))
    lapply(combos, function(cb) {
      list(name = sprintf("boot_s%d_l%.2f", cb[1L], cb[2L]),
           filter = filter_on, zero = list(method = "half_min"), normalize = "clr",
           smooth = NULL, bootstrap_samples = cb[1L], bootstrap_level = cb[2L])
    })
  },
  stop("unknown config_type")
)
config_names <- vapply(configs, `[[`, character(1L), "name")

run_one <- function(ds, rep_str, cfg) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, rep_str))
  tryCatch({
    pa <- as_panel_array(read_panel_table(arch))
    pp <- preprocess_array(pa$abundance_array, filter = cfg$filter,
                           zero = cfg$zero, normalize = cfg$normalize)
    data_in <- pp$abundance_array
    if (!is.null(cfg$smooth)) {
      sargs <- cfg$smooth
      sargs$abundance_array <- data_in
      sargs$time_meta <- pa$time_meta
      if ("breakpoints" %in% names(sargs) && identical(sargs$breakpoints, "t0")) {
        sargs$breakpoints <- t0
      }
      data_in <- do.call(smooth_array, sargs)
    }
    est <- panel_causal_estimate(feature ~ label | observation + time,
                                 data = data_in, method = "scm",
                                 scm_constraint = "simplex", scm_offset = TRUE,
                                 bootstrap = if ("bootstrap_samples" %in% names(cfg)) TRUE else FALSE,
                                 bootstrap_samples = cfg$bootstrap_samples %||% 500,
                                 bootstrap_level = cfg$bootstrap_level %||% 0.95,
                                 intervention_index = t0,
                                 observation_meta = pa$observation_meta,
                                 time_meta = pa$time_meta)
    e <- est$effects
    data.frame(dataset = ds, rep = rep_str, config = cfg$name, feature = e$feature,
               estimate = e$estimate, standard_error = e$standard_error,
               statistic = e$statistic, p_value = e$p_value, fdr = e$fdr,
               bootstrap_lower = if ("bootstrap_lower" %in% names(e)) e$bootstrap_lower else NA_real_,
               bootstrap_upper = if ("bootstrap_upper" %in% names(e)) e$bootstrap_upper else NA_real_,
               n_kept = length(pp$kept_features), error = "", stringsAsFactors = FALSE)
  }, error = function(err) {
    data.frame(dataset = ds, rep = rep_str, config = cfg$name,
               feature = sprintf("taxa_%03d", 1:50),
               estimate = NA_real_, standard_error = NA_real_, statistic = NA_real_,
               p_value = NA_real_, fdr = NA_real_,
               bootstrap_lower = NA_real_, bootstrap_upper = NA_real_,
               n_kept = NA_integer_, error = conditionMessage(err), stringsAsFactors = FALSE)
  })
}

flatten_df <- function(x) {
  if (is.data.frame(x)) return(list(x))
  unlist(lapply(x, flatten_df), recursive = FALSE)
}

reps_use <- if (config_type == "bootstrap") boot_reps else all_reps
rows <- lapply(datasets, function(ds) {
  lapply(reps_use, function(r) {
    lapply(configs, function(cfg) run_one(ds, r, cfg))
  })
})
perfeature <- do.call(rbind, flatten_df(rows))
pf_path <- file.path(res_dir, sprintf("sen_%s_%s_perfeature.csv", out_tag, config_type))
write.csv(perfeature, pf_path, row.names = FALSE)
cat("WROTE", pf_path, "rows:", nrow(perfeature), "\n")
cat("errors:\n")
print(table(perfeature$error[perfeature$error != ""]))
