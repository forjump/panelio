# 06_diagnose_locreg.R
# Which locreg smoothers stay finite and preserve the intervention jump on
# short post segments (A1: T=5 post=2; A3: T=7 post=3)?
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

check <- function(ds, rep_str) {
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  t0 <- as.integer(tr$t0)
  arch <- file.path(bench_dir, ds, sprintf("%s_rep%s_v1.tar.gz", ds, rep_str))
  pa <- as_panel_array(read_panel_table(arch))
  pp <- preprocess_array(pa$abundance_array,
                         filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = list(method = "half_min"), normalize = "clr")$abundance_array
  eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
  cat(sprintf("=== %s rep%s (t0=%d, post=%d) ===\n", ds, rep_str, t0,
              nrow(pa$time_meta) - t0 + 1L))
  variants <- list(
    gaussian = list(method = "gaussian", sigma = 2),
    savgol = list(method = "savgol", halfwindow = 2),
    ma = list(method = "moving_average", halfwindow = 2),
    ewma = list(method = "ewma", alpha = 0.3),
    spline = list(method = "spline", lambda = 1.0),
    loess = list(method = "loess", span = 0.75)
  )
  for (nm in names(variants)) {
    args <- c(variants[[nm]], list(breakpoints = t0, time_meta = pa$time_meta))
    sp <- tryCatch(do.call(smooth_array, c(list(abundance_array = pp), args)),
                   error = function(e) NA)
    if (is.logical(sp) && is.na(sp)) {
      cat(sprintf("  %-10s ERROR\n", nm))
      next
    }
    fin <- mean(is.finite(sp))
    pre_post_jump <- mean(sp[, eff, t0, drop = FALSE]) -
                     mean(sp[, eff, t0 - 1L, drop = FALSE])
    cat(sprintf("  %-10s finite=%.3f range=[%.3f, %.3f]  effect-jump(pre->post)=%+.3f\n",
                nm, fin, min(sp), max(sp), pre_post_jump))
  }
}

check("A1", "01")
check("A3", "01")
check("A5", "01")
