# =============================================================================
# truth/build_truth.R
# Builds the standalone ground-truth table for the A1..A12 benchmark from the
# source .rds files. The table records scenario-level truth (generative params,
# effect taxa/signs, treatment split, t0, time, per-replicate seeds) plus a
# content fingerprint (md5 of the relab_true / Y tensors) so any archive can be
# verified independently of the archives themselves.
# =============================================================================
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))
rds_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark"
truth_dir <- file.path(rds_dir, "truth")
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20260918L

fingerprint <- \(obj) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(obj, f, compress = "xz")
  unname(tools::md5sum(f))
}

rows <- vector("list", 12L)
for (i in 1:12) {
  nm <- sprintf("A%d", i)
  sim <- readRDS(file.path(rds_dir, sprintf("%s.rds", nm)))
  p <- sim$params
  rows[[i]] <- data.frame(
    dataset           = nm,
    base_seed         = base_seed,
    N                 = p$N,
    n_treated         = p$n_treated,
    n_donor           = p$N - p$n_treated,
    D                 = p$D,
    T                 = p$T,
    t0                = sim$t0,
    pre               = sim$t0 - 1L,
    post              = p$T - sim$t0 + 1L,
    n_rep             = dim(sim$Y)[4L],
    time              = paste(sim$time, collapse = ";"),
    treatment         = sprintf("first %d subjects treated, rest donor", p$n_treated),
    effect_size       = p$effect_size,
    n_effect_taxa     = length(sim$effect_taxa),
    effect_taxa       = paste(sim$effect_taxa, collapse = ";"),
    effect_signs      = paste(sim$effect_signs, collapse = ";"),
    baseline_sd       = p$baseline_sd,
    trend_sd          = p$trend_sd,
    innov_sd          = p$innov_sd,
    ref_dt            = p$ref_dt,
    rho               = p$rho,
    phi_range         = paste(p$phi_range, collapse = ":"),
    base_depth        = p$base_depth,
    depth_cv          = p$depth_cv,
    structural_zero   = p$structural_zero_prob,
    trend_diff        = p$trend_diff,
    effect_hetero     = p$effect_hetero,
    twin_trend        = p$twin_trend,
    relab_true_md5    = fingerprint(sim$relab_true),
    Y_md5             = fingerprint(sim$Y),
    rep_seeds         = paste(sim$seeds, collapse = ";"),
    stringsAsFactors  = FALSE
  )
}

truth <- do.call(rbind, rows)
out <- file.path(truth_dir, "truth_table.csv")
write.csv(truth, out, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out, nrow(truth)))
cat(sprintf("columns: %s\n", paste(names(truth), collapse = ", ")))
