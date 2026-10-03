# =============================================================================
# benchmark/generate_A1_A12.R
# -----------------------------------------------------------------------------
# Generates the A1..A12 benchmark datasets from the CLSC simulation engine.
# Each dataset is a list with the same layout as the original sim_effect.rds:
#   Y             (N,D,T,R) integer count array   [raw counts, NOT normalized]
#   mu_true       (N,D,T,R) double  latent log-ratio ground truth
#   relab_true    (N,D,T,R) double  true relative abundances (closed to 1)
#   depth         (N,T,R)   double  sequencing depths
#   treatment     (N)              treated(1)/donor(0) indicator
#   structural_zero (N,D)          structural absence mask
#   effect        (D)              per-taxon effect on latent CLR scale
#   phi           (D)              NB dispersion
#   time          (T)              observation times
#   t0            (T index)        first post-intervention index
#   effect_taxa / effect_signs / params / seeds
#
# "先不做均一化": the primary observed signal is the raw negative-binomial
# count array Y; relative abundances are kept separately as ground truth.
# =============================================================================

args      <- commandArgs(trailingOnly = TRUE)
outdir    <- if (length(args) >= 1) args[1] else "."
base_seed <- as.integer(if (length(args) >= 2) args[2] else 20260918)
n_rep     <- 30L

source(file.path(outdir, "simulate_microbiome.R"))

# ---- baseline engine defaults (mirror generate_sim_data.R) ------------------
BASE <- list(
  N = 40L, n_treated = 10L, D = 50L,
  base_depth = 12000, depth_cv = 0.3,
  baseline_sd = 1.5, trend_sd = 0.25, innov_sd = 0.4, ref_dt = 7, rho = 0.4,
  phi_range = c(0.3, 3.5), structural_zero_prob = 0.15,
  n_effect_taxa = 6L, trend_diff = 0, effect_hetero = 0, twin_trend = FALSE
)

# ---- time grids --------------------------------------------------------------
# regular grids (roughly evenly spaced)
REG_5  <- c(0, 7, 14, 42, 70)
REG_7  <- c(0, 3, 7, 14, 21, 42, 70)
REG_9  <- c(0, 3, 7, 14, 21, 28, 42, 56, 70)
# irregular grids with 2-27 day gaps
IRR_9a <- c(0, 2, 9, 17, 26, 38, 51, 66, 83)
IRR_9b <- c(0, 5, 13, 22, 32, 47, 63, 80, 98)
IRR_9c <- c(0, 3, 12, 20, 27, 40, 55, 71, 89)

# ---- the twelve benchmark scenarios ------------------------------------------
# T, t0 and time define the series; other axes diversify the microbiome settings.
CONFIGS <- list(
  A1  = list(T = 5L, t0 = 4L, time = REG_5,  effect_size = 0.4),                       # baseline, short, small effect
  A2  = list(T = 5L, t0 = 4L, time = REG_5,  effect_size = 1.2),                       # baseline, short, large effect
  A3  = list(T = 7L, t0 = 5L, time = REG_7,  effect_size = 0.4),                       # medium series (pre4/post3), small
  A4  = list(T = 7L, t0 = 5L, time = REG_7,  effect_size = 1.2),                       # medium series, large
  A5  = list(T = 9L, t0 = 4L, time = REG_9,  effect_size = 0.4),                       # long series (pre3/post6), small
  A6  = list(T = 9L, t0 = 4L, time = IRR_9a, effect_size = 1.2),                       # long + irregular sampling, large
  A7  = list(T = 9L, t0 = 7L, time = IRR_9b, effect_size = 0.4),                       # long pre6/post3 + irregular, small
  A8  = list(T = 9L, t0 = 7L, time = IRR_9c, effect_size = 1.2, innov_sd = 0.8),       # + doubled innovation noise
  A9  = list(T = 9L, t0 = 4L, time = IRR_9a, effect_size = 0.4, innov_sd = 0.8,
             base_depth = 4000, structural_zero_prob = 0.4),                           # + strong zero inflation, reduced depth
  A10 = list(T = 9L, t0 = 4L, time = IRR_9b, effect_size = 1.2, innov_sd = 0.8,
             base_depth = 4000, structural_zero_prob = 0.4, trend_sd = 0.6),           # + strengthened succession
  A11 = list(T = 9L, t0 = 7L, time = IRR_9c, effect_size = 0,   innov_sd = 0.8,
             base_depth = 4000, structural_zero_prob = 0.4, trend_sd = 0.6),           # null (Type I error)
  A12 = list(T = 9L, t0 = 7L, time = IRR_9b, effect_size = 1.2,
             base_depth = 4000, structural_zero_prob = 0.4,
             trend_diff = 0.5, effect_hetero = 0.6)                                    # + confounding, heterogeneous effect
)

# ---- generation ---------------------------------------------------------------
generate_benchmark <- \(name, overrides, seeds) {
  cfg <- modifyList(BASE, overrides)
  cfg$n_effect_taxa <- 6L
  sims <- lapply(seeds, \(s) do.call(simulate_microbiome, c(cfg, list(seed = s))))
  params <- sims[[1L]]$params
  pick <- \(tag) lapply(sims, \(s) s[[tag]])
  N <- params$N; D <- params$D; T <- params$T; R <- length(seeds)
  Y <- pick("Y") |> unlist() |> array(dim = c(N, D, T, R))
  storage.mode(Y) <- "integer"   # counts are integer; keep the array integer-typed
  list(
    Y               = Y,
    mu_true         = pick("mu_true")    |> unlist() |> array(dim = c(N, D, T, R)),
    relab_true      = pick("relab_true") |> unlist() |> array(dim = c(N, D, T, R)),
    depth           = pick("depth")      |> unlist() |> array(dim = c(N, T, R)),
    treatment       = pick("treatment")  |> simplify2array(),
    structural_zero = pick("structural_zero") |> simplify2array(),
    effect          = pick("effect")     |> simplify2array(),
    phi             = pick("phi")        |> simplify2array(),
    time            = sims[[1L]]$time,
    t0              = sims[[1L]]$t0,
    effect_taxa     = if (cfg$effect_size != 0) sims[[1L]]$effect_taxa else integer(0),
    effect_signs    = if (cfg$effect_size != 0) sims[[1L]]$effect_signs else integer(0),
    effect_mult     = pick("effect_mult") |> simplify2array(),
    params          = params,
    seeds           = seeds
  )
}

cat("== CLSC benchmark generation (A1..A12) ==\n")
cat(sprintf("  base seed: %d   replicates per dataset: %d\n", base_seed, n_rep))

for (nm in names(CONFIGS)) {
  overrides <- CONFIGS[[nm]]
  seeds <- base_seed + as.integer(substring(nm, 2L)) * 100000L + seq_len(n_rep)
  sim <- generate_benchmark(nm, overrides, seeds)
  path <- file.path(outdir, sprintf("%s.rds", nm))
  saveRDS(sim, path, compress = "xz")
  p <- sim$params
  cat(sprintf(
    "  %-3s N=%d D=%d T=%d t0=%d pre/post=%d/%d | gaps=[%d,%d] | innov=%.1f trend=%.2f ",
    nm, p$N, p$D, p$T, p$t0, p$t0 - 1L, p$T - p$t0 + 1L,
    min(diff(p$time)), max(diff(p$time)), p$innov_sd, p$trend_sd))
  cat(sprintf("depth=%.0f zero=%.2f effect=%s | Y zeros=%.3f -> %s\n",
              p$base_depth, p$structural_zero_prob,
              format(p$effect_size, digits = 2), mean(sim$Y == 0), path))
}

cat("\n== done ==\n")
