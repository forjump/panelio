# simulate_microbiome.R
# ==============================================================================
# Compositional Longitudinal Synthetic Control (CLSC)
# Simulation engine for short, sparse, compositional longitudinal microbiome
#
# Generative model (latent log-ratio scale):
#   mu[i,d,t] = a[i,d] + b[i,d]*s[t] + eps[i,d,t]     (s = scaled window time)
#   a ~ N(0, baseline_sd)                     inter-individual baseline offsets
#   b ~ N(0, trend_sd)                        individual spontaneous succession
#   eps: AR(1), autocorrelation rho, innovation sd scaled by sqrt(dt / ref_dt)
#   causal effect: treated units, t >= t0, fixed effect_taxa + effect_signs
#   relab[i,.,t] = softmax(mu); structural zeros forced to 0
#   lambda[i,d,t] = depth[i,t] * relab[i,d,t]
#   y[i,d,t] ~ NB(size = phi_d, mu = lambda)   (overdispersion / dropout)
# ==============================================================================

simulate_microbiome <- \(
    N = 40L,                       # subjects
    n_treated = 10L,               # treated subjects (the rest are donors)
    D = 50L,                       # taxa
    T = 5L,                        # time points (pre = t0-1, post = T-t0+1)
    t0 = 4L,                       # first post-intervention index
    time = c(0, 7, 14, 42, 70),    # observation times (days)
    base_depth = 12000,            # typical sequencing depth per sample
    depth_cv = 0.3,                # log-normal CV of depth
    baseline_sd = 1.5,             # sd of inter-individual baseline offsets
    trend_sd = 0.25,               # sd of individual drift over the window
    innov_sd = 0.4,                # AR(1) innovation sd (per ref_dt)
    ref_dt = 7,                    # reference time step for innovation sd
    rho = 0.4,                     # AR(1) autocorrelation
    phi_range = c(0.3, 3.5),       # NB dispersion range (per taxon)
    effect_size = 0.8,             # causal effect magnitude (log-ratio scale)
    n_effect_taxa = 6L,            # number of taxa carrying a real effect
    effect_taxa = NULL,            # fixed effect taxa indices
    effect_signs = NULL,           # fixed signs for the effect taxa
    structural_zero_prob = 0.15,   # P(taxon absent in a subject)
    trend_diff = 0,                # extra post/pre drift for treated units (confounding)
    effect_hetero = 0,             # per-treated-unit effect heterogeneity (sd of multiplier)
    twin_trend = FALSE,            # treated units copy a random donor's baseline+trend
    seed = NULL) {                 # RNG seed

  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(time) == T, t0 >= 2, t0 <= T,
            n_treated >= 1, n_treated <= N, phi_range[1] > 0)

  treatment        <- c(rep(1L, n_treated), rep(0L, N - n_treated))
  phi              <- runif(D, phi_range[1], phi_range[2])
  a                <- matrix(rnorm(N * D, 0, baseline_sd), N, D)
  b                <- matrix(rnorm(N * D, 0, trend_sd),   N, D)
  structural_zero  <- matrix(runif(N * D) < structural_zero_prob, N, D)
  open             <- !structural_zero

  # --- optional: treated units are "twins" of random donors (matchable pre-trend)
  if (twin_trend) {
    don <- which(treatment == 0L); tr <- which(treatment == 1L)
    twin <- sample(don, length(tr), replace = TRUE)
    a[tr, ] <- a[twin, ]; b[tr, ] <- b[twin, ]
  }

  # --- AR(1) temporal noise, vectorized via a lower-triangular rho matrix ----
  s      <- (time - time[1L]) / max(time - time[1L])
  sd_step <- innov_sd * sqrt(c(1, diff(time)) / ref_dt)
  u      <- array(rnorm(N * D * T), c(N, D, T)) * rep(sd_step, each = N * D)
  R_ar   <- outer(seq_len(T), seq_len(T), \(tk, sk) ifelse(tk >= sk, rho^(tk - sk), 0))
  ar_ma  <- \(R) \(m) m %*% R                       # pipe-safe matrix product
  eps    <- u |> matrix(nrow = N * D) |> ar_ma(t(R_ar))() |> array(dim = c(N, D, T))

  # --- ground-truth latent states: baseline + drift + AR(1) noise -------------
  mu <- array(a, c(N, D, T)) +
        array(b, c(N, D, T)) * rep(s, each = N * D) + eps

  treated <- which(treatment == 1L)

  # --- optional: treated units carry an extra time drift (time-varying confound)
  if (trend_diff != 0) {
    mu[treated, , ] <- mu[treated, , ] + rep(trend_diff * s, each = length(treated) * D)
  }

  # --- causal effect: fixed taxa & signs, treated units, t >= t0 ----------------
  effect <- rep(0, D)
  mult   <- rep(1, length(treated))
  if (effect_size != 0) {
    if (is.null(effect_taxa)) effect_taxa <- sort(sample.int(D, min(n_effect_taxa, D)))
    if (is.null(effect_signs)) effect_signs <- sample(c(-1L, 1L), length(effect_taxa), replace = TRUE)
    stopifnot(length(effect_signs) == length(effect_taxa))
    effect[effect_taxa] <- effect_size * effect_signs
    if (t0 <= T) {
      if (effect_hetero == 0) {
        mu[treated, , t0:T] <- mu[treated, , t0:T] + rep(effect, each = length(treated))
      } else {
        mult <- 1 + effect_hetero * stats::rnorm(length(treated))
        for (k in seq_along(treated)) {
          mu[treated[k], effect_taxa, t0:T] <- mu[treated[k], effect_taxa, t0:T] +
            rep(effect_size * mult[k] * effect_signs, times = length(t0:T))
        }
      }
    }
  }

  # --- softmax -> relative abundances; structural zeros stay at 0 ----------------
  e     <- exp(sweep(mu, c(1L, 3L), apply(mu, c(1L, 3L), max), "-"))
  relab <- sweep(e, c(1L, 3L), apply(e, c(1L, 3L), sum), "/")
  relab[rep(structural_zero, T)] <- 0
  # renormalize so each composition is closed on the simplex (sum = 1 over all D)
  scl   <- apply(relab, c(1L, 3L), sum); scl[scl == 0] <- 1
  relab <- sweep(relab, c(1L, 3L), scl, "/")

  # --- NB observations -----------------------------------------------------------------
  depth  <- matrix(rlnorm(N * T, log(base_depth), depth_cv), N, T)
  lambda <- relab |> sweep(c(1L, 3L), depth, "*")
  open_ar   <- array(rep(open, T), c(N, D, T))
  phi_ar    <- array(rep(phi, each = N), c(N, D, T))
  Y         <- array(0L, c(N, D, T))
  Y[open_ar] <- rnbinom(sum(open_ar), size = phi_ar[open_ar], mu = lambda[open_ar])

  list(
    Y               = Y,                 # (N,D,T) observed counts
    mu_true         = mu,                # (N,D,T) latent log-ratio states
    relab_true      = relab,             # (N,D,T) true relative abundances
    depth           = depth,             # (N,T)   sequencing depths
    structural_zero = structural_zero,   # (N,D)   structural absence mask
    treatment       = treatment,         # (N)     treated/donor indicator
    time            = time,              # (T)
    t0              = t0,
    effect          = effect,            # (D)     per-taxon effect (log-ratio)
    effect_taxa     = effect_taxa,
    effect_signs    = effect_signs,
    effect_mult     = mult,              # (n_treated) per-unit effect multiplier
    phi             = phi,               # (D)     NB dispersion
    seed            = seed,
    params = list(
      N = N, n_treated = n_treated, D = D, T = T, t0 = t0, time = time,
      base_depth = base_depth, depth_cv = depth_cv,
      baseline_sd = baseline_sd, trend_sd = trend_sd,
      innov_sd = innov_sd, ref_dt = ref_dt, rho = rho,
      phi_range = phi_range, effect_size = effect_size,
      n_effect_taxa = n_effect_taxa, structural_zero_prob = structural_zero_prob,
      trend_diff = trend_diff, effect_hetero = effect_hetero, twin_trend = twin_trend
    )
  )
}
