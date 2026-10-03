# =============================================================================
# tests/testthat/test-causal.R
# Regression tests for panel_causal_estimate(): formula parsing, data
# ingestion (panel_array / 3D array / long table), all four estimators against
# closed-form or lm references, KKT optimality of the Fortran NNLS for the
# synthetic control, deterministic output, and the reported per-feature FDR
# and per-time effect path.
# =============================================================================

make_synth_panel <- function() {
  set.seed(42)
  nobs <- 20L
  ntime <- 9L
  nfeat <- 3L
  t0 <- 5L
  treated <- seq_len(5L)
  control <- 6:20
  treat <- logical(nobs)
  treat[treated] <- TRUE
  mu <- rnorm(nobs, 0, 0.5)
  lam <- rnorm(ntime, 0, 0.3)
  feat_base <- c(1.0, 2.0, 0.5)
  beta <- c(0.8, -0.5, 0.2)
  val <- array(0, c(nobs, nfeat, ntime))
  for (f in seq_len(nfeat)) {
    for (t in seq_len(ntime)) {
      for (i in seq_len(nobs)) {
        val[i, f, t] <- feat_base[f] + mu[i] + lam[t] +
          (treat[i] && t >= t0) * beta[f] + rnorm(1, 0, 0.2)
      }
    }
  }
  dimnames(val) <- list(paste0("s", seq_len(nobs)), paste0("f", seq_len(nfeat)),
                        sprintf("time_%03d", seq_len(ntime)))
  obs <- data.frame(observation_id = dimnames(val)[[1]],
                    label = factor(ifelse(treat, "treat", "ctrl"),
                                   levels = c("ctrl", "treat")))
  tim <- data.frame(time_identifier = dimnames(val)[[3]],
                    elapsed_time = seq_len(ntime),
                    delta_t = 1,
                    is_intervention_timepoint = seq_len(ntime) >= t0)
  list(val = val, obs = obs, tim = tim, t0 = t0, beta = beta,
       treat = treat, treated = treated, control = control)
}

test_that("panel_causal_estimate returns the full reported structure", {
  s <- make_synth_panel()
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  r <- panel_causal_estimate(feature ~ label | subject + time, data = pa)
  expect_s3_class(r, "panel_causal_result")
  expect_named(r, c("call", "formula", "method", "t0", "t0_time_identifier",
                    "treatment_level", "control_level", "effects",
                    "time_effects", "weights", "weight_sum", "placebo",
                    "bootstrap", "pre_fit"))
  expect_equal(nrow(r$effects), 3L)
  expect_true(all(c("estimate", "standard_error", "statistic", "p_value",
                    "fdr", "effect_range") %in% names(r$effects)))
  expect_true(all(r$effects$fdr >= 0 & r$effects$fdr <= 1))
  expect_true(nrow(r$time_effects) >= 3L)
  expect_equal(r$t0, 5L)
  expect_equal(r$treatment_level, "treat")
  expect_equal(r$control_level, "ctrl")
})

test_that("data is accepted as panel_array, 3D array, and long table", {
  s <- make_synth_panel()
  es <- exp(s$val)
  sums <- apply(es, c(1L, 3L), sum)                     # (O, T)
  comp <- aperm(sweep(aperm(es, c(2L, 1L, 3L)), c(2L, 3L), sums, "/"),
                c(2L, 1L, 3L))
  dimnames(comp) <- dimnames(s$val)
  pa <- structure(list(abundance_array = comp, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  r1 <- panel_causal_estimate(feature ~ label | subject + time, data = pa)
  # 3D array with explicit meta
  r2 <- panel_causal_estimate(feature ~ label | subject + time, data = comp,
                              observation_meta = s$obs, time_meta = s$tim)
  # long table (observation, time_identifier, elapsed_time, feature,
  # abundance_value, label, is_intervention)
  tt <- rep(seq_len(9L), each = 20L)
  long <- data.frame(
    observation = rep(dimnames(s$val)[[1]], times = 9L),
    time_identifier = rep(dimnames(s$val)[[3]], each = 20L),
    elapsed_time = tt,
    feature = rep(dimnames(s$val)[[2]], each = 180L),
    abundance_value = as.vector(aperm(comp, c(1L, 3L, 2L))),
    label = factor(rep(rep(as.character(s$obs$label), 9L)), levels = c("ctrl", "treat")),
    is_intervention = tt >= s$t0
  )
  r3 <- panel_causal_estimate(feature ~ label | subject + time, data = long)
  expect_equal(r1$effects$estimate, r2$effects$estimate)
  expect_equal(r1$effects$estimate, r3$effects$estimate)
})

test_that("did matches the closed-form difference-in-differences", {
  s <- make_synth_panel()
  f <- 1L
  pre <- seq_len(s$t0 - 1L)
  post <- s$t0:9L
  ref <- (mean(s$val[s$treated, f, post]) - mean(s$val[s$treated, f, pre])) -
    (mean(s$val[s$control, f, post]) - mean(s$val[s$control, f, pre]))
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  r <- panel_causal_estimate(feature ~ label | subject + time, data = pa,
                             method = "did")
  expect_lt(abs(r$effects$estimate[f] - ref), 1e-9)
})

test_that("pooled_ols matches lm with subject and time fixed effects", {
  s <- make_synth_panel()
  f <- 1L
  nobs <- 20L
  ntime <- 9L
  y <- as.vector(s$val[, f, ])
  subj <- factor(rep(dimnames(s$val)[[1]], ntime))
  timi <- factor(rep(dimnames(s$val)[[3]], each = nobs))
  trt <- rep(s$treat, ntime)
  post <- rep(seq_len(ntime) >= s$t0, each = nobs)
  tp <- as.numeric(trt & post)                 # additive Treat x Post dummy
  ref <- coef(lm(y ~ subj + timi + tp))[["tp"]]
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  r <- panel_causal_estimate(feature ~ label | subject + time, data = pa,
                             method = "pooled_ols")
  expect_lt(abs(r$effects$estimate[f] - ref), 1e-6)
})

test_that("event_study coefficients match the joint K-dummy regression", {
  s <- make_synth_panel()
  f <- 1L
  ntime <- 9L
  nobs <- 20L
  t0 <- s$t0
  y <- as.vector(s$val[, f, ])
  subj <- factor(rep(dimnames(s$val)[[1]], ntime))
  timi <- factor(rep(dimnames(s$val)[[3]], each = nobs))
  rel <- seq_len(ntime) - t0
  ks <- sort(unique(rel))
  ks <- ks[ks != -1L]
  ev_rows <- outer(rep(rel, each = nobs), ks, "==") * rep(as.numeric(s$treat), ntime)
  Xev <- matrix(as.numeric(ev_rows), nrow = nobs * ntime, ncol = length(ks))
  fit <- lm(y ~ subj + timi + Xev)
  ref <- coef(fit)[grepl("Xev", names(coef(fit)))]
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  r <- panel_causal_estimate(feature ~ label | subject + time, data = pa,
                             method = "event_study")
  ev <- r$time_effects[r$time_effects$feature == dimnames(s$val)[[2]][f], ]
  expect_equal(nrow(ev), length(ks))
  expect_lt(max(abs(ev$estimate - ref)), 1e-6)
})

test_that("scm weights satisfy the NNLS KKT optimality conditions", {
  s <- make_synth_panel()
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  sc <- panel_causal_estimate(feature ~ label | subject + time, data = pa,
                              method = "scm")
  expect_equal(dim(sc$weights), c(5L, 15L, 3L))
  for (ff in seq_len(3L)) {
    ctrl_m <- s$val[s$control, ff, ]
    tret_m <- s$val[s$treated, ff, ]
    wi <- as.matrix(sc$weights[, , ff])
    for (i in seq_len(5L)) {
      A <- t(ctrl_m[, seq_len(4L)])
      b <- tret_m[i, seq_len(4L)]
      r <- as.vector(A %*% wi[i, ] - b)
      grad <- as.vector(t(A) %*% r)
      act <- wi[i, ] > 1e-9
      gmax_act <- if (any(act)) max(abs(grad[act])) else 0
      gmin_in <- if (any(!act)) min(grad[!act]) else 0
      expect_lt(gmax_act, 1e-6)
      expect_gt(gmin_in, -1e-6)
    }
  }
})

test_that("scm output is deterministic across repeated calls", {
  s <- make_synth_panel()
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  a <- panel_causal_estimate(feature ~ label | subject + time, data = pa, method = "scm")
  b <- panel_causal_estimate(feature ~ label | subject + time, data = pa, method = "scm")
  expect_equal(a$effects$estimate, b$effects$estimate)
  expect_equal(a$weights, b$weights)
})

test_that("invalid inputs are rejected", {
  s <- make_synth_panel()
  expect_error(panel_causal_estimate(value ~ label | subject + time,
                                     data = s$val))  # array without meta
  pa <- structure(list(abundance_array = s$val, observation_meta = s$obs,
                       time_meta = s$tim), class = "panel_array")
  expect_error(panel_causal_estimate(feature ~ label | subject + time, data = pa,
                                     method = "bogus"))
  bad <- pa$abundance_array
  dimnames(bad) <- NULL
  bad_pa <- structure(list(abundance_array = bad,
                           observation_meta = s$obs,
                           time_meta = s$tim),
                      class = "panel_array")
  expect_error(panel_causal_estimate(feature ~ label | subject + time,
                                     data = bad_pa))
})
