# =============================================================================
# R/causal.R
# Per-feature causal effect estimation for a balanced longitudinal panel.
# The entry point panel_causal_estimate() accepts a three-part formula
#   response ~ treatment | subject + time
# (fixest style), a panel array / panel_array / long-table data source, and one
# of several causal estimators: did, pooled_ols (two-way fixed effects),
# event_study (dynamic event-time effects) and scm (synthetic control).
#
# Treatment assignment is the label column; the intervention onset t0 is the
# first intervention timepoint. Two-way fixed-effect demeaning and the
# synthetic-control non-negative least squares run inside the Fortran kernels
# under src/panelio_causal.f90; the R layer only dispatches one .Fortran call
# per operation and carries no for/while loop. Per-feature p-values and
# Benjamini-Hochberg FDR, plus a per-time (time-fluctuation) effect path, are
# always reported.
# =============================================================================

# ---- formula parsing ----------------------------------------------------------
# Parse a three-part formula "Y ~ X | FE1 + FE2" into its named components.
parse_causal_formula <- \(formula) {
  text <- gsub(" ", "", paste(deparse(formula), collapse = ""))
  if (!grepl("~", text, fixed = TRUE)) {
    panel_abort("A causal formula must contain '~', e.g. feature ~ label | subject + time.")
  }
  lhs <- sub("~.*", "", text)
  rhs_all <- sub("^.*~", "", text)
  rhs_parts <- strsplit(rhs_all, "\\|", fixed = FALSE)[[1L]]
  treatment <- rhs_parts[1L]
  fixed <- character(0L)
  if (length(rhs_parts) > 1L) {
    fixed <- strsplit(rhs_parts[2L], "+", fixed = TRUE)[[1L]]
  }
  list(response = lhs, treatment = treatment, fixed = fixed)
}

# Map formula names onto the canonical panelio variables. Aliases are accepted
# so the user can write subject/observation, time/time_identifier, and so on.
resolve_causal_variable <- \(name) {
  switch(
    tolower(name),
    feature = "feature",
    abundance = "abundance_value",
    abundance_value = "abundance_value",
    value = "abundance_value",
    label = "label",
    treatment = "label",
    arm = "label",
    group = "label",
    subject = "observation",
    observation = "observation",
    obs = "observation",
    individual = "observation",
    time = "time_identifier",
    time_identifier = "time_identifier",
    timepoint = "time_identifier",
    t = "time_identifier",
    name
  )
}

validate_causal_formula <- \(parsed) {
  response <- resolve_causal_variable(parsed$response)
  treatment <- resolve_causal_variable(parsed$treatment)
  if (!(response %in% c("feature", "abundance_value"))) {
    panel_abort(sprintf(
      "The causal response '%s' must resolve to a feature abundance.",
      parsed$response
    ))
  }
  if (treatment != "label") {
    panel_abort(sprintf(
      "The causal treatment '%s' must resolve to the label column.",
      parsed$treatment
    ))
  }
  fixed_resolved <- vapply(parsed$fixed, resolve_causal_variable, character(1L))
  if (length(fixed_resolved) > 0L &&
        !all(fixed_resolved %in% c("observation", "time_identifier"))) {
    panel_abort("The fixed-effect terms after '|' must resolve to subject and/or time.")
  }
  invisible(TRUE)
}

# ---- data ingestion -----------------------------------------------------------
# Normalize the data argument (panel_array, 3D array, or long table) onto the
# canonical panel structure used by every estimator.
as_causal_panel <- \(data, observation_meta, time_meta) {
  if (inherits(data, "panel_array")) {
    return(data)
  }
  if (is.data.frame(data)) {
    return(as_panel_array(data))
  }
  if (is.array(data) && length(dim(data)) == 3L) {
    if (is.null(observation_meta) || is.null(time_meta)) {
      panel_abort(
        "When data is a 3D array, both observation_meta (with a label column) ",
        "and time_meta (with is_intervention_timepoint) must be supplied."
      )
    }
    return(structure(
      list(
        abundance_array = data,
        observation_meta = observation_meta,
        time_meta = time_meta
      ),
      class = "panel_array"
    ))
  }
  panel_abort(
    "data must be a panel_array, a three-dimensional abundance array, ",
    "or a long-format panel table."
  )
}

resolve_t0 <- \(time_meta, intervention_index, ntime) {
  if (!is.null(intervention_index)) {
    t0 <- as.integer(intervention_index)
    if (t0 < 2L || t0 > ntime) {
      panel_abort("intervention_index must satisfy 2 <= t0 <= ntime (need at least one pre time).")
    }
    return(t0)
  }
  if (!is.null(time_meta) && "is_intervention_timepoint" %in% names(time_meta)) {
    onset <- which(as.logical(time_meta$is_intervention_timepoint))
    if (length(onset) > 0L) {
      return(as.integer(min(onset)))
    }
  }
  panel_abort(
    "Cannot determine the intervention onset t0. Supply intervention_index ",
    "or time_meta with an is_intervention_timepoint column."
  )
}

resolve_treatment <- \(observation_meta, treatment_level, control_level) {
  label <- observation_meta$label
  label_chr <- as.character(label)
  # Prefer the declared factor level order (usually c("ctrl","treat")), falling
  # back to order of first appearance for character labels.
  label_levels <- if (is.factor(label)) levels(label) else unique(label_chr)
  if (length(label_levels) < 2L) {
    panel_abort("The label column must contain at least two levels (treated and control).")
  }
  if (is.null(treatment_level)) {
    if (length(label_levels) == 2L) {
      treatment_level <- label_levels[2L]
    } else {
      panel_abort(
        "label has more than two levels; supply treatment_level explicitly."
      )
    }
  }
  if (!(treatment_level %in% label_levels)) {
    panel_abort(sprintf("treatment_level '%s' is not a label level.", treatment_level))
  }
  control_set <- if (is.null(control_level)) {
    setdiff(label_levels, treatment_level)
  } else {
    as.character(control_level)
  }
  if (length(setdiff(control_set, label_levels)) > 0L) {
    panel_abort("control_level contains levels absent from the label column.")
  }
  treat <- label_chr %in% treatment_level
  control <- label_chr %in% control_set
  list(treat = treat, control = control,
       treatment_level = treatment_level, control_level = control_set)
}

# ---- shared helpers -----------------------------------------------------------
# Mean over a selected set of observations and time indices, per feature.
panel_mean <- \(value, obs_idx, time_idx) {
  if (length(obs_idx) == 0L || length(time_idx) == 0L) {
    return(rep(NA_real_, dim(value)[3L]))
  }
  apply(value[obs_idx, time_idx, , drop = FALSE], 3L, mean)
}

two_way_demean <- \(array3d) {
  d <- dim(array3d)
  dm <- .Fortran(
    "demean2",
    mat = array3d, nobs = d[1L], ntime = d[2L],
    nfeat = d[3L], maxiter = as.integer(CAUSAL_MAX_DEMEAN_ITER),
    tol = as.numeric(CAUSAL_DEMEAN_TOLERANCE)
  )$mat
  dimnames(dm) <- dimnames(array3d)
  dm
}

# Two-sided p-value from the standard-normal tail.
p_from_t <- \(statistic) 2 * pnorm(-abs(statistic))

# ---- estimators ---------------------------------------------------------------
# Each estimator receives the canonical panel list and returns a named list
#   effects       per-feature summary data.frame
#   time_effects  per-feature x per-time effect path data.frame
# The R layer performs only vectorized matrix operations (no for/while loops);
# the heavy per-feature iteration lives in the Fortran kernels.

did_estimator <- \(pan) {
  value <- pan$value
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  t0 <- pan$t0
  treated_idx <- which(pan$treat)
  control_idx <- which(pan$control)
  n_tr <- length(treated_idx)
  n_co <- length(control_idx)
  pre <- seq_len(t0 - 1L)
  post <- t0:ntime

  pre_tr <- panel_mean(value, treated_idx, pre)
  post_tr <- panel_mean(value, treated_idx, post)
  pre_co <- panel_mean(value, control_idx, pre)
  post_co <- panel_mean(value, control_idx, post)
  estimate <- (post_tr - pre_tr) - (post_co - pre_co)

  pre_tr_obs <- apply(value[treated_idx, pre, , drop = FALSE], c(1L, 3L), mean)
  post_tr_obs <- apply(value[treated_idx, post, , drop = FALSE], c(1L, 3L), mean)
  pre_co_obs <- apply(value[control_idx, pre, , drop = FALSE], c(1L, 3L), mean)
  post_co_obs <- apply(value[control_idx, post, , drop = FALSE], c(1L, 3L), mean)
  var_tr <- apply(post_tr_obs - pre_tr_obs, 2L, var)
  var_co <- apply(post_co_obs - pre_co_obs, 2L, var)
  standard_error <- sqrt(var_tr / n_tr + var_co / n_co)
  statistic <- estimate / standard_error
  p_value <- p_from_t(statistic)

  # per-time treated-control difference path (time fluctuation)
  tr_series <- apply(value[treated_idx, , , drop = FALSE], c(2L, 3L), mean)
  co_series <- apply(value[control_idx, , , drop = FALSE], c(2L, 3L), mean)
  time_diff <- tr_series - co_series
  tr_var <- apply(value[treated_idx, , , drop = FALSE], c(2L, 3L), var)
  co_var <- apply(value[control_idx, , , drop = FALSE], c(2L, 3L), var)
  time_se <- sqrt(tr_var / n_tr + co_var / n_co)
  time_stat <- time_diff / time_se
  time_p <- p_from_t(time_stat)
  effect_range <- apply(time_diff, 2L, function(x) diff(range(x)))

  list(
    effects = data.frame(
      feature = pan$feature_ids,
      method = rep(CAUSAL_DID, nfeat),
      estimate = estimate,
      standard_error = standard_error,
      statistic = statistic,
      p_value = p_value,
      mean_pre_treated = pre_tr,
      mean_post_treated = post_tr,
      mean_pre_control = pre_co,
      mean_post_control = post_co,
      n_treated = rep(n_tr, nfeat),
      n_control = rep(n_co, nfeat),
      effect_range = effect_range,
      stringsAsFactors = FALSE
    ),
    time_effects = build_time_effects(pan, time_diff, time_se, time_p,
                                      method = CAUSAL_DID)
  )
}

pooled_ols_estimator <- \(pan) {
  value <- pan$value
  nobs <- pan$nobs
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  t0 <- pan$t0
  pre <- seq_len(t0 - 1L)
  post <- t0:ntime

  # design column Treat x Post on the balanced [obs, time] grid
  design <- outer(pan$treat, as.numeric(seq_len(ntime) >= t0))
  design_arr <- array(as.numeric(design), c(nobs, ntime, 1L))
  design_dm <- two_way_demean(design_arr)
  dvec <- as.vector(design_dm[, , 1L])

  value_dm <- two_way_demean(value)
  ymat <- matrix(value_dm, nrow = nobs * ntime, ncol = nfeat)

  dd <- sum(dvec ^ 2)
  estimate <- as.vector(crossprod(dvec, ymat) / dd)
  resid <- ymat - outer(dvec, estimate)
  rss <- colSums(resid ^ 2)
  df <- nobs * ntime - (nobs + ntime - 1L) - 1L
  standard_error <- sqrt((rss / df) / dd)
  statistic <- estimate / standard_error
  p_value <- p_from_t(statistic)

  treated_idx <- which(pan$treat)
  control_idx <- which(pan$control)
  pre_tr <- panel_mean(value, treated_idx, pre)
  post_tr <- panel_mean(value, treated_idx, post)
  pre_co <- panel_mean(value, control_idx, pre)
  post_co <- panel_mean(value, control_idx, post)
  tr_series <- apply(value[treated_idx, , , drop = FALSE], c(2L, 3L), mean)
  co_series <- apply(value[control_idx, , , drop = FALSE], c(2L, 3L), mean)
  time_diff <- tr_series - co_series
  tr_var <- apply(value[treated_idx, , , drop = FALSE], c(2L, 3L), var)
  co_var <- apply(value[control_idx, , , drop = FALSE], c(2L, 3L), var)
  time_se <- sqrt(tr_var / length(treated_idx) + co_var / length(control_idx))
  time_stat <- time_diff / time_se
  time_p <- p_from_t(time_stat)
  effect_range <- apply(time_diff, 2L, function(x) diff(range(x)))

  list(
    effects = data.frame(
      feature = pan$feature_ids,
      method = rep(CAUSAL_POOLED, nfeat),
      estimate = estimate,
      standard_error = standard_error,
      statistic = statistic,
      p_value = p_value,
      mean_pre_treated = pre_tr,
      mean_post_treated = post_tr,
      mean_pre_control = pre_co,
      mean_post_control = post_co,
      n_treated = rep(length(treated_idx), nfeat),
      n_control = rep(length(control_idx), nfeat),
      effect_range = effect_range,
      stringsAsFactors = FALSE
    ),
    time_effects = build_time_effects(pan, time_diff, time_se, time_p,
                                      method = CAUSAL_POOLED)
  )
}

event_study_estimator <- \(pan) {
  value <- pan$value
  nobs <- pan$nobs
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  t0 <- pan$t0
  base_rel <- CAUSAL_EVENT_BASE_REL

  rel <- seq_len(ntime) - t0
  ks <- sort(unique(rel))
  ks <- ks[ks != base_rel]
  nk <- length(ks)

  # rows of the (nobs*ntime) stack are obs-fastest, time-slowest
  treat_r <- rep(pan$treat, times = ntime)
  rel_r <- rep(rel, each = nobs)
  # (rel_r == ks) * treat_r as an explicit outer comparison avoids recycling.
  xmat <- matrix(as.numeric(outer(rel_r, ks, "==") * treat_r),
                 nrow = nobs * ntime, ncol = nk)
  dim(xmat) <- c(nobs, ntime, nk)
  xmat_dm <- two_way_demean(xmat)
  xm <- matrix(xmat_dm, nrow = nobs * ntime, ncol = nk)

  value_dm <- two_way_demean(value)
  ymat <- matrix(value_dm, nrow = nobs * ntime, ncol = nfeat)

  gram <- crossprod(xm) + CAUSAL_GRAM_RIDGE * diag(nk)
  inv_gram <- solve(gram)
  beta_kf <- inv_gram %*% crossprod(xm, ymat)   # nk x nfeat
  resid <- ymat - xm %*% beta_kf
  rss <- colSums(resid ^ 2)
  df <- nobs * ntime - (nobs + ntime - 1L) - nk
  se_kf <- sqrt(outer(diag(inv_gram), rss / df))  # nk x nfeat
  stat_kf <- beta_kf / se_kf
  p_kf <- p_from_t(stat_kf)

  treated_idx <- which(pan$treat)
  control_idx <- which(pan$control)
  pre <- seq_len(t0 - 1L)
  post <- t0:ntime
  pre_tr <- panel_mean(value, treated_idx, pre)
  post_tr <- panel_mean(value, treated_idx, post)
  pre_co <- panel_mean(value, control_idx, pre)
  post_co <- panel_mean(value, control_idx, post)

  # summary effect = mean of post-intervention event-time coefficients (k >= 0)
  post_k <- which(ks >= 0L)
  estimate <- if (length(post_k) > 0L) {
    colMeans(beta_kf[post_k, , drop = FALSE])
  } else {
    rep(NA_real_, nfeat)
  }
  se_post <- if (length(post_k) > 0L) {
    sqrt(colMeans(se_kf[post_k, , drop = FALSE] ^ 2))
  } else {
    rep(NA_real_, nfeat)
  }
  statistic <- estimate / se_post
  p_value <- p_from_t(statistic)
  effect_range <- apply(beta_kf, 2L, function(x) diff(range(x)))

  # per-relative-time effect path
  time_effects <- data.frame(
    feature = rep(pan$feature_ids, each = nk),
    time_identifier = rep(pan$time_ids[ks + t0], times = nfeat),
    relative_time = rep(ks, times = nfeat),
    method = rep(CAUSAL_EVENT, nk * nfeat),
    estimate = as.vector(beta_kf),
    standard_error = as.vector(se_kf),
    statistic = as.vector(stat_kf),
    p_value = as.vector(p_kf),
    stringsAsFactors = FALSE
  )

  list(
    effects = data.frame(
      feature = pan$feature_ids,
      method = rep(CAUSAL_EVENT, nfeat),
      estimate = estimate,
      standard_error = se_post,
      statistic = statistic,
      p_value = p_value,
      mean_pre_treated = pre_tr,
      mean_post_treated = post_tr,
      mean_pre_control = pre_co,
      mean_post_control = post_co,
      n_treated = rep(length(treated_idx), nfeat),
      n_control = rep(length(control_idx), nfeat),
      effect_range = effect_range,
      stringsAsFactors = FALSE
    ),
    time_effects = time_effects
  )
}

scm_estimator <- \(pan, constraint = CAUSAL_SCM_NNLS, offset = CAUSAL_SCM_OFFSET_DEFAULT) {
  value <- pan$value
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  t0 <- pan$t0
  treated_idx <- which(pan$treat)
  control_idx <- which(pan$control)
  n_tr <- length(treated_idx)
  n_co <- length(control_idx)
  tpost <- ntime - t0 + 1L

  ctrl <- value[control_idx, , , drop = FALSE]
  tret <- value[treated_idx, , , drop = FALSE]

  if (constraint == CAUSAL_SCM_SIMPLEX) {
    out <- .Fortran(
      "scm_simplex",
      ctrl = ctrl,
      tret = tret,
      nobs = pan$nobs,
      ncontrol = as.integer(n_co),
      ntreated = as.integer(n_tr),
      nfeat = as.integer(nfeat),
      ntime = as.integer(ntime),
      t0 = as.integer(t0),
      offset = as.integer(isTRUE(offset)),
      maxiter = as.integer(CAUSAL_SCM_MAX_ITER),
      tol = as.numeric(CAUSAL_SCM_TOLERANCE),
      weights = array(0, c(n_tr, n_co, nfeat)),
      gap = array(0, c(n_tr, tpost, nfeat)),
      att = numeric(nfeat),
      wsum = array(0, c(n_tr, nfeat))
    )
    gap <- out$gap
    att <- out$att
    weight_sum <- out$wsum
  } else {
    out <- .Fortran(
      "scm_effects",
      ctrl = ctrl,
      tret = tret,
      nobs = pan$nobs,
      ncontrol = as.integer(n_co),
      ntreated = as.integer(n_tr),
      nfeat = as.integer(nfeat),
      ntime = as.integer(ntime),
      t0 = as.integer(t0),
      weights = array(0, c(n_tr, n_co, nfeat)),
      gap = array(0, c(n_tr, tpost, nfeat)),
      att = numeric(nfeat)
    )
    gap <- out$gap
    att <- out$att
    weight_sum <- apply(out$weights, c(1L, 3L), sum)
  }

  pre <- seq_len(t0 - 1L)
  post <- t0:ntime
  pre_tr <- panel_mean(value, treated_idx, pre)
  post_tr <- panel_mean(value, treated_idx, post)
  pre_co <- panel_mean(value, control_idx, pre)
  post_co <- panel_mean(value, control_idx, post)

  unit_att <- apply(gap, c(1L, 3L), mean)          # (n_tr, nfeat) post gap per unit
  var_att <- apply(unit_att, 2L, var)
  standard_error <- sqrt(var_att / n_tr)
  statistic <- att / standard_error
  p_value <- p_from_t(statistic)

  time_gap <- apply(gap, c(2L, 3L), mean)          # (tpost, nfeat)
  time_var <- apply(gap, c(2L, 3L), var)
  time_se <- sqrt(time_var / n_tr)
  time_stat <- time_gap / time_se
  time_p <- p_from_t(time_stat)
  effect_range <- apply(time_gap, 2L, function(x) diff(range(x)))

  time_effects <- data.frame(
    feature = rep(pan$feature_ids, each = tpost),
    time_identifier = rep(pan$time_ids[post], times = nfeat),
    relative_time = rep(post - t0, times = nfeat),
    method = rep(CAUSAL_SCM, tpost * nfeat),
    estimate = as.vector(time_gap),
    standard_error = as.vector(time_se),
    statistic = as.vector(time_stat),
    p_value = as.vector(time_p),
    stringsAsFactors = FALSE
  )

  list(
    effects = data.frame(
      feature = pan$feature_ids,
      method = rep(CAUSAL_SCM, nfeat),
      estimate = att,
      standard_error = standard_error,
      statistic = statistic,
      p_value = p_value,
      mean_pre_treated = pre_tr,
      mean_post_treated = post_tr,
      mean_pre_control = pre_co,
      mean_post_control = post_co,
      n_treated = rep(n_tr, nfeat),
      n_control = rep(n_co, nfeat),
      effect_range = effect_range,
      stringsAsFactors = FALSE
    ),
    time_effects = time_effects,
    weights = out$weights,
    weight_sum = weight_sum
  )
}

# ---- inference depth (paper: placebo / bootstrap / pre-fit diagnostics) ------
# Lower-level synthetic-control point estimate: returns att (nfeat), weights
# (n_tr, n_co, nfeat) and post gap (n_tr, tpost, nfeat) for an arbitrary
# treated/control split. Reused by the placebo permutation and the subject-level
# bootstrap so those layers recompute only the point estimate, not the whole
# result object.
scm_point_att <- \(pan, treated_idx, control_idx, constraint, offset) {
  n_tr <- length(treated_idx)
  n_co <- length(control_idx)
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  t0 <- pan$t0
  tpost <- ntime - t0 + 1L
  ctrl <- pan$value[control_idx, , , drop = FALSE]
  tret <- pan$value[treated_idx, , , drop = FALSE]
  if (constraint == CAUSAL_SCM_SIMPLEX) {
    out <- .Fortran(
      "scm_simplex",
      ctrl = ctrl, tret = tret, nobs = pan$nobs,
      ncontrol = as.integer(n_co), ntreated = as.integer(n_tr),
      nfeat = as.integer(nfeat), ntime = as.integer(ntime), t0 = as.integer(t0),
      offset = as.integer(isTRUE(offset)),
      maxiter = as.integer(CAUSAL_SCM_MAX_ITER),
      tol = as.numeric(CAUSAL_SCM_TOLERANCE),
      weights = array(0, c(n_tr, n_co, nfeat)),
      gap = array(0, c(n_tr, tpost, nfeat)),
      att = numeric(nfeat), wsum = array(0, c(n_tr, nfeat))
    )
  } else {
    out <- .Fortran(
      "scm_effects",
      ctrl = ctrl, tret = tret, nobs = pan$nobs,
      ncontrol = as.integer(n_co), ntreated = as.integer(n_tr),
      nfeat = as.integer(nfeat), ntime = as.integer(ntime), t0 = as.integer(t0),
      weights = array(0, c(n_tr, n_co, nfeat)),
      gap = array(0, c(n_tr, tpost, nfeat)), att = numeric(nfeat)
    )
  }
  list(att = out$att, weights = out$weights, gap = out$gap)
}

# Per-feature point estimate for a (possibly resampled) panel; used by the
# subject-level bootstrap. Fully vectorized: matrix/apply only, no for/while.
panel_point_estimate <- \(pan, method, constraint, offset) {
  if (method == CAUSAL_DID) {
    pre <- seq_len(pan$t0 - 1L)
    post <- pan$t0:pan$ntime
    tr <- which(pan$treat)
    co <- which(pan$control)
    (panel_mean(pan$value, tr, post) - panel_mean(pan$value, tr, pre)) -
      (panel_mean(pan$value, co, post) - panel_mean(pan$value, co, pre))
  } else if (method == CAUSAL_POOLED) {
    design <- outer(pan$treat, as.numeric(seq_len(pan$ntime) >= pan$t0))
    design_arr <- array(as.numeric(design), c(pan$nobs, pan$ntime, 1L))
    dvec <- as.vector(two_way_demean(design_arr)[, , 1L])
    ymat <- matrix(two_way_demean(pan$value),
                   nrow = pan$nobs * pan$ntime, ncol = pan$nfeat)
    as.vector(crossprod(dvec, ymat) / sum(dvec ^ 2))
  } else if (method == CAUSAL_EVENT) {
    rel <- seq_len(pan$ntime) - pan$t0
    ks <- sort(unique(rel))
    ks <- ks[ks != CAUSAL_EVENT_BASE_REL]
    treat_r <- rep(pan$treat, times = pan$ntime)
    rel_r <- rep(rel, each = pan$nobs)
    xmat <- matrix(as.numeric(outer(rel_r, ks, "==") * treat_r),
                   nrow = pan$nobs * pan$ntime, ncol = length(ks))
    dim(xmat) <- c(pan$nobs, pan$ntime, length(ks))
    xm <- matrix(two_way_demean(xmat),
                 nrow = pan$nobs * pan$ntime, ncol = length(ks))
    ymat <- matrix(two_way_demean(pan$value),
                   nrow = pan$nobs * pan$ntime, ncol = pan$nfeat)
    inv_gram <- solve(crossprod(xm) + CAUSAL_GRAM_RIDGE * diag(length(ks)))
    beta_kf <- inv_gram %*% crossprod(xm, ymat)
    post_k <- which(ks >= 0L)
    if (length(post_k) > 0L) {
      colMeans(beta_kf[post_k, , drop = FALSE])
    } else {
      rep(NA_real_, pan$nfeat)
    }
  } else {
    scm_point_att(pan, which(pan$treat), which(pan$control),
                  constraint, offset)$att
  }
}

# Subject-level bootstrap percentile confidence interval on the per-feature
# estimate: resample treated and control observations with replacement and
# recompute the point estimate n_samples times (vectorized via vapply).
bootstrap_ci <- \(pan, method, constraint, offset, n_samples, level) {
  tr_idx <- which(pan$treat)
  co_idx <- which(pan$control)
  ntr <- length(tr_idx)
  nco <- length(co_idx)
  alpha <- 1 - level
  draws <- vapply(seq_len(n_samples), function(b) {
    s_tr <- sample(tr_idx, ntr, replace = TRUE)
    s_co <- sample(co_idx, nco, replace = TRUE)
    pan_b <- pan
    pan_b$value <- pan$value[c(s_tr, s_co), , , drop = FALSE]
    pan_b$treat <- c(rep(TRUE, ntr), rep(FALSE, nco))
    pan_b$control <- !pan_b$treat
    pan_b$nobs <- ntr + nco
    panel_point_estimate(pan_b, method, constraint, offset)
  }, numeric(pan$nfeat))                       # (nfeat, n_samples)
  data.frame(
    feature = pan$feature_ids,
    bootstrap_lower = apply(draws, 1L, quantile, probs = alpha, na.rm = TRUE),
    bootstrap_upper = apply(draws, 1L, quantile, probs = 1 - alpha, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Placebo permutation test: every control unit in turn acts as a pseudo-treated
# unit with the remaining control units as its donor pool. The per-feature
# p-value is the fraction of |placebo effect| >= |real effect| across control
# units (a Fisher-style placebo test), followed by Benjamini-Hochberg FDR.
placebo_test <- \(pan, constraint, offset) {
  tr <- which(pan$treat)
  co <- which(pan$control)
  real <- scm_point_att(pan, tr, co, constraint, offset)$att
  nco <- length(co)
  plc <- vapply(seq_len(nco), function(d) {
    scm_point_att(pan, co[d], co[-d], constraint, offset)$att
  }, numeric(pan$nfeat))                       # (nfeat, nco)
  p_value <- apply(abs(plc) >= abs(real), 1L,
                   function(x) (sum(x) + 1L) / (nco + 1L))
  data.frame(
    feature = pan$feature_ids,
    placebo_p_value = p_value,
    placebo_fdr = p.adjust(p_value, method = "BH"),
    stringsAsFactors = FALSE
  )
}

# Abadie pre-fit diagnostics per treated unit and feature: pre-intervention
# root-mean-square error of the synthetic control against the treated trajectory
# and the signal-to-noise ratio mean(|pre gap|) / pre-RMSE. Only meaningful for
# scm, whose weights drive the pre counterfactual.
add_pre_fit <- \(pan, weights, offset) {
  tr <- which(pan$treat)
  co <- which(pan$control)
  tpre <- pan$t0 - 1L
  ntr <- length(tr)
  nfeat <- pan$nfeat
  X <- pan$value[co, seq_len(tpre), , drop = FALSE]    # (n_co, tpre, nfeat)
  Y <- pan$value[tr, seq_len(tpre), , drop = FALSE]    # (n_tr, tpre, nfeat)
  cf <- array(vapply(seq_len(nfeat), function(f) weights[, , f] %*% X[, , f],
                     matrix(0, ntr, tpre)), c(ntr, tpre, nfeat))
  if (isTRUE(offset)) {
    ym <- apply(Y, c(1L, 3L), mean)                    # (ntr, nfeat)
    xm <- apply(X, c(1L, 3L), mean)                    # (n_co, nfeat)
    wmu <- vapply(seq_len(nfeat), function(f) weights[, , f] %*% xm[, f],
                  numeric(ntr))
    mu <- ym - wmu
    cf <- cf + array(rep(as.vector(mu), times = tpre), c(ntr, tpre, nfeat))
  }
  gap_pre <- Y - cf
  pre_rmse <- sqrt(apply(gap_pre ^ 2, c(1L, 3L), mean))
  pre_snr <- apply(abs(gap_pre), c(1L, 3L), mean) / (pre_rmse + CAUSAL_SNR_EPS)
  data.frame(
    feature = rep(pan$feature_ids, each = ntr),
    treated_unit = rep(pan$observation_ids[tr], times = nfeat),
    pre_rmse = as.vector(pre_rmse),
    pre_snr = as.vector(pre_snr),
    stringsAsFactors = FALSE
  )
}

# Assemble the generic per-time difference path for methods without a native
# per-time coefficient (did / pooled_ols).
build_time_effects <- \(pan, time_diff, time_se, time_p, method) {
  ntime <- pan$ntime
  nfeat <- pan$nfeat
  data.frame(
    feature = rep(pan$feature_ids, each = ntime),
    time_identifier = rep(pan$time_ids, times = nfeat),
    relative_time = rep(seq_len(ntime) - pan$t0, times = nfeat),
    method = rep(method, ntime * nfeat),
    estimate = as.vector(time_diff),
    standard_error = as.vector(time_se),
    statistic = as.vector(time_diff / time_se),
    p_value = as.vector(time_p),
    stringsAsFactors = FALSE
  )
}

# ---- entry point --------------------------------------------------------------
#' Estimate per-feature causal effects of a panel intervention
#'
#' Runs a causal estimator on a balanced longitudinal microbiome panel. The
#' response is a feature abundance; the treatment is the \code{label} column;
#' the intervention onset is the first intervention timepoint. Estimators:
#' \code{"did"} (difference-in-differences, default), \code{"pooled_ols"}
#' (two-way fixed-effect OLS), \code{"event_study"} (dynamic event-time
#' effects) and \code{"scm"} (synthetic control). Two-way fixed-effect demeaning
#' and the synthetic-control non-negative least squares are compiled Fortran
#' kernels; the R layer dispatches one \code{.Fortran} call per operation.
#'
#' @param formula a three-part formula \code{response ~ treatment | subject + time}
#'   (fixest style). \code{response} resolves to a feature abundance,
#'   \code{treatment} to the label column, and the fixed-effect terms after
#'   \code{|} to subject and/or time (aliases accepted).
#' @param data a panel_array, a three-dimensional abundance array
#'   [observation, feature, time_identifier], or a long-format panel table.
#' @param method one of \code{CAUSAL_METHODS}: \code{"did"}, \code{"pooled_ols"},
#'   \code{"event_study"}, \code{"scm"}.
#' @param treatment_level label level treated by the intervention; when NULL it
#'   defaults to the second label level when exactly two levels exist.
#' @param control_level label level(s) used as control; NULL pools every level
#'   other than the treatment level.
#' @param intervention_index optional 1-based index of the intervention onset;
#'   used when time_meta carries no intervention flag.
#' @param scm_constraint for \code{method = "scm"}, \code{"nnls"} (classic
#'   non-negative least squares, w >= 0) or \code{"simplex"} (paper-depth
#'   exactly-convex weights, w >= 0 and sum(w) = 1, solved by an exact
#'   active-set KKT solver so the counterfactual stays in the convex hull of donors).
#' @param scm_offset logical; when \code{scm_constraint = "simplex"}, include a
#'   scalar level offset (centred log-ratio translation) absorbing baseline
#'   differences between the treated unit and the convex donor mixture.
#' @param placebo logical; when TRUE run the paper's placebo permutation test:
#'   every control unit in turn acts as a pseudo-treated unit and the per-feature
#'   p-value is the fraction of |placebo effect| >= |real effect|, with
#'   Benjamini-Hochberg FDR appended as placebo_p_value / placebo_fdr.
#' @param bootstrap logical; when TRUE compute a subject-level bootstrap
#'   percentile confidence interval on the per-feature estimate (resampling
#'   treated and control observations with replacement).
#' @param bootstrap_samples number of bootstrap resamples (default 500).
#' @param bootstrap_level confidence level for the bootstrap interval (default
#'   0.95).
#' @param pre_fit logical; when TRUE (and method = "scm") report Abadie pre-fit
#'   diagnostics per treated unit and feature: pre-intervention RMSE of the
#'   synthetic control and its signal-to-noise ratio mean(|pre gap|)/pre-RMSE.
#' @param observation_meta optional observation metadata with a label column,
#'   needed when data is a plain 3D array.
#' @param time_meta optional time metadata with an is_intervention_timepoint
#'   column, needed when data is a plain 3D array.
#'
#' @return an object of class \code{panel_causal_result} with
#'   \code{$effects} (per-feature estimate, standard error, statistic, p-value,
#'   Benjamini-Hochberg FDR, pre/post means and effect range),
#'   \code{$time_effects} (per-feature by per-time effect path) and, for
#'   synthetic control, \code{$weights} and \code{$weight_sum}.
#' @importFrom stats p.adjust pnorm var
#' @export
panel_causal_estimate <- \(formula,
                           data,
                           method = CAUSAL_DEFAULT_METHOD,
                           treatment_level = NULL,
                           control_level = NULL,
                           intervention_index = NULL,
                           scm_constraint = NULL,
                           scm_offset = NULL,
                           placebo = NULL,
                           bootstrap = NULL,
                           bootstrap_samples = NULL,
                           bootstrap_level = NULL,
                           pre_fit = NULL,
                           observation_meta = NULL,
                           time_meta = NULL) {
  parsed <- parse_causal_formula(formula)
  validate_causal_formula(parsed)
  method <- method %||% CAUSAL_DEFAULT_METHOD
  if (!(method %in% CAUSAL_METHODS)) {
    panel_abort(sprintf(
      "Unknown causal method '%s'; expected one of: %s.",
      method,
      paste(CAUSAL_METHODS, collapse = ", ")
    ))
  }

  panel <- as_causal_panel(data, observation_meta, time_meta)
  array_dims <- dim(panel$abundance_array)
  if (is.null(array_dims) || length(array_dims) != 3L) {
    panel_abort("The abundance data must be a three-dimensional array.")
  }
  dimension_names <- dimnames(panel$abundance_array)
  if (is.null(dimension_names) || any(vapply(dimension_names, is.null, logical(1L)))) {
    panel_abort("The abundance array must carry complete dimnames on all three dimensions.")
  }

  t0 <- resolve_t0(panel$time_meta, intervention_index, array_dims[3L])
  assignment <- resolve_treatment(panel$observation_meta,
                                  treatment_level, control_level)

  # The abundance array is [observation, feature, time]; estimators and the
  # Fortran kernels work on [observation, time, feature], so permute once here.
  value_otf <- aperm(panel$abundance_array, c(1L, 3L, 2L))

  panel_input <- list(
    value = value_otf,
    nobs = array_dims[1L],
    ntime = array_dims[3L],
    nfeat = array_dims[2L],
    t0 = t0,
    treat = assignment$treat,
    control = assignment$control,
    feature_ids = dimension_names[[2L]],
    time_ids = dimension_names[[3L]],
    observation_ids = dimension_names[[1L]]
  )

  # Named dispatch table: method name -> estimator. [[method]] replaces a switch.
  scm_constraint <- scm_constraint %||% CAUSAL_SCM_CONSTRAINT_DEFAULT
  scm_offset <- scm_offset %||% CAUSAL_SCM_OFFSET_DEFAULT
  placebo <- placebo %||% CAUSAL_PLACEBO_DEFAULT
  bootstrap <- bootstrap %||% CAUSAL_BOOTSTRAP_DEFAULT
  bootstrap_samples <- bootstrap_samples %||% CAUSAL_BOOTSTRAP_SAMPLES
  bootstrap_level <- bootstrap_level %||% CAUSAL_BOOTSTRAP_LEVEL
  pre_fit <- pre_fit %||% CAUSAL_PRE_FIT_DEFAULT
  if (!(scm_constraint %in% c(CAUSAL_SCM_NNLS, CAUSAL_SCM_SIMPLEX))) {
    panel_abort(sprintf(
      "Unknown scm_constraint '%s'; expected 'nnls' or 'simplex'.",
      scm_constraint
    ))
  }
  estimators <- list(
    did = \() did_estimator(panel_input),
    pooled_ols = \() pooled_ols_estimator(panel_input),
    event_study = \() event_study_estimator(panel_input),
    scm = \() scm_estimator(panel_input, scm_constraint, scm_offset)
  )
  result <- estimators[[method]]()

  # Benjamini-Hochberg FDR across features on the per-feature p-values.
  result$effects$fdr <- p.adjust(result$effects$p_value, method = "BH")
  result$effects <- result$effects[, CAUSAL_EFFECT_COLUMNS]

  # Inference depth (paper): placebo permutation, subject bootstrap CI and
  # Abadie pre-fit diagnostics, each optional and off by default. New columns
  # are appended by feature id so the canonical EFFECT_COLUMNS stay unchanged.
  placebo_out <- NULL
  if (isTRUE(placebo)) {
    if (sum(panel_input$control) >= 2L) {
      plc <- placebo_test(panel_input, scm_constraint, scm_offset)
      result$effects$placebo_p_value <- plc$placebo_p_value[match(result$effects$feature, plc$feature)]
      result$effects$placebo_fdr <- plc$placebo_fdr[match(result$effects$feature, plc$feature)]
      placebo_out <- plc
    } else {
      warning("placebo inference needs at least two control units; skipped.")
    }
  }
  bootstrap_out <- NULL
  if (isTRUE(bootstrap)) {
    bst <- bootstrap_ci(panel_input, method, scm_constraint, scm_offset,
                        bootstrap_samples, bootstrap_level)
    result$effects$bootstrap_lower <- bst$bootstrap_lower[match(result$effects$feature, bst$feature)]
    result$effects$bootstrap_upper <- bst$bootstrap_upper[match(result$effects$feature, bst$feature)]
    bootstrap_out <- bst
  }
  pre_fit_out <- NULL
  if (isTRUE(pre_fit) && method == CAUSAL_SCM) {
    pre_fit_out <- add_pre_fit(panel_input, result$weights, scm_offset)
  } else if (isTRUE(pre_fit)) {
    warning("pre_fit diagnostics are only defined for method = 'scm'; skipped.")
  }

  structure(
    list(
      call = match.call(),
      formula = formula,
      method = method,
      t0 = t0,
      t0_time_identifier = panel$time_meta$time_identifier[t0],
      treatment_level = assignment$treatment_level,
      control_level = assignment$control_level,
      effects = result$effects,
      time_effects = result$time_effects,
      weights = result$weights %||% NULL,
      weight_sum = result$weight_sum %||% NULL,
      placebo = placebo_out,
      bootstrap = bootstrap_out,
      pre_fit = pre_fit_out
    ),
    class = "panel_causal_result"
  )
}

#' Print a panel_causal_result
#'
#' @param x an object of class panel_causal_result.
#' @param ... unused, for compatibility.
#' @export
print.panel_causal_result <- \(x, ...) {
  cat("panelio causal estimate (method:", x$method, ")\n", sep = "")
  cat("intervention onset t0 =", x$t0, "(", x$t0_time_identifier, ")\n", sep = "")
  cat(sprintf("%d features, %.3f effective share below FDR 0.05\n",
              nrow(x$effects),
              if (nrow(x$effects) > 0L) {
                mean(x$effects$fdr < CAUSAL_TAIL_ALPHA, na.rm = TRUE)
              } else {
                NA_real_
              }))
  print(x$effects)
  invisible(x)
}
