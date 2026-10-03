# =============================================================================
# R/smooth.R
# Time-axis smoothing for the three-dimensional abundance array
# [observation, feature, time_identifier]. Smoothing operates on the time axis,
# dimension 3, and never changes the observation or feature axes. Every method
# runs its iteration inside the Fortran kernels under src/panelio_smoothing.f90;
# the R layer only dispatches a single .Fortran call per method (no for/while
# loops, purely vectorized). All tunables are named constants and every method
# parameter defaults to NULL (the named default applies).
#
# Breakpoint semantics: the Kalman family (EKF / RTS) smooths the whole series
# and intentionally ignores breakpoints; the benchmark methods (gaussian, loess,
# savgol, moving average, ewma, spline) honour breakpoints so an intervention
# effect at the breakpoint is never smoothed away. A breakpoint marks the first
# time index of a new post-intervention segment.
# =============================================================================

# Resolve the 1-based breakpoint indices for a smoothing call. Kalman methods
# ignore breakpoints entirely. Explicit breakpoints are validated and used
# as-is; otherwise, when time_meta carries the intervention flag, the first
# intervention timepoint becomes the single breakpoint.
resolve_breakpoints <- \(breakpoints, time_meta, ntime, method) {
  if (method %in% SMOOTH_KALMAN_METHODS) {
    return(integer(0L))
  }

  if (!is.null(breakpoints)) {
    breakpoint_indices <- as.integer(breakpoints)
    if (any(breakpoint_indices < 1L | breakpoint_indices > ntime)) {
      panel_abort(
        "breakpoints must be 1-based time indices within the time axis ",
        "[1, ntime]."
      )
    }
    return(sort(unique(breakpoint_indices)))
  }

  if (!is.null(time_meta) &&
        "is_intervention_timepoint" %in% names(time_meta)) {
    intervention_indices <- which(as.logical(time_meta$is_intervention_timepoint))
    if (length(intervention_indices) > 0L) {
      # Once an intervention starts it persists; the single meaningful break
      # is the onset (the earliest intervention timepoint).
      return(as.integer(min(intervention_indices)))
    }
  }

  integer(0L)
}

# Derive the per (observation, time) sequencing depth as the row sum of the
# feature axis (dimension 2). Only used by the negative-binomial observation
# model, which expects raw counts.
derive_library_size <- \(abundance_array) {
  apply(abundance_array, c(1L, 3L), sum)
}

# Derive the sampling interval vector (length ntime) used to scale the Kalman
# process noise for irregular schedules. Prefers the delta_t column of the time
# metadata, then the elapsed_time differences, then unit intervals.
derive_sampling_intervals <- \(time_meta, ntime) {
  if (!is.null(time_meta) && "delta_t" %in% names(time_meta)) {
    return(as.numeric(time_meta$delta_t))
  }
  if (!is.null(time_meta) && "elapsed_time" %in% names(time_meta)) {
    elapsed <- as.numeric(time_meta$elapsed_time)
    return(c(0, pmax(diff(elapsed), 0)))
  }
  rep(1, ntime)
}

# ---- smoothing (time axis) ---------------------------------------------------
#' Smooth an abundance array along its time axis
#'
#' Applies one smoothing method along the time axis (dimension 3) of an
#' abundance array \code{[observation, feature, time_identifier]}. Every
#' (observation, feature) series is smoothed independently inside the Fortran
#' kernel; the R layer makes a single \code{.Fortran} dispatch and performs no
#' explicit loops.
#'
#' @param abundance_array a three-dimensional array
#'   [observation, feature, time_identifier] with complete dimnames.
#' @param method one of \code{SMOOTH_METHODS}: \code{"ekf"} (default, extended /
#'   scalar Kalman forward filter), \code{"rts"} (Rauch-Tung-Striebel backward
#'   smoother), \code{"gaussian"}, \code{"loess"}, \code{"savgol"},
#'   \code{"moving_average"}, \code{"ewma"}, \code{"spline"}.
#' @param observation_model for the Kalman family, \code{"gaussian"} (default,
#'   scalar local level) or \code{"nb"} (paper-depth negative-binomial
#'   observation over raw counts: latent log relative abundance, expected count
#'   = library size times softmax, delta-method measurement variance).
#' @param hierarchical logical; when the observation model is \code{"nb"},
#'   estimate the per-feature innovation variance and dispersion shared across
#'   observations by a hierarchical EM loop.
#' @param library_size per (observation, time) sequencing depth as an
#'   \code{nobs x ntime} matrix; NULL derives it as the row sum of the feature
#'   axis (only relevant for \code{observation_model = "nb"}).
#' @param dispersion optional per-feature negative-binomial dispersion (scalar
#'   or length-nfeature); NULL uses the named default (and the EM when
#'   \code{hierarchical = TRUE}).
#' @param em_iterations EM iterations for the hierarchical estimator (paper uses
#'   ~80); only used when \code{hierarchical = TRUE}.
#' @param time_scaling logical; when TRUE the Kalman process noise is scaled by
#'   the sampling interval (from \code{time_meta}) for irregular schedules.
#' @param breakpoints optional 1-based time indices marking the first time of
#'   each post-intervention segment. Ignored by \code{"ekf"} and \code{"rts"};
#'   passed to the Fortran kernel for the other methods so an intervention
#'   effect is not smoothed away.
#' @param time_meta optional time metadata with an
#'   \code{is_intervention_timepoint} column; used to derive the single
#'   breakpoint (the first intervention time) when \code{breakpoints} is NULL,
#'   and the sampling intervals when \code{time_scaling = TRUE}.
#' @param qratio Q/R ratio of the Kalman family (process vs observation noise).
#' @param initvar scale of the initial state variance (P0 = initvar * R).
#' @param sigma Gaussian kernel bandwidth in time-index units.
#' @param span tricube loess span (fraction of a segment).
#' @param halfwindow half-window in time indices for moving average and savgol.
#' @param degree polynomial degree for savgol.
#' @param alpha EWMA smoothing factor in (0, 1).
#' @param lambda second-difference penalty of the spline smoother.
#'
#' @return an array with the same shape and dimnames as \code{abundance_array},
#'   smoothed along the time axis.
#' @export
smooth_array <- \(abundance_array,
                  method = SMOOTH_DEFAULT_METHOD,
                  observation_model = NULL,
                  hierarchical = NULL,
                  library_size = NULL,
                  dispersion = NULL,
                  em_iterations = NULL,
                  time_scaling = NULL,
                  breakpoints = NULL,
                  time_meta = NULL,
                  qratio = NULL,
                  initvar = NULL,
                  sigma = NULL,
                  span = NULL,
                  halfwindow = NULL,
                  degree = NULL,
                  alpha = NULL,
                  lambda = NULL) {
  assert_preprocess_array(abundance_array)
  dimension_names <- dimnames(abundance_array)
  array_dims <- dim(abundance_array)

  method <- method %||% SMOOTH_DEFAULT_METHOD
  if (!(method %in% SMOOTH_METHODS)) {
    panel_abort(sprintf(
      "Unknown smoothing method '%s'; expected one of: %s.",
      method,
      paste(SMOOTH_METHODS, collapse = ", ")
    ))
  }

  observation_model <- observation_model %||% SMOOTH_OBSERVATION_DEFAULT
  if (!(observation_model %in% c(SMOOTH_OBSERVATION_GAUSSIAN, SMOOTH_OBSERVATION_NB))) {
    panel_abort(sprintf(
      "Unknown observation_model '%s'; expected 'gaussian' or 'nb'.",
      observation_model
    ))
  }
  if (observation_model == SMOOTH_OBSERVATION_NB &&
        !(method %in% SMOOTH_KALMAN_METHODS)) {
    panel_abort(
      "observation_model = 'nb' applies only to the Kalman family (ekf / rts)."
    )
  }

  breakpoint_indices <- resolve_breakpoints(
    breakpoints, time_meta, array_dims[3L], method
  )
  nbreakpoints <- length(breakpoint_indices)
  # Fortran assumed-size arrays need at least one element; a length-one vector
  # with nbreakpoints == 0 is never dereferenced by the kernel.
  breakpoint_vector <- if (nbreakpoints > 0L) breakpoint_indices else integer(1L)

  smoothed <- array(0, dim = array_dims)

  # Per-kernel call helpers: each performs exactly one .Fortran dispatch along
  # the time axis and returns the smoothed array. Arguments are named so the
  # returned list exposes the output slot as $y.
  kalman_call <- \(smooth_mode, qratio_value, initvar_value) {
    .Fortran(
      "smooth_kalman",
      x = abundance_array, nobs = array_dims[1L], nfeat = array_dims[2L],
      ntime = array_dims[3L],
      qratio = as.numeric(qratio_value), initvar = as.numeric(initvar_value),
      mode = as.integer(smooth_mode), y = smoothed
    )$y
  }
  nb_kalman_call <- \(smooth_mode) {
    lib_matrix <- library_size %||% derive_library_size(abundance_array)
    if (!is.matrix(lib_matrix) ||
          !identical(dim(lib_matrix), c(array_dims[1L], array_dims[3L]))) {
      panel_abort(sprintf(
        "library_size must be a %d x %d matrix (observation x time).",
        array_dims[1L], array_dims[3L]
      ))
    }
    nfeat <- array_dims[2L]
    phi_default <- if (is.null(dispersion)) {
      rep(SMOOTH_NB_DEFAULT_DISPERSION, nfeat)
    } else if (length(dispersion) == 1L) {
      rep(as.numeric(dispersion), nfeat)
    } else if (length(dispersion) == nfeat) {
      as.numeric(dispersion)
    } else {
      panel_abort(sprintf(
        "dispersion must be scalar or length %d (nfeature).", nfeat
      ))
    }
    q_default <- rep(as.numeric(qratio %||% SMOOTH_KALMAN_Q_RATIO), nfeat)
    use_dt <- as.integer(isTRUE(time_scaling))
    dt <- if (use_dt == 1L) {
      derive_sampling_intervals(time_meta, array_dims[3L])
    } else {
      rep(1, array_dims[3L])
    }
    n_em <- if (isTRUE(hierarchical)) {
      as.integer(em_iterations %||% SMOOTH_NB_EM_ITERATIONS)
    } else {
      0L
    }
    .Fortran(
      "smooth_kalman_nb",
      x = abundance_array, nobs = array_dims[1L], nfeat = nfeat,
      ntime = array_dims[3L],
      lib = as.matrix(lib_matrix),
      phi = phi_default, q = q_default,
      mode = as.integer(smooth_mode),
      n_em = n_em,
      dt = as.numeric(dt),
      use_dt = use_dt,
      y = smoothed
    )$y
  }
  locreg_call <- \(kernel_code, poly_degree, band_value) {
    .Fortran(
      "smooth_locreg",
      x = abundance_array, nobs = array_dims[1L], nfeat = array_dims[2L],
      ntime = array_dims[3L],
      kcode = as.integer(kernel_code), degree = as.integer(poly_degree),
      band = as.numeric(band_value),
      nbrk = as.integer(nbreakpoints), brk = as.integer(breakpoint_vector),
      y = smoothed
    )$y
  }
  ewma_call <- \(alpha_value) {
    .Fortran(
      "smooth_ewma",
      x = abundance_array, nobs = array_dims[1L], nfeat = array_dims[2L],
      ntime = array_dims[3L],
      alpha = as.numeric(alpha_value),
      nbrk = as.integer(nbreakpoints), brk = as.integer(breakpoint_vector),
      y = smoothed
    )$y
  }
  spline_call <- \(lambda_value) {
    .Fortran(
      "smooth_spline",
      x = abundance_array, nobs = array_dims[1L], nfeat = array_dims[2L],
      ntime = array_dims[3L],
      lambda = as.numeric(lambda_value),
      nbrk = as.integer(nbreakpoints), brk = as.integer(breakpoint_vector),
      y = smoothed
    )$y
  }

  # Named dispatch table: method name -> one Fortran kernel call. Selecting with
  # [[method]] replaces a switch so the R layer carries no branch on the method.
  smoothing_handlers <- list(
    ekf = \() {
      if (observation_model == SMOOTH_OBSERVATION_NB) {
        nb_kalman_call(0L)
      } else {
        kalman_call(
          0L,
          qratio %||% SMOOTH_KALMAN_Q_RATIO,
          initvar %||% SMOOTH_KALMAN_INITVAR_SCALE
        )
      }
    },
    rts = \() {
      if (observation_model == SMOOTH_OBSERVATION_NB) {
        nb_kalman_call(1L)
      } else {
        kalman_call(
          1L,
          qratio %||% SMOOTH_KALMAN_Q_RATIO,
          initvar %||% SMOOTH_KALMAN_INITVAR_SCALE
        )
      }
    },
    gaussian = \() locreg_call(1L, 0L, sigma %||% SMOOTH_GAUSSIAN_SIGMA),
    loess = \() locreg_call(2L, SMOOTH_LOESS_DEGREE, span %||% SMOOTH_LOESS_SPAN),
    savgol = \() {
      locreg_call(
        0L, degree %||% SMOOTH_SAVGOL_DEGREE,
        halfwindow %||% SMOOTH_SAVGOL_HALFWINDOW
      )
    },
    moving_average = \() locreg_call(0L, 0L, halfwindow %||% SMOOTH_MA_HALFWINDOW),
    ewma = \() ewma_call(alpha %||% SMOOTH_EWMA_ALPHA),
    spline = \() spline_call(lambda %||% SMOOTH_SPLINE_LAMBDA)
  )

  smoothed <- smoothing_handlers[[method]]()
  dimnames(smoothed) <- dimension_names
  smoothed
}

#' Smooth a panel array along its time axis
#'
#' Applies \code{smooth_array} to the abundance array of a panel_array object
#' and returns a new panel_array. When \code{breakpoints} is NULL the first
#' intervention timepoint in the panel's time metadata is used as the breakpoint
#' (except for the Kalman family, which ignores breakpoints).
#'
#' @param panel_array a panel_array object, as returned by as_panel_array().
#' @param ... arguments passed to smooth_array().
#'
#' @return a panel_array object whose abundance_array has been smoothed along
#'   the time axis; the observation_meta and time_meta are carried over.
#' @export
smooth_panel <- \(panel_array, ...) {
  if (!inherits(panel_array, "panel_array")) {
    panel_abort("panel_array must be a panel_array object, as returned by as_panel_array().")
  }

  smoothed_abundance <- smooth_array(
    panel_array$abundance_array,
    time_meta = panel_array$time_meta,
    ...
  )

  structure(
    list(
      abundance_array = smoothed_abundance,
      observation_meta = panel_array$observation_meta,
      time_meta = panel_array$time_meta
    ),
    class = "panel_array"
  )
}
