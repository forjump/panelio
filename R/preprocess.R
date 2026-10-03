# =============================================================================
# R/preprocess.R
# Feature-axis preprocessing for the three-dimensional abundance array
# [observation, feature, time_identifier]. Filtering, zero replacement and
# normalization all operate on the feature (species) axis, dimension 2, and
# never change the observation or time axes. All operations are vectorized
# (no for/while loops) and every tunable value is a named constant.
# =============================================================================

# ---- local helpers ----------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

assert_preprocess_array <- \(abundance_array) {
  if (!is.array(abundance_array) || length(dim(abundance_array)) != 3L) {
    panel_abort(
      "abundance_array must be a three-dimensional array ",
      "[observation, feature, time_identifier]."
    )
  }

  dimension_names <- dimnames(abundance_array)

  if (is.null(dimension_names) || length(dimension_names) != 3L ||
        any(vapply(dimension_names, is.null, logical(1L)))) {
    panel_abort(
      "abundance_array must carry complete dimnames on all three dimensions ",
      "[observation, feature, time_identifier]."
    )
  }

  invisible(TRUE)
}

# Apply a per-sample function along the feature axis of every (observation,
# time) column and return an array with the exact original shape and dimnames.
# The apply+aperm pair is the single, dimension-safe entry point used by all
# zero-fill and normalization routines so the feature axis can never leak into
# the observation or time axes.
map_per_sample <- \(abundance_array, sample_fun) {
  dimension_names <- dimnames(abundance_array)
  per_sample <- apply(abundance_array, c(1L, 3L), sample_fun)
  # apply returns [feature, observation, time]; permute back to the canonical
  # [observation, feature, time] order.
  mapped <- aperm(per_sample, c(2L, 1L, 3L))
  dimnames(mapped) <- dimension_names
  mapped
}

# ---- filtering (feature axis) -----------------------------------------------
#' Filter features (species) of an abundance array
#'
#' Drops features from dimension 2 according to one or more prevalence /
#' abundance criteria. All criteria are evaluated per feature across every
#' (observation x time) sample, so the same features are dropped everywhere and
#' the panel stays balanced. Each criterion is optional: leave it NULL to skip
#' it. If every criterion is NULL the array is returned unchanged.
#'
#' @param abundance_array a three-dimensional array
#'   [observation, feature, time_identifier] with complete dimnames.
#' @param min_prevalence keep a feature present (abundance > min_abundance) in
#'   at least this fraction of samples; NULL disables the prevalence filter.
#' @param prevalence_threshold threshold used when min_prevalence is given.
#' @param min_abundance a feature is considered present in a sample when its
#'   abundance exceeds this value.
#' @param min_mean_abundance keep features whose mean abundance across samples
#'   is at least this value; NULL disables.
#' @param min_total_abundance keep features whose total abundance across
#'   samples is at least this value; NULL disables.
#' @param min_max_abundance keep features whose maximum abundance across
#'   samples is at least this value; this is the direct "drop features whose
#'   abundance never reaches this threshold" criterion; NULL disables.
#' @param max_zero_proportion keep features whose zero proportion across
#'   samples is at most this value; NULL disables.
#' @return the abundance array with a subset of the feature dimension.
#' @export
filter_features <- \(abundance_array,
                     min_prevalence = NULL,
                     prevalence_threshold = FILTER_PREVALENCE_THRESHOLD,
                     min_abundance = FILTER_MIN_ABUNDANCE,
                     min_mean_abundance = NULL,
                     min_total_abundance = NULL,
                     min_max_abundance = NULL,
                     max_zero_proportion = NULL) {
  assert_preprocess_array(abundance_array)
  feature_count <- dim(abundance_array)[2L]
  keep_mask <- rep(TRUE, feature_count)

  if (!is.null(min_prevalence)) {
    feature_prevalence <- apply(abundance_array > min_abundance, 2L, mean)
    keep_mask <- keep_mask & (feature_prevalence >= prevalence_threshold)
  }

  if (!is.null(min_mean_abundance)) {
    feature_mean <- apply(abundance_array, 2L, mean)
    keep_mask <- keep_mask & (feature_mean >= min_mean_abundance)
  }

  if (!is.null(min_total_abundance)) {
    feature_total <- apply(abundance_array, 2L, sum)
    keep_mask <- keep_mask & (feature_total >= min_total_abundance)
  }

  if (!is.null(min_max_abundance)) {
    feature_max <- apply(abundance_array, 2L, max)
    keep_mask <- keep_mask & (feature_max >= min_max_abundance)
  }

  if (!is.null(max_zero_proportion)) {
    zero_proportion <- 1 - apply(abundance_array > 0, 2L, mean)
    keep_mask <- keep_mask & (zero_proportion <= max_zero_proportion)
  }

  abundance_array[, keep_mask, , drop = FALSE]
}

# ---- zero replacement (feature axis) ----------------------------------------
#' Replace zeros of an abundance array
#'
#' Replaces zero abundances along the feature axis of every (observation, time)
#' sample. This is required before log based normalization (CLR, RLE) because
#' the logarithm of zero is undefined. The array shape and dimnames are
#' preserved.
#'
#' @param abundance_array a three-dimensional array
#'   [observation, feature, time_identifier].
#' @param method one of "none", "constant", "multiplicative", "half_min",
#'   "sqrt_min", "bayesian"; NULL selects the default method. "none" performs no
#'   operation. Every other method touches only the zero entries (the
#'   "multiplicative" and "bayesian" methods also rescale the positive entries
#'   so the sample composition sum is preserved); no entry is ever replaced by
#'   zero.
#' @param constant pseudo-count for the "constant" method; NULL uses the
#'   default.
#' @param delta replacement value for the "multiplicative" method; NULL uses
#'   the default.
#' @param fraction replacement fraction for the "half_min" and "sqrt_min"
#'   methods; NULL uses the default (0.5).
#' @return the abundance array with zeros replaced; same shape and dimnames.
#' @export
zero_fill_features <- \(abundance_array,
                        method = NULL,
                        constant = NULL,
                        delta = NULL,
                        fraction = ZERO_FILL_FRACTION) {
  assert_preprocess_array(abundance_array)
  method <- method %||% ZERO_FILL_DEFAULT_METHOD

  if (!(method %in% ZERO_FILL_METHODS)) {
    panel_abort(sprintf(
      "Unknown zero-fill method '%s'; expected one of: %s.",
      method,
      paste(ZERO_FILL_METHODS, collapse = ", ")
    ))
  }

  # Named dispatch table: method name -> per-sample handler. Selecting with
  # [[method]] replaces a switch so the body carries no branch on the method.
  zero_fill_handlers <- list(
    none = identity,
    constant = \(values) {
      values[values == 0] <- constant %||% ZERO_FILL_CONSTANT
      values
    },
    half_min = \(values) {
      positive_values <- values[values > 0]
      if (length(positive_values) == 0L) return(values)
      values[values == 0] <- fraction * min(positive_values)
      values
    },
    sqrt_min = \(values) {
      positive_values <- values[values > 0]
      if (length(positive_values) == 0L) return(values)
      values[values == 0] <- fraction * sqrt(min(positive_values))
      values
    },
    multiplicative = \(values) {
      zero_count <- sum(values == 0)
      if (zero_count == 0L) return(values)
      replacement <- delta %||% ZERO_FILL_DELTA
      rescaled <- values
      rescaled[values == 0] <- replacement
      positive_sum <- sum(values[values > 0])
      scale <- (1 - zero_count * replacement) / positive_sum
      rescaled[values > 0] <- values[values > 0] * scale
      rescaled
    },
    bayesian = \(values) {
      zero_count <- sum(values == 0)
      if (zero_count == 0L) return(values)
      if (zero_count == length(values)) return(values)  # all-zero sample
      sample_sum <- sum(values)
      positive_values <- values[values > 0]
      # Perks/Dirichlet-prior estimate of the missing mass: the geometric mean
      # of the positive proportions is the non-zero share; its complement delta
      # is split equally among the zeros and the positives are rescaled by
      # (1 - delta) so the composition sum is preserved.
      positive_share <- exp(mean(log(positive_values / sample_sum)))
      delta_share <- 1 - positive_share
      rescaled <- values
      rescaled[values == 0] <- sample_sum * delta_share / zero_count
      rescaled[values > 0] <- positive_values * (1 - delta_share)
      rescaled
    }
  )

  map_per_sample(abundance_array, zero_fill_handlers[[method]])
}

# ---- normalization (feature axis) -------------------------------------------
# Geometric-mean reference used by the RLE method: one positive geometric mean
# per feature across all samples.
rle_reference <- \(abundance_array) {
  log_values <- log(abundance_array)
  log_values[!is.finite(log_values)] <- NA_real_
  reference <- exp(colMeans(log_values, na.rm = TRUE))
  reference[!is.finite(reference)] <- NA_real_
  reference
}

# Cumulative-sum scaling factor per sample (metagenomeSeq style).
css_scale <- \(values, quantile_value) {
  positive_values <- values[values > 0]
  if (length(positive_values) == 0L) return(NA_real_)
  ordered_values <- sort(positive_values)
  index_count <- max(1L, floor(quantile_value * length(ordered_values)))
  sum(ordered_values[seq_len(index_count)])
}

# Full-rank orthonormal basis used by the ILR transform: the D x D matrix whose
# first D-1 columns are the (column-normalized) Helmert contrasts and whose last
# column is the constant vector 1/sqrt(D). contr.helmert(D) already returns the
# D x (D-1) contrast columns (orthogonal to the constant and to each other), so
# only normalization is required before appending the constant column. Applying
# the basis keeps the feature dimension D (an isometry of the CLR space) while
# remaining an equidistant log-ratio map.
ilr_basis <- \(dimension) {
  helmert_columns <- stats::contr.helmert(dimension)
  basis <- apply(helmert_columns, 2L, \(column) column / sqrt(sum(column^2)))
  cbind(basis, rep(1 / sqrt(dimension), dimension))
}

#' Normalize an abundance array
#'
#' Applies a normalization to every (observation, time) sample along the
#' feature axis. Zeros are replaced first whenever the selected method requires
#' positive values; pass zero_method = "none" to disable that automatic step.
#'
#' @param abundance_array a three-dimensional array
#'   [observation, feature, time_identifier].
#' @param method one of "none", "clr", "tss", "rle", "css", "log1p", "ilr",
#'   "rclr", "z_score", "min_max"; NULL selects the default (CLR). "ilr" uses a
#'   full-rank Helmert orthonormal basis so the feature dimension is preserved;
#'   "rclr", "z_score" and "min_max" handle zeros internally and never trigger
#'   the automatic zero-fill.
#' @param zero_method zero replacement method applied first for methods that
#'   need positive values; NULL uses the default.
#' @param constant pseudo-count passed to the zero-fill step.
#' @param delta replacement value passed to the zero-fill step.
#' @param fraction replacement fraction passed to the zero-fill step.
#' @param reference optional per-feature reference for the "rle" method; NULL
#'   computes the geometric mean reference automatically.
#' @param quantile quantile for the "css" method; NULL uses the default.
#' @return the normalized abundance array; same shape and dimnames.
#' @export
normalize_features <- \(abundance_array,
                        method = NORMALIZATION_DEFAULT_METHOD,
                        zero_method = NULL,
                        constant = NULL,
                        delta = NULL,
                        fraction = ZERO_FILL_FRACTION,
                        reference = NULL,
                        quantile = NULL) {
  assert_preprocess_array(abundance_array)
  method <- method %||% NORMALIZATION_DEFAULT_METHOD

  if (!(method %in% NORMALIZE_METHODS)) {
    panel_abort(sprintf(
      "Unknown normalization method '%s'; expected one of: %s.",
      method,
      paste(NORMALIZE_METHODS, collapse = ", ")
    ))
  }

  # Auto zero-fill only for methods that need positive values, and only unless
  # the caller opts out via zero_method = "none". Data-driven flags replace the
  # previous fill-method branch.
  requires_positive <- method %in% NORMALIZE_REQUIRES_POSITIVE
  zero_fill_active <- !identical(zero_method %||% ZERO_FILL_DEFAULT_METHOD, ZERO_NONE)
  if (requires_positive && any(abundance_array == 0) && zero_fill_active) {
    abundance_array <- zero_fill_features(
      abundance_array,
      method = zero_method,
      constant = constant,
      delta = delta,
      fraction = fraction
    )
  }

  # Named dispatch table: method name -> per-sample handler. A single [[method]]
  # lookup replaces the sequential per-method branches. The RLE reference is
  # computed once across all samples (a per-feature vector shared by every
  # sample), the CSS quantile is fixed before the per-sample pass, and the ILR
  # orthonormal basis is built once for the feature dimension.
  quantile_value <- quantile %||% CSS_DEFAULT_QUANTILE
  feature_reference <- reference %||% rle_reference(abundance_array)
  ilr_basis_value <- ilr_basis(dim(abundance_array)[2L])
  normalize_handlers <- list(
    none = identity,
    clr = \(values) log(values) - mean(log(values)),
    tss = \(values) {
      sample_sum <- sum(values)
      safe_sum <- sample_sum + as.integer(sample_sum == 0) # never 0
      values / safe_sum
    },
    log1p = log1p,
    css = \(values) {
      scale_factor <- css_scale(values, quantile_value)
      scale_factor[is.na(scale_factor) | scale_factor == 0] <- 1
      values / scale_factor
    },
    rle = \(values) {
      ratio <- ifelse(values > 0 & feature_reference > 0,
                      values / feature_reference, NA_real_)
      scale_factor <- stats::median(ratio, na.rm = TRUE)
      scale_factor[is.na(scale_factor) | scale_factor == 0] <- 1
      values / scale_factor
    },
    ilr = \(values) {
      centered <- log(values) - mean(log(values))  # clr, sums to zero
      as.numeric(centered %*% ilr_basis_value)      # equidistant log-ratio map
    },
    rclr = \(values) {
      positive_mask <- values > 0
      if (!any(positive_mask)) return(values)
      transformed <- values
      transformed[positive_mask] <-
        log(values[positive_mask]) - mean(log(values[positive_mask]))
      transformed[!positive_mask] <- 0
      transformed
    },
    z_score = \(values) {
      sample_mean <- mean(values)
      sample_sd <- stats::sd(values)
      safe_sd <- sample_sd + as.integer(sample_sd == 0) # never 0
      (values - sample_mean) / safe_sd
    },
    min_max = \(values) {
      sample_min <- min(values)
      sample_range <- max(values) - sample_min
      safe_range <- sample_range + as.integer(sample_range == 0) # never 0
      (values - sample_min) / safe_range
    }
  )

  map_per_sample(abundance_array, normalize_handlers[[method]])
}

# ---- combined preprocessing --------------------------------------------------
#' Preprocess an abundance array
#'
#' Pipeline wrapper: optional feature filtering, zero replacement and
#' normalization, applied in that order along the feature axis. The original
#' observation and time axes are never modified.
#'
#' @param abundance_array a three-dimensional array
#'   [observation, feature, time_identifier].
#' @param filter NULL or "none" (skip filtering) or a named list of arguments
#'   passed to filter_features().
#' @param zero NULL (auto zero-fill only when the normalization needs it) or a
#'   named list of arguments passed to zero_fill_features().
#' @param normalize NULL or "clr" (default), a method name, or a named list of
#'   arguments passed to normalize_features().
#' @return a list with the preprocessed abundance_array, the kept/removed
#'   feature names, and the applied filter/zero/normalize spec for recording.
#' @export
preprocess_array <- \(abundance_array,
                      filter = NULL,
                      zero = NULL,
                      normalize = "clr") {
  assert_preprocess_array(abundance_array)
  original_features <- dimnames(abundance_array)[[2L]]

  skip_filter <- is.null(filter) || identical(filter, "none")
  filtered <- if (skip_filter) {
    abundance_array
  } else {
    do.call(filter_features, c(list(abundance_array = abundance_array), filter))
  }

  kept_features <- dimnames(filtered)[[2L]]
  removed_features <- setdiff(original_features, kept_features)

  zeroed <- if (is.list(zero)) {
    do.call(zero_fill_features, c(list(abundance_array = filtered), zero))
  } else {
    filtered
  }

  normalized <- if (is.list(normalize)) {
    do.call(normalize_features, c(list(abundance_array = zeroed), normalize))
  } else {
    normalize_features(zeroed, method = normalize %||% NORMALIZATION_DEFAULT_METHOD)
  }

  list(
    abundance_array = normalized,
    kept_features = kept_features,
    removed_features = removed_features,
    filter = filter,
    zero = zero,
    normalize = normalize
  )
}

#' Preprocess a panel_array object
#'
#' Applies the feature-axis preprocessing pipeline to the abundance array of a
#' panel_array and returns a new panel_array. Because only the feature axis is
#' touched, observation_meta and time_meta are carried over unchanged.
#'
#' @param panel_array a panel_array object (as returned by as_panel_array()).
#' @param filter NULL or "none" (skip filtering) or a named list of
#'   filter_features() arguments.
#' @param zero NULL or a named list of zero_fill_features() arguments.
#' @param normalize NULL, a method name, or a named list of
#'   normalize_features() arguments.
#' @return a new panel_array whose abundance array has been preprocessed.
#' @export
preprocess_panel <- \(panel_array,
                      filter = NULL,
                      zero = NULL,
                      normalize = "clr") {
  if (!inherits(panel_array, "panel_array")) {
    panel_abort("panel_array must be a panel_array object, as returned by as_panel_array().")
  }

  result <- preprocess_array(panel_array$abundance_array,
                             filter = filter,
                             zero = zero,
                             normalize = normalize)

  structure(
    list(
      abundance_array = result$abundance_array,
      observation_meta = panel_array$observation_meta,
      time_meta = panel_array$time_meta
    ),
    class = "panel_array"
  )
}
