# =============================================================================
# R/convert_panel.R
# Bidirectional object model conversions.
#
#   long table (data.frame)  --as_panel_array()-->  panel_array
#   panel_array components   --as_panel_table()-->  long table
#
#   long table / manifest raw text lines --as_observation_meta()--> meta table
#   long table / manifest raw text lines --as_time_meta()--------> meta table
# =============================================================================

new_panel_array <- \(abundance_array, observation_meta, time_meta) {
  structure(
    list(
      abundance_array = abundance_array,
      observation_meta = observation_meta,
      time_meta = time_meta
    ),
    class = "panel_array"
  )
}

# ---- as_observation_meta ----------------------------------------------------
#' Extract the observation metadata
#'
#' Reads the observation metadata from a long-format table, from a `panel_array`
#' object, or from raw manifest text lines. The returned table uses the canonical
#' manifest column names `observation_id` and `label`.
#'
#' @param x A long-format data.frame, a `panel_array` object, or a character
#'   vector of raw manifest text lines.
#' @param ... Ignored, kept for S3 method dispatch.
#'
#' @return A data.frame of observation metadata with the columns
#'   `observation_id` and `label`.
#' @export
as_observation_meta <- \(x, ...) {
  UseMethod("as_observation_meta")
}

#' @rdname as_observation_meta
#' @export
as_observation_meta.default <- \(x, ...) {
  panel_abort(
    "as_observation_meta accepts only a long-format data.frame, a panel_array ",
    "object, or manifest text lines."
  )
}

#' @rdname as_observation_meta
#' @export
as_observation_meta.data.frame <- \(x, ...) {
  assert_long_table(x, "as_observation_meta")
  first_flag <- !duplicated(x[[OBSERVATION_COLUMN]])
  observation_meta <- data.frame(
    observation_id = as.character(x[[OBSERVATION_COLUMN]][first_flag]),
    label = x[[LABEL_COLUMN]][first_flag],
    stringsAsFactors = FALSE
  )
  rownames(observation_meta) <- NULL
  observation_meta
}
#' @rdname as_observation_meta
#' @export
as_observation_meta.character <- \(x, ...) {
  block_lines <- parse_manifest_blocks(x)
  block_name <- OBSERVATION_META_BLOCK

  if (!block_name %in% names(block_lines)) {
    panel_abort(sprintf("The manifest text is missing the block [%s].", block_name))
  }

  observation_meta <- coerce_observation_meta_frame(
    parse_manifest_pairs(block_lines[[block_name]])
  )
  assert_observation_meta(observation_meta)
  observation_meta[, OBSERVATION_META_COLUMNS, drop = FALSE]
}

#' @rdname as_observation_meta
#' @export
as_observation_meta.panel_array <- \(x, ...) {
  assert_observation_meta(x$observation_meta)
  x$observation_meta
}

# ---- as_time_meta -----------------------------------------------------------
#' Extract the time layer metadata
#'
#' Reads the time layer metadata from a long-format table, from a `panel_array`
#' object, or from raw manifest text lines. The returned table uses the
#' canonical manifest column names `time_identifier`, `elapsed_time` and
#' `is_intervention_timepoint`.
#'
#' @param x A long-format table, a `panel_array` object, or a character vector
#'   of raw manifest text lines.
#' @param ... Ignored, kept for S3 method dispatch.
#'
#' @return A data.frame of time layer metadata.
#' @export
as_time_meta <- \(x, ...) {
  UseMethod("as_time_meta")
}

#' @rdname as_time_meta
#' @export
as_time_meta.default <- \(x, ...) {
  panel_abort(
    "as_time_meta accepts only a long-format data.frame, a panel_array ",
    "object, or manifest text lines."
  )
}

#' @rdname as_time_meta
#' @export
as_time_meta.data.frame <- \(x, ...) {
  assert_long_table(x, "as_time_meta")
  first_flag <- !duplicated(x[[TIME_IDENTIFIER_COLUMN]])
  elapsed_time <- as.numeric(x[[ELAPSED_TIME_COLUMN]][first_flag])
  time_meta <- data.frame(
    time_identifier = as.character(x[[TIME_IDENTIFIER_COLUMN]][first_flag]),
    elapsed_time = elapsed_time,
    # delta_t is the interval since the previous time layer. The long table does
    # not carry it, so it is derived from elapsed_time, which is unique by
    # construction. The earliest layer has no predecessor and gets 0.
    delta_t = derive_delta_t(elapsed_time),
    is_intervention_timepoint = as.logical(x[[IS_INTERVENTION_COLUMN]][first_flag]),
    stringsAsFactors = FALSE
  )
  rownames(time_meta) <- NULL
  time_meta
}

#' @rdname as_time_meta
#' @export
as_time_meta.character <- \(x, ...) {
  block_lines <- parse_manifest_blocks(x)
  block_name <- TIME_META_BLOCK

  if (!block_name %in% names(block_lines)) {
    panel_abort(sprintf("The manifest text is missing the block [%s].", block_name))
  }

  time_meta <- parse_manifest_pairs(block_lines[[block_name]])
  assert_time_meta(time_meta)
  time_meta[, TIME_META_COLUMNS, drop = FALSE]
}
#' @rdname as_time_meta
#' @export
as_time_meta.panel_array <- \(x, ...) {
  assert_time_meta(x$time_meta)
  x$time_meta
}

# ---- as_panel_array ---------------------------------------------------------
#' Convert a long-format table into a `panel_array` object
#'
#' Validates that the input is a balanced longitudinal panel, splits out the
#' three dimensional abundance array together with the observation metadata and
#' the time layer metadata, and returns a complete `panel_array` S3 object. An
#' input that is already a `panel_array` is only validated and returned as is.
#'
#' @param x A long-format data.frame or a `panel_array` object.
#'
#' @return A `panel_array` object holding `abundance_array`, `observation_meta`
#'   and `time_meta`.
#' @export
as_panel_array <- \(x) {
  if (inherits(x, "panel_array")) {
    assert_array_meta_alignment(
      x$abundance_array,
      x$observation_meta,
      x$time_meta
    )
    return(x)
  }

  assert_long_table(x, "as_panel_array")
  assert_unique_panel_keys(x)
  observation_levels <- unique(as.character(x[[OBSERVATION_COLUMN]]))
  time_levels <- unique(as.character(x[[TIME_IDENTIFIER_COLUMN]]))
  feature_levels <- unique(as.character(x[[FEATURE_COLUMN]]))
  assert_balanced_panel(x, observation_levels, time_levels, feature_levels)
  abundance_array <- long_to_panel_array(x, observation_levels, time_levels, feature_levels)
  observation_meta <- as_observation_meta(x)
  time_meta <- as_time_meta(x)
  assert_array_meta_alignment(abundance_array, observation_meta, time_meta)
  new_panel_array(abundance_array, observation_meta, time_meta)
}

# ---- as_panel_table ---------------------------------------------------------
#' Assemble or read a long-format table
#'
#' A multi-mode entry point supporting two calling conventions.
#'
#' Mode one: pass the archive file path, read the archive and emit the
#' long-format table.
#'
#' Mode two: pass `x = NULL` together with `abundance_array`,
#' `observation_meta` and `time_meta`; the array is expanded and the two
#' metadata tables are merged back into a long-format table.
#'
#' @param x Path to a tar.gz archive as a single string, or `NULL`.
#' @param abundance_array Optional three dimensional base R array with
#'   dimensions `[observation, feature, time_identifier]`.
#' @param observation_meta Optional observation metadata data.frame.
#' @param time_meta Optional time layer metadata data.frame.
#'
#' @return A long-format data.frame with the fixed column order `observation`,
#'   `time_identifier`, `elapsed_time`, `feature`, `abundance_value`, `label`,
#'   `is_intervention`.
#' @export
as_panel_table <- \(x, abundance_array = NULL, observation_meta = NULL, time_meta = NULL) {
  if (is.character(x)) {
    assert_scalar_text(x, "x")
    return(read_panel_table(x))
  }

  if (!is.null(x)) {
    panel_abort("The x argument of as_panel_table must be an archive path or NULL.")
  }

  missing_arguments <- c(
    if (is.null(abundance_array)) "abundance_array",
    if (is.null(observation_meta)) "observation_meta",
    if (is.null(time_meta)) "time_meta"
  )

  if (length(missing_arguments) > 0L) {
    panel_abort(sprintf(
      paste0(
        "When x is NULL, abundance_array, observation_meta and time_meta must ",
        "all be supplied; missing: %s."
      ),
      paste(missing_arguments, collapse = ", ")
    ))
  }

  assert_array_meta_alignment(abundance_array, observation_meta, time_meta)
  long_table <- panel_array_to_long(abundance_array)
  attach_metadata(long_table, observation_meta, time_meta)
}

# ---- panel_array methods ----------------------------------------------------
#' @param x A `panel_array` object.
#' @param ... Ignored.
#' @rdname as_panel_array
#' @export
print.panel_array <- \(x, ...) {
  cat("<panel_array>\n")
  cat("  abundance_array : dim = ", paste(dim(x$abundance_array), collapse = " x "), "\n", sep = "")
  cat("  observation_meta: ", nrow(x$observation_meta), " rows\n", sep = "")
  cat("  time_meta       : ", nrow(x$time_meta), " rows\n", sep = "")
  invisible(x)
}

#' @rdname as_panel_array
#' @export
dim.panel_array <- \(x) {
  dim(x$abundance_array)
}

#' @rdname as_panel_array
#' @export
dimnames.panel_array <- \(x) {
  dimnames(x$abundance_array)
}
