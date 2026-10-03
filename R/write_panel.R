# =============================================================================
# R/write_panel.R
# All write_* entry points of the panel archive.
# Packaging never calls setwd(): the external tar receives an explicit
# -C <staging_dir> argument instead.
# =============================================================================

#' Write a balanced longitudinal panel archive
#'
#' The input is either a long-format data.frame or a `panel_array` object,
#' which is converted internally by `as_panel_array()`. The three dimensional
#' array is split into one `time_xxx.csv` wide table per time layer, the
#' `observation_meta`, `time_meta` and `file_index` blocks are assembled and
#' handed to `write_manifest()` to emit the machine readable `manifest.txt`,
#' and the result is finally packed as
#' `\{dataset_name\}_\{version\}.tar.gz`.
#'
#' @param x A long-format data.frame or a `panel_array` object.
#' @param dataset_name Dataset name used to build the archive file name.
#' @param version Dataset version used to build the archive file name.
#' @param out_dir Output directory of the archive, `tempdir()` by default.
#' @param manifest_comments Optional extra comment lines written at the top of
#'   the manifest, a character vector defaulting to `character(0)`.
#'
#' @return The full file path of the generated tar.gz archive.
#' @export
write_panel_archive <- \(x, dataset_name, version, out_dir = tempdir(), manifest_comments = character(0)) {
  assert_scalar_text(dataset_name, "dataset_name")
  assert_scalar_text(version, "version")

  if (!is.character(manifest_comments)) {
    panel_abort("manifest_comments must be a character vector.")
  }

  if (anyNA(manifest_comments)) {
    panel_abort("manifest_comments must not contain missing values.")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  panel_object <- normalize_panel_input(x)
  time_levels <- as.character(panel_object$time_meta[[TIME_IDENTIFIER_COLUMN]])
  observation_levels <- as.character(
    panel_object$observation_meta[[OBSERVATION_ID_COLUMN]]
  )
  feature_levels <- dimnames(panel_object$abundance_array)[[2L]]
  archive_name <- sprintf(
    ARCHIVE_FILE_TEMPLATE,
    dataset_name = dataset_name,
    version = version
  )
  archive_path <- file.path(out_dir, archive_name)
  staging_dir <- tempfile("panelio_write_")
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
  file_index <- time_file_name(time_levels)
  layer_tables <- build_wide_layers(
    abundance_array = panel_object$abundance_array,
    observation_levels = observation_levels,
    feature_levels = feature_levels,
    time_levels = time_levels
  )
  invisible(Map(
    write_csv_table,
    layer_tables,
    file.path(staging_dir, file_index)
  ))
  write_manifest(
    observation_meta = panel_object$observation_meta,
    time_meta = panel_object$time_meta,
    file_index = file_index,
    comments = manifest_comments,
    outfile = file.path(staging_dir, MANIFEST_FILE_NAME)
  )
  create_archive(archive_path, staging_dir)
}

#' Generate a machine readable manifest.txt
#'
#' Emits the machine parsable manifest text as the `[observation_meta]`,
#' `[time_meta]` and `[file_index]` blocks. Key value pairs use `key = value`
#' and several fields on the same line are separated by semicolons; a leading
#' `#` marks a comment.
#'
#' @param observation_meta Observation metadata data.frame holding
#'   `observation_id` and `label`.
#' @param time_meta Time layer metadata data.frame holding `time_identifier`,
#'   `elapsed_time` and `is_intervention_timepoint`.
#' @param file_index Character vector of the archive file index.
#' @param comments Extra comment lines, a character vector defaulting to
#'   `character(0)`.
#' @param outfile Output file path.
#'
#' @return The output file path, invisibly.
#' @export
write_manifest <- \(observation_meta, time_meta, file_index, comments = character(0), outfile) {
  assert_observation_meta(observation_meta)
  assert_time_meta(time_meta)
  assert_scalar_text(outfile, "outfile")

  if (!is.character(file_index) || length(file_index) == 0L || anyNA(file_index)) {
    panel_abort("file_index must be a non-empty character vector.")
  }

  if (!is.character(comments)) {
    panel_abort("comments must be a character vector.")
  }

  observation_lines <- format_meta_rows(
    observation_meta[, OBSERVATION_META_COLUMNS, drop = FALSE]
  )
  time_lines <- format_meta_rows(
    time_meta[, TIME_META_COLUMNS, drop = FALSE]
  )
  comment_lines <- if (length(comments) > 0L) {
    paste0(COMMENT_PREFIX, " ", trimws(comments))
  } else {
    character(0)
  }
  header_lines <- c(
    paste0(COMMENT_PREFIX, " ", MANIFEST_HEADER_COMMENT),
    comment_lines,
    character(0)
  )
  manifest_lines <- c(
    header_lines,
    paste0("[", OBSERVATION_META_BLOCK, "]"),
    observation_lines,
    character(0),
    paste0("[", TIME_META_BLOCK, "]"),
    time_lines,
    character(0),
    paste0("[", FILE_INDEX_BLOCK, "]"),
    trimws(file_index)
  )
  con <- file(outfile, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(
    enc2utf8(manifest_lines),
    con = con,
    useBytes = TRUE
  )
  invisible(outfile)
}

# ---- internals --------------------------------------------------------------
normalize_panel_input <- \(x) {
  if (inherits(x, "panel_array")) {
    assert_array_meta_alignment(x$abundance_array, x$observation_meta, x$time_meta)
    return(x)
  }

  as_panel_array(x)
}

format_meta_rows <- \(meta_table) {
  column_values <- as.list(meta_table)
  column_names <- names(meta_table)
  names(column_values) <- column_names

  if (nrow(meta_table) == 0L) {
    return(character(0))
  }

  vapply(
    seq_len(nrow(meta_table)),
    \(row_index) {
      row_values <- lapply(
        column_values,
        \(column_value) format_manifest_value(column_value[[row_index]])
      )
      names(row_values) <- column_names
      paste0(
        paste0(
          names(row_values),
          " ",
          KEY_VALUE_SEPARATOR,
          " ",
          unlist(row_values, use.names = FALSE)
        ),
        collapse = MANIFEST_FIELD_JOIN
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

format_abundance_layer <- \(abundance_matrix) {
  flat_value <- as.vector(abundance_matrix)
  formatted_value <- format_numeric_value(flat_value)
  matrix(
    formatted_value,
    nrow = nrow(abundance_matrix),
    ncol = ncol(abundance_matrix),
    dimnames = dimnames(abundance_matrix)
  )
}

build_wide_layers <- \(abundance_array, observation_levels, feature_levels, time_levels) {
  observation_count <- length(observation_levels)
  feature_count <- length(feature_levels)
  observation_frame <- data.frame(
    observation_levels,
    stringsAsFactors = FALSE
  )
  names(observation_frame) <- OBSERVATION_ID_COLUMN

  Map(
    \(time_index) {
      layer_slice <- abundance_array[, , time_index, drop = FALSE]
      layer_matrix <- matrix(
        as.vector(layer_slice),
        nrow = observation_count,
        ncol = feature_count,
        dimnames = list(NULL, feature_levels)
      )
      abundance_frame <- as.data.frame(
        format_abundance_layer(layer_matrix),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      cbind(observation_frame, abundance_frame)
    },
    seq_along(time_levels)
  )
}
