# =============================================================================
# R/read_panel.R
# All read_* entry points of the panel archive.
# =============================================================================

#' Read a balanced longitudinal panel archive
#'
#' Extracts the tar.gz archive, machine parses `manifest.txt` to obtain
#' `observation_meta` and `time_meta`, reads every `time_xxx.csv` wide table,
#' converts the wide tables to the long format, merges the two metadata tables
#' and returns the long-format data.frame.
#'
#' Archive level contract:
#' * the `observation_id` sets of all `time_*.csv` must be identical;
#' * that set must match the `observation_meta` of the manifest exactly;
#' * the feature column names and their order must be identical across layers;
#' * any parse or validation failure aborts immediately.
#'
#' @param archive_path Path to the tar.gz archive.
#'
#' @return A long-format data.frame with the columns `observation`,
#'   `time_identifier`, `elapsed_time`, `feature`, `abundance_value`, `label`
#'   and `is_intervention`. Rows follow the observation x time_identifier x
#'   feature order with feature varying fastest.
#' @export
read_panel_table <- \(archive_path) {
  assert_scalar_text(archive_path, "archive_path")

  if (!file.exists(archive_path)) {
    panel_abort(sprintf("The archive file does not exist: %s", archive_path))
  }

  staging_dir <- tempfile("panelio_read_")
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
  extract_archive(archive_path, staging_dir)
  manifest_path <- file.path(staging_dir, MANIFEST_FILE_NAME)

  if (!file.exists(manifest_path)) {
    panel_abort(sprintf("The archive is missing %s.", MANIFEST_FILE_NAME))
  }

  manifest_lines <- readLines(
    manifest_path,
    warn = FALSE,
    encoding = CSV_FILE_ENCODING
  )
  manifest <- parse_manifest(manifest_lines)
  observation_meta <- manifest$observation_meta
  time_meta <- manifest$time_meta
  assert_observation_meta(observation_meta)
  assert_time_meta(time_meta)
  layer_paths <- resolve_layer_paths(
    staging_dir = staging_dir,
    time_meta = time_meta,
    file_index = manifest$file_index
  )
  layer_tables <- lapply(layer_paths, read_csv_table)
  assemble_long_table(layer_tables, observation_meta, time_meta)
}

# ---- internals --------------------------------------------------------------
resolve_layer_paths <- \(staging_dir, time_meta, file_index) {
  time_levels <- as.character(time_meta[[TIME_IDENTIFIER_COLUMN]])
  expected_index <- time_file_name(time_levels)
  index_mismatch <- which(!expected_index %in% file_index)

  if (length(index_mismatch) > 0L) {
    panel_abort(sprintf(
      "The file_index block is missing these file entries: %s.",
      truncate_report(
        expected_index[index_mismatch],
        "unknown_id_report_count"
      )
    ))
  }

  orphan_index <- setdiff(file_index, expected_index)

  if (length(orphan_index) > 0L) {
    panel_abort(sprintf(
      "The file_index block lists files that are not registered in time_meta: %s.",
      truncate_report(orphan_index, "unknown_id_report_count")
    ))
  }

  file.path(staging_dir, expected_index)
}

validate_layer_tables <- \(layer_tables, observation_levels) {
  reference_names <- names(layer_tables[[1L]])

  if (length(reference_names) < 2L) {
    panel_abort(
      "A time layer wide table needs at least the observation_id column ",
      "plus one feature column."
    )
  }

  if (!identical(reference_names[1L], OBSERVATION_ID_COLUMN)) {
    panel_abort(sprintf(
      "The first column of a time layer wide table must be %s but was %s.",
      OBSERVATION_ID_COLUMN,
      reference_names[1L]
    ))
  }

  feature_names <- reference_names[-1L]

  if (anyDuplicated(reference_names) || any(!nzchar(feature_names))) {
    panel_abort("Time layer wide table column names must be non-empty and unique.")
  }

  inconsistent_layer <- which(!vapply(
    layer_tables,
    \(layer_table) identical(names(layer_table), reference_names),
    logical(1L)
  ))

  if (length(inconsistent_layer) > 0L) {
    panel_abort(sprintf(
      paste0(
        "Unbalanced panel detected: every time layer must share the same feature ",
        "column names and order; offending layer indices: %s."
      ),
      truncate_report(inconsistent_layer, "missing_cell_report_count")
    ))
  }

  layer_observation_levels <- lapply(
    layer_tables,
    \(layer_table) as.character(layer_table[[OBSERVATION_ID_COLUMN]])
  )
  unbalanced_layer <- which(!vapply(
    layer_observation_levels,
    \(layer_levels) setequal(layer_levels, observation_levels),
    logical(1L)
  ))

  if (length(unbalanced_layer) > 0L) {
    panel_abort(sprintf(
      paste0(
        "Unbalanced panel detected: every time layer must share the same ",
        "observation_id set; offending layer indices: %s."
      ),
      truncate_report(unbalanced_layer, "missing_cell_report_count")
    ))
  }

  duplicated_layer <- which(vapply(
    layer_observation_levels,
    \(layer_levels) anyDuplicated(layer_levels) > 0L,
    logical(1L)
  ))

  if (length(duplicated_layer) > 0L) {
    panel_abort(sprintf(
      "observation_id must be unique inside a time layer wide table; offending layer indices: %s.",
      truncate_report(duplicated_layer, "missing_cell_report_count")
    ))
  }

  observation_levels
}

assemble_long_table <- \(layer_tables, observation_meta, time_meta) {
  observation_levels <- as.character(observation_meta[[OBSERVATION_ID_COLUMN]])
  time_levels <- as.character(time_meta[[TIME_IDENTIFIER_COLUMN]])
  validate_layer_tables(layer_tables, observation_levels)
  aligned_layer_tables <- Map(
    \(layer_table) {
      row_index <- match(observation_levels, as.character(layer_table[[OBSERVATION_ID_COLUMN]]))
      aligned_table <- layer_table[row_index, , drop = FALSE]
      rownames(aligned_table) <- NULL
      aligned_table
    },
    layer_tables
  )
  layer_long_tables <- Map(panel_wide_to_long, aligned_layer_tables, time_levels)
  long_table <- do.call(rbind, layer_long_tables)
  rownames(long_table) <- NULL
  attach_metadata(long_table, observation_meta, time_meta)
}
