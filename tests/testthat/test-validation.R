test_that("A non data.frame input is rejected", {
  expect_error(as_panel_array(list(a = 1)), "requires a long-format data.frame")
  expect_error(as_panel_array(NULL), "requires a long-format data.frame")
})

test_that("A missing required column is rejected", {
  panel_table <- make_test_panel_table()
  incomplete_table <- panel_table[, setdiff(names(panel_table), "elapsed_time"), drop = FALSE]
  expect_error(as_panel_array(incomplete_table), "is missing required columns")
  expect_error(as_observation_meta(incomplete_table), "is missing required columns")
  expect_error(as_time_meta(incomplete_table), "is missing required columns")
})

test_that("A non numeric abundance column is rejected", {
  panel_table <- make_test_panel_table()
  panel_table$abundance_value <- as.character(panel_table$abundance_value)
  expect_error(as_panel_array(panel_table), "must be numeric")
})

test_that("Missing values are rejected", {
  panel_table <- make_test_panel_table()
  panel_table$abundance_value[1L] <- NA_real_
  expect_error(as_panel_array(panel_table), "must not contain missing values")
})

test_that("A duplicated composite key is rejected", {
  panel_table <- make_test_panel_table()
  duplicated_table <- rbind(panel_table, panel_table[1L, , drop = FALSE])
  expect_error(as_panel_array(duplicated_table), "Duplicate composite key")
})

test_that("An unbalanced panel with missing cells is rejected", {
  panel_table <- make_test_panel_table()
  expect_error(
    as_panel_array(panel_table[-1L, , drop = FALSE]),
    "Unbalanced panel"
  )
})

test_that("An unbalanced panel missing an observation is rejected", {
  panel_table <- make_test_panel_table(time_point_count = 3L)
  drop_flag <- panel_table$time_identifier == "time_002" &
    panel_table$observation == "sample_003"
  expect_error(
    as_panel_array(panel_table[!drop_flag, , drop = FALSE]),
    "Unbalanced panel"
  )
})

test_that("Misaligned dimnames are rejected", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  shifted_observation_meta <- panel_object$observation_meta
  shifted_observation_meta$observation_id <- rev(shifted_observation_meta$observation_id)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = shifted_observation_meta,
      time_meta = panel_object$time_meta
    ),
    "must align exactly"
  )
  shifted_time_meta <- panel_object$time_meta
  shifted_time_meta$time_identifier <- rev(shifted_time_meta$time_identifier)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = panel_object$observation_meta,
      time_meta = shifted_time_meta
    ),
    "must align exactly"
  )
})

test_that("Array dimensions inconsistent with the metadata are rejected", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  trimmed_meta <- panel_object$observation_meta[1L, , drop = FALSE]
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = trimmed_meta,
      time_meta = panel_object$time_meta
    ),
    "Array dimensions disagree with the metadata key lengths"
  )
})

test_that("A non three dimensional array is rejected", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = matrix(1, 2L, 2L),
      observation_meta = panel_object$observation_meta,
      time_meta = panel_object$time_meta
    ),
    "three-dimensional array"
  )
})

test_that("A malformed time_identifier is rejected", {
  panel_table <- make_test_panel_table()
  panel_table$time_identifier <- sub("^time_", "T", panel_table$time_identifier)
  expect_error(as_panel_array(panel_table), "time_identifier must match")
})

test_that("A duplicated metadata key is rejected", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  duplicated_meta <- rbind(
    panel_object$observation_meta,
    panel_object$observation_meta[1L, , drop = FALSE]
  )
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = duplicated_meta,
      time_meta = panel_object$time_meta
    ),
    "must be unique"
  )
})

test_that("Wrong metadata column types are rejected", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  layer_count <- nrow(panel_object$time_meta)
  numeric_intervention_meta <- panel_object$time_meta
  numeric_intervention_meta$is_intervention_timepoint <- rep(c(0, 1), length.out = layer_count)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = panel_object$observation_meta,
      time_meta = numeric_intervention_meta
    ),
    "must be logical"
  )
  character_elapsed_meta <- panel_object$time_meta
  character_elapsed_meta$elapsed_time <- rep(c("0", "7"), length.out = layer_count)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = panel_object$observation_meta,
      time_meta = character_elapsed_meta
    ),
    "must be numeric"
  )
  character_delta_meta <- panel_object$time_meta
  character_delta_meta$delta_t <- rep(c("0", "7"), length.out = layer_count)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = panel_object$observation_meta,
      time_meta = character_delta_meta
    ),
    "delta_t must be numeric"
  )
})

test_that("A missing archive file is rejected", {
  expect_error(
    read_panel_table(file.path(tempdir(), "definitely_missing_archive.tar.gz")),
    "The archive file does not exist"
  )
})

test_that("An archive without manifest.txt is rejected", {
  staging_dir <- tempfile("panelio_no_manifest_")
  dir.create(staging_dir)
  system2("tar", c("-czf", file.path(tempdir(), "no_manifest.tar.gz"), "-C", staging_dir, "."))
  expect_error(
    read_panel_table(file.path(tempdir(), "no_manifest.tar.gz")),
    "is missing manifest.txt"
  )
  unlink(staging_dir, recursive = TRUE)
})

test_that("Inconsistent observation_id sets are rejected", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = 0",
    "observation_id = sample_002; label = 1",
    "observation_id = sample_003; label = 2",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "time_identifier = time_002; elapsed_time = 7; delta_t = 7; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv",
    "time_002.csv"
  )
  layer_tables <- list(
    "time_001.csv" = make_composition_block(
      c("sample_001", "sample_002", "sample_003"),
      "feature_001",
      c(1, 2, 3)
    ),
    "time_002.csv" = make_composition_block(
      c("sample_001", "sample_003"),
      "feature_001",
      c(4, 5)
    )
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "unbalanced_v1.tar.gz")
  )
  expect_error(read_panel_table(archive_path), "Unbalanced panel")
})

test_that("Inconsistent feature columns are rejected", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = 0",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "time_identifier = time_002; elapsed_time = 7; delta_t = 7; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv",
    "time_002.csv"
  )
  layer_tables <- list(
    "time_001.csv" = make_composition_block("sample_001", "feature_001", 1),
    "time_002.csv" = make_composition_block("sample_001", "feature_002", 2)
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "feature_mismatch_v1.tar.gz")
  )
  expect_error(read_panel_table(archive_path), "Unbalanced panel")
})

test_that("A file_index that disagrees with time_meta is rejected", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = 0",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_002.csv"
  )
  layer_tables <- list(
    "time_002.csv" = make_composition_block("sample_001", "feature_001", 1)
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "index_mismatch_v1.tar.gz")
  )
  expect_error(read_panel_table(archive_path), "The file_index block is missing these file entries")
})

test_that("A manifest missing a required field is rejected", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv"
  )
  layer_tables <- list(
    "time_001.csv" = make_composition_block("sample_001", "feature_001", 1)
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "missing_label_v1.tar.gz")
  )
  expect_error(read_panel_table(archive_path), "is missing required columns: label")
})

test_that("An observation absent from the metadata is rejected", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = 0",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv"
  )
  layer_tables <- list(
    "time_001.csv" = make_composition_block(
      c("sample_001", "sample_999"),
      "feature_001",
      c(1, 2)
    )
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "unknown_observation_v1.tar.gz")
  )
  expect_error(read_panel_table(archive_path), "Unbalanced panel")
})

test_that("Argument validation", {
  panel_table <- make_test_panel_table()
  valid_observation_meta <- data.frame(
    observation_id = "sample_001",
    label = factor("ctrl", levels = "ctrl"),
    stringsAsFactors = FALSE
  )
  valid_time_meta <- data.frame(
    time_identifier = "time_001",
    elapsed_time = 0,
    delta_t = 0,
    is_intervention_timepoint = FALSE,
    stringsAsFactors = FALSE
  )
  expect_error(write_panel_archive(panel_table, "", "v1"), "dataset_name")
  expect_error(write_panel_archive(panel_table, "name", ""), "version")
  expect_error(
    write_panel_archive(panel_table, "name", "v1", manifest_comments = NA_character_),
    "manifest_comments"
  )
  expect_error(
    write_manifest(
      observation_meta = valid_observation_meta,
      time_meta = valid_time_meta,
      file_index = "time_001.csv",
      outfile = ""
    ),
    "outfile"
  )
  expect_error(
    write_manifest(
      observation_meta = data.frame(observation_id = "sample_001"),
      time_meta = valid_time_meta,
      file_index = "time_001.csv",
      outfile = tempfile()
    ),
    "is missing required columns: label"
  )
  expect_error(
    write_manifest(
      observation_meta = valid_observation_meta,
      time_meta = valid_time_meta,
      file_index = character(0),
      outfile = tempfile()
    ),
    "file_index"
  )
})
