test_that("Read write loop: long table to archive to long table is lossless", {
  panel_table <- make_test_panel_table(
    observation_count = 4L,
    feature_count = 5L,
    time_point_count = 3L
  )
  archive_path <- write_panel_archive(
    panel_table,
    dataset_name = "roundtrip",
    version = "v1"
  )
  expect_true(file.exists(archive_path))
  expect_identical(basename(archive_path), "roundtrip_v1.tar.gz")
  restored_table <- read_panel_table(archive_path)
  expect_identical(restored_table, panel_table)
  expect_identical(names(restored_table), names(panel_table))
})

test_that("The archive file names match the manifest", {
  panel_table <- make_test_panel_table()
  archive_path <- write_panel_archive(panel_table, "layout", "v3")
  extract_dir <- tempfile("panelio_inspect_")
  dir.create(extract_dir)
  system2("tar", c("-xzf", archive_path, "-C", extract_dir))
  archive_entries <- list.files(extract_dir, recursive = TRUE)
  expect_setequal(
    archive_entries,
    c(
      "manifest.txt",
      sprintf("time_%0*d.csv", 3L, 1:3)
    )
  )
  manifest_lines <- readLines(file.path(extract_dir, "manifest.txt"))
  expect_true(any(manifest_lines == "[file_index]"))
  expect_true(any(manifest_lines == "time_001.csv"))
  unlink(extract_dir, recursive = TRUE)
})

test_that("Long table and panel_array inputs produce identical archives", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  from_table <- write_panel_archive(panel_table, "same", "v1")
  from_object <- write_panel_archive(panel_object, "same", "v1")
  expect_identical(read_panel_table(from_table), read_panel_table(from_object))
  expect_identical(read_panel_table(from_table), panel_table)
})

test_that("Rewriting the same path overwrites the previous archive", {
  panel_table <- make_test_panel_table()
  out_dir <- tempfile("panelio_overwrite_")
  dir.create(out_dir)
  first_path <- write_panel_archive(panel_table, "overwrite", "v1", out_dir = out_dir)
  second_path <- write_panel_archive(panel_table, "overwrite", "v1", out_dir = out_dir)
  expect_identical(first_path, second_path)
  expect_identical(read_panel_table(second_path), panel_table)
})

test_that("Manifest comments are stored inside the archive", {
  panel_table <- make_test_panel_table()
  archive_path <- write_panel_archive(
    panel_table,
    "commented",
    "v1",
    manifest_comments = c("operator = alice", "instrument = illumina")
  )
  extract_dir <- tempfile("panelio_comment_")
  dir.create(extract_dir)
  system2("tar", c("-xzf", archive_path, "-C", extract_dir))
  manifest_lines <- readLines(file.path(extract_dir, "manifest.txt"))
  expect_true(any(manifest_lines == "# operator = alice"))
  expect_true(any(manifest_lines == "# instrument = illumina"))
  unlink(extract_dir, recursive = TRUE)
})

test_that("A missing out_dir is created automatically", {
  panel_table <- make_test_panel_table()
  out_dir <- file.path(tempfile("panelio_root_"), "nested", "out")
  archive_path <- write_panel_archive(panel_table, "created", "v1", out_dir = out_dir)
  expect_true(dir.exists(out_dir))
  expect_identical(read_panel_table(archive_path), panel_table)
})

test_that("Floating point abundances keep full precision through the archive", {
  observation_count <- 4L
  feature_count <- 4L
  time_point_count <- 3L
  panel_table <- make_test_panel_table(
    observation_count = observation_count,
    feature_count = feature_count,
    time_point_count = time_point_count
  )
  # The long table is ordered time, then feature fastest, then observation, so
  # one time block holds observation_count * feature_count rows.
  tricky_patterns <- rbind(
    c(1 / 3, 2 / 3, 0, 0),
    c(1e-20, 1, 0, 0),
    c(0.1, 0.9, 1e-300, 0),
    c(pi / 4, 1 - pi / 4, 0, 0)
  )
  block_size <- observation_count * feature_count
  row_offset <- seq_len(nrow(panel_table)) - 1L
  row_position <- row_offset %% block_size
  row_feature <- row_position %% observation_count + 1L
  row_observation <- row_position %/% observation_count + 1L
  cell_index <- (row_offset %/% block_size) * observation_count + row_observation
  cell_index <- (cell_index - 1L) %% nrow(tricky_patterns) + 1L
  panel_table$abundance_value <- tricky_patterns[cbind(cell_index, row_feature)]
  expect_true(is.double(panel_table$abundance_value))
  cell_sum <- tapply(
    panel_table$abundance_value,
    list(panel_table$observation, panel_table$time_identifier),
    sum
  )
  expect_true(all(abs(cell_sum - 1) < 1e-9))
  archive_path <- write_panel_archive(panel_table, "precision", "v1")
  restored_table <- read_panel_table(archive_path)
  expect_identical(restored_table$abundance_value, panel_table$abundance_value)
})

test_that("Reading a hand built reference archive", {
  manifest_lines <- make_reference_manifest_lines()
  feature_levels <- c("feature_001", "feature_002")
  layer_tables <- list(
    "time_001.csv" = make_composition_block(
      c("sample_001", "sample_002"),
      feature_levels,
      c(1, 2, 3, 4)
    ),
    "time_002.csv" = make_composition_block(
      c("sample_001", "sample_002"),
      feature_levels,
      c(5, 6, 7, 8)
    )
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "reference_v1.tar.gz")
  )
  restored_table <- read_panel_table(archive_path)
  expect_identical(
    names(restored_table),
    c(
      "observation",
      "time_identifier",
      "elapsed_time",
      "feature",
      "abundance_value",
      "label",
      "is_intervention"
    )
  )
  expect_identical(nrow(restored_table), 8L)
  expect_identical(
    unique(restored_table$observation),
    c("sample_001", "sample_002")
  )
  expect_identical(restored_table$elapsed_time, rep(c(0, 7), each = 4L))
  expect_identical(
    as.character(restored_table$label),
    rep(c("ctrl", "treat"), each = 2L, times = 2L)
  )
  expect_identical(restored_table$is_intervention, rep(c(FALSE, TRUE), each = 4L))
  expect_identical(
    restored_table$abundance_value,
    c(1 / 3, 2 / 3, 3 / 7, 4 / 7, 5 / 11, 6 / 11, 7 / 15, 8 / 15)
  )
  time_meta <- as_time_meta(restored_table)
  expect_identical(time_meta$delta_t, c(0, 7))
})

test_that("The archive observation order follows the manifest", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_002; label = treat",
    "observation_id = sample_001; label = ctrl",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv"
  )
  layer_tables <- list(
    "time_001.csv" = make_composition_block(
      c("sample_001", "sample_002"),
      c("feature_001", "feature_002"),
      c(11, 11, 5, 15)
    )
  )
  archive_path <- write_raw_archive(
    manifest_lines,
    layer_tables,
    file.path(tempdir(), "reordered_v1.tar.gz")
  )
  restored_table <- read_panel_table(archive_path)
  expect_identical(
    restored_table$observation,
    c("sample_002", "sample_002", "sample_001", "sample_001")
  )
  expect_identical(restored_table$abundance_value, c(0.25, 0.75, 0.5, 0.5))
  expect_identical(as.character(restored_table$label), c("treat", "treat", "ctrl", "ctrl"))
})
