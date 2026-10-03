test_that("Converting a long table to panel_array preserves dim and dimnames", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)

  expect_s3_class(panel_object, "panel_array")
  expect_identical(dim(panel_object$abundance_array), c(4L, 4L, 3L))
  expect_identical(
    dimnames(panel_object$abundance_array)[[1L]],
    sprintf("sample_%0*d", 3L, 1:4)
  )
  expect_identical(
    dimnames(panel_object$abundance_array)[[2L]],
    sprintf("feature_%0*d", 3L, 1:4)
  )
  expect_identical(
    dimnames(panel_object$abundance_array)[[3L]],
    sprintf("time_%0*d", 3L, 1:3)
  )
  expect_identical(
    panel_object$observation_meta$observation_id,
    sprintf("sample_%0*d", 3L, 1:4)
  )
  expect_identical(panel_object$time_meta$elapsed_time, c(0, 7, 14))
  expect_identical(panel_object$time_meta$delta_t, c(0, 7, 7))
  expect_identical(panel_object$time_meta$is_intervention_timepoint, c(FALSE, FALSE, TRUE))
})

test_that("The bidirectional equivalence contract holds", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  rebuilt_table <- as_panel_table(
    x = NULL,
    abundance_array = panel_object$abundance_array,
    observation_meta = as_observation_meta(panel_table),
    time_meta = as_time_meta(panel_table)
  )
  expect_identical(rebuilt_table, panel_table)
})

test_that("A shuffled row order is restored to the canonical order", {
  panel_table <- make_test_panel_table()
  shuffled_index <- rev(seq_len(nrow(panel_table)))
  shuffled_table <- panel_table[shuffled_index, , drop = FALSE]
  rownames(shuffled_table) <- NULL
  panel_object <- as_panel_array(shuffled_table)
  rebuilt_table <- as_panel_table(
    x = NULL,
    abundance_array = panel_object$abundance_array,
    observation_meta = as_observation_meta(shuffled_table),
    time_meta = as_time_meta(shuffled_table)
  )
  expect_identical(
    panelio:::sort_panel_table(rebuilt_table),
    panelio:::sort_panel_table(panel_table)
  )
})

test_that("as_panel_table supports its two entry modes", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  archive_path <- write_panel_archive(
    panel_object,
    dataset_name = "unit",
    version = "v1"
  )
  from_archive <- as_panel_table(archive_path)
  from_array <- as_panel_table(
    x = NULL,
    abundance_array = panel_object$abundance_array,
    observation_meta = as_observation_meta(panel_table),
    time_meta = as_time_meta(panel_table)
  )
  expect_identical(from_archive, from_array)
  expect_identical(from_archive, panel_table)
})

test_that("as_panel_table errors when required arguments are missing", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  expect_error(
    as_panel_table(x = NULL, abundance_array = panel_object$abundance_array),
    "must all be supplied"
  )
  expect_error(as_panel_table(x = 42), "must be an archive path or NULL")
})

test_that("as_observation_meta and as_time_meta accept long tables, panel_array objects and manifest text", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  observation_meta <- as_observation_meta(panel_table)
  time_meta <- as_time_meta(panel_table)

  expect_identical(
    observation_meta$observation_id,
    sprintf("sample_%0*d", 3L, 1:4)
  )
  expect_true(is.factor(observation_meta$label))
  expect_identical(levels(observation_meta$label), c("ctrl", "treat"))
  expect_identical(
    as.integer(observation_meta$label),
    c(1L, 2L, 1L, 2L)
  )
  expect_identical(time_meta$time_identifier, sprintf("time_%0*d", 3L, 1:3))
  expect_identical(time_meta$delta_t, c(0, 7, 7))
  expect_identical(as_observation_meta(panel_object), observation_meta)
  expect_identical(as_time_meta(panel_object), time_meta)

  manifest_lines <- make_reference_manifest_lines()
  expect_identical(as_observation_meta(manifest_lines)$observation_id, c("sample_001", "sample_002"))
  expect_identical(as_time_meta(manifest_lines)$elapsed_time, c(0, 7))
  expect_error(as_observation_meta(42), "accepts only")
  expect_error(as_time_meta(42), "accepts only")
})

test_that("A character label is rejected in a long table but accepted in a manifest", {
  character_table <- make_test_panel_table(
    observation_count = 3L,
    feature_count = 3L,
    time_point_count = 2L,
    label_values = c("control", "treated", "control")
  )
  character_table$label <- as.character(character_table$label)
  expect_error(as_panel_array(character_table), "must be an R factor")
  expect_error(as_observation_meta(character_table), "must be an R factor")

  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = control",
    "observation_id = sample_002; label = treated",
    "observation_id = sample_003; label = control"
  )
  observation_meta <- as_observation_meta(manifest_lines)
  expect_true(is.factor(observation_meta$label))
  expect_identical(levels(observation_meta$label), c("control", "treated"))
  expect_identical(nlevels(droplevels(observation_meta$label)), 2L)
})

test_that("panel_array provides dim, dimnames and print methods", {
  panel_object <- as_panel_array(make_test_panel_table())
  expect_identical(dim(panel_object), c(4L, 4L, 3L))
  expect_identical(
    dimnames(panel_object),
    list(
      sprintf("sample_%0*d", 3L, 1:4),
      sprintf("feature_%0*d", 3L, 1:4),
      sprintf("time_%0*d", 3L, 1:3)
    )
  )
  expect_output(print(panel_object), "panel_array")
})
