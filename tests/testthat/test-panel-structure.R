test_that("label is the grouping variable and must be an R factor", {
  panel_table <- make_test_panel_table()
  expect_true(is.factor(panel_table$label))

  character_table <- panel_table
  character_table$label <- as.character(character_table$label)
  expect_error(as_panel_array(character_table), "must be an R factor")

  numeric_table <- panel_table
  numeric_table$label <- rep(c(1, 2), length.out = nrow(numeric_table))
  expect_error(as_panel_array(numeric_table), "must be an R factor")
})

test_that("The number of label groups is capped", {
  two_group_table <- make_test_panel_table(
    observation_count = 4L,
    feature_count = 3L,
    time_point_count = 2L,
    group_count = 2L
  )
  expect_silent(as_panel_array(two_group_table))
  expect_identical(nlevels(droplevels(two_group_table$label)), 2L)

  four_group_table <- make_test_panel_table(
    observation_count = 4L,
    feature_count = 3L,
    time_point_count = 2L,
    group_count = 4L
  )
  expect_silent(as_panel_array(four_group_table))
  expect_identical(nlevels(droplevels(four_group_table$label)), 4L)

  five_group_table <- make_test_panel_table(
    observation_count = 5L,
    feature_count = 3L,
    time_point_count = 2L,
    group_count = 5L
  )
  expect_error(as_panel_array(five_group_table), "at most 4 groups are allowed")
})

test_that("A manifest label becomes a grouping factor with few levels", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = ctrl",
    "observation_id = sample_002; label = treat",
    "observation_id = sample_003; label = ctrl",
    "observation_id = sample_004; label = treat"
  )
  observation_meta <- as_observation_meta(manifest_lines)
  expect_true(is.factor(observation_meta$label))
  expect_identical(levels(observation_meta$label), c("ctrl", "treat"))
  expect_identical(
    as.integer(observation_meta$label),
    c(1L, 2L, 1L, 2L)
  )
})

test_that("A non monotone intervention flag is rejected", {
  panel_table <- make_test_panel_table(
    time_point_count = 4L,
    intervention_time_index = 2L
  )
  time_layer_flag <- !duplicated(panel_table$time_identifier)
  time_levels <- unique(panel_table$time_identifier)
  expect_identical(
    panel_table$is_intervention[time_layer_flag],
    c(FALSE, TRUE, TRUE, TRUE)
  )
  expect_silent(as_panel_array(panel_table))

  # TRUE must never fall back to FALSE.
  broken_table <- panel_table
  broken_table$is_intervention[broken_table$time_identifier == time_levels[4L]] <- FALSE
  expect_error(
    as_panel_array(broken_table),
    "must stay TRUE once it turns TRUE"
  )

  # The earliest time layer must already be FALSE.
  early_true_table <- panel_table
  early_true_table$is_intervention[early_true_table$time_identifier == time_levels[1L]] <- TRUE
  expect_error(
    as_panel_array(early_true_table),
    "must be FALSE for the earliest time layer"
  )
})

test_that("A panel with no intervention at all is accepted", {
  panel_table <- make_test_panel_table(time_point_count = 3L, intervention_time_index = 99L)
  time_layer_flag <- !duplicated(panel_table$time_identifier)
  expect_false(any(panel_table$is_intervention[time_layer_flag]))
  expect_silent(as_panel_array(panel_table))
})

test_that("Every observation and time layer must be a composition", {
  panel_table <- make_test_panel_table()
  composition_sum <- tapply(
    panel_table$abundance_value,
    list(panel_table$observation, panel_table$time_identifier),
    sum
  )
  expect_true(all(abs(composition_sum - 1) < 1e-9))

  doubled_table <- panel_table
  doubled_table$abundance_value[doubled_table$observation == "sample_001"] <-
    doubled_table$abundance_value[doubled_table$observation == "sample_001"] * 2
  expect_error(as_panel_array(doubled_table), "must be compositional")

  shifted_table <- panel_table
  shifted_table$abundance_value <- shifted_table$abundance_value + 0.1
  expect_error(as_panel_array(shifted_table), "must be compositional")
})

test_that("The array pathway enforces the composition as well", {
  panel_table <- make_test_panel_table()
  panel_object <- as_panel_array(panel_table)
  broken_array <- panel_object$abundance_array
  broken_array[1L, , 1L] <- broken_array[1L, , 1L] + 0.5
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = broken_array,
      observation_meta = panel_object$observation_meta,
      time_meta = panel_object$time_meta
    ),
    "must be compositional"
  )
})

test_that("A sparse composition is accepted and stays sparse", {
  panel_table <- make_test_panel_table(
    feature_count = 6L,
    nonzero_feature_count = 2L
  )
  zero_count <- tapply(
    panel_table$abundance_value == 0,
    list(panel_table$observation, panel_table$time_identifier),
    sum
  )
  expect_true(all(zero_count == 4L))
  panel_object <- as_panel_array(panel_table)
  expect_silent(as_panel_table(
    x = NULL,
    abundance_array = panel_object$abundance_array,
    observation_meta = panel_object$observation_meta,
    time_meta = panel_object$time_meta
  ))
})

test_that("delta_t is derived from elapsed_time and may be non uniform", {
  panel_table <- make_test_panel_table(
    time_point_count = 4L,
    elapsed_time_values = c(0, 3, 7, 21),
    intervention_time_index = 2L
  )
  time_meta <- as_time_meta(panel_table)
  expect_identical(
    time_meta$elapsed_time,
    c(0, 3, 7, 21)
  )
  expect_identical(time_meta$delta_t, c(0, 3, 4, 14))
  expect_true(length(unique(time_meta$delta_t)) > 1L)
  expect_identical(
    names(time_meta),
    c("time_identifier", "elapsed_time", "delta_t", "is_intervention_timepoint")
  )
  expect_silent(as_panel_array(panel_table))
})

test_that("delta_t is derived in elapsed_time order even for a shuffled table", {
  panel_table <- make_test_panel_table(
    time_point_count = 4L,
    elapsed_time_values = c(0, 3, 7, 21),
    intervention_time_index = 2L
  )
  shuffled_table <- panel_table[rev(seq_len(nrow(panel_table))), , drop = FALSE]
  rownames(shuffled_table) <- NULL
  shuffled_meta <- as_time_meta(shuffled_table)
  shuffled_meta <- shuffled_meta[order(shuffled_meta$time_identifier), , drop = FALSE]
  rownames(shuffled_meta) <- NULL
  expect_identical(shuffled_meta, as_time_meta(panel_table))
})

test_that("An inconsistent delta_t is rejected", {
  panel_table <- make_test_panel_table(time_point_count = 3L)
  panel_object <- as_panel_array(panel_table)
  broken_time_meta <- panel_object$time_meta
  broken_time_meta$delta_t <- c(0, 7, 99)
  expect_error(
    as_panel_table(
      x = NULL,
      abundance_array = panel_object$abundance_array,
      observation_meta = panel_object$observation_meta,
      time_meta = broken_time_meta
    ),
    "delta_t must equal the elapsed_time difference"
  )
})

test_that("A manifest time_meta carries delta_t", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = ctrl",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; delta_t = 0; is_intervention_timepoint = FALSE",
    "time_identifier = time_002; elapsed_time = 3; delta_t = 3; is_intervention_timepoint = FALSE",
    "time_identifier = time_003; elapsed_time = 7; delta_t = 4; is_intervention_timepoint = TRUE",
    "[file_index]",
    "time_001.csv",
    "time_002.csv",
    "time_003.csv"
  )
  time_meta <- as_time_meta(manifest_lines)
  expect_identical(time_meta$delta_t, c(0, 3, 4))
  expect_identical(time_meta$is_intervention_timepoint, c(FALSE, FALSE, TRUE))
})
