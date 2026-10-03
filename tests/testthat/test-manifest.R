test_that("Manifest blocks parse into canonical key value data frames", {
  manifest_lines <- make_reference_manifest_lines()
  blocks <- panelio:::parse_manifest_blocks(manifest_lines)

  expect_named(blocks, c("observation_meta", "time_meta", "file_index"))
  expect_identical(blocks$file_index, c("time_001.csv", "time_002.csv"))

  observation_meta <- panelio:::parse_manifest_pairs(blocks$observation_meta)
  expect_identical(names(observation_meta), c("observation_id", "label"))
  expect_type(observation_meta$observation_id, "character")
  expect_type(observation_meta$label, "character")

  time_meta <- panelio:::parse_manifest_pairs(blocks$time_meta)
  expect_identical(
    names(time_meta),
    c("time_identifier", "elapsed_time", "delta_t", "is_intervention_timepoint")
  )
  expect_type(time_meta$elapsed_time, "double")
  expect_type(time_meta$delta_t, "double")
  expect_type(time_meta$is_intervention_timepoint, "logical")
  expect_identical(time_meta$elapsed_time, c(0, 7))
  expect_identical(time_meta$delta_t, c(0, 7))
  expect_identical(time_meta$is_intervention_timepoint, c(FALSE, TRUE))
})

test_that("Comment lines and blank lines are ignored", {
  manifest_lines <- c(
    "# leading comment",
    "   ",
    "[time_meta]",
    "# indented comment is still a comment",
    "time_identifier = time_001; elapsed_time = 0; is_intervention_timepoint = FALSE",
    "",
    "[file_index]",
    "time_001.csv"
  )
  blocks <- panelio:::parse_manifest_blocks(manifest_lines)
  time_meta <- panelio:::parse_manifest_pairs(blocks$time_meta)
  expect_identical(nrow(time_meta), 1L)
  expect_identical(time_meta$time_identifier, "time_001")
})

test_that("A character label stays character typed and becomes a grouping factor", {
  manifest_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = control",
    "observation_id = sample_002; label = treatment"
  )
  blocks <- panelio:::parse_manifest_blocks(manifest_lines)
  observation_meta <- panelio:::parse_manifest_pairs(blocks$observation_meta)
  expect_type(observation_meta$label, "character")
  expect_identical(observation_meta$label, c("control", "treatment"))

  grouped_meta <- as_observation_meta(manifest_lines)
  expect_true(is.factor(grouped_meta$label))
  expect_identical(levels(grouped_meta$label), c("control", "treatment"))
})

test_that("A missing block aborts immediately", {
  broken_lines <- c("[observation_meta]", "observation_id = sample_001; label = 0")
  expect_error(
    panelio:::parse_manifest(broken_lines),
    "missing block"
  )
})

test_that("A missing block marker aborts immediately", {
  expect_error(
    panelio:::parse_manifest_blocks(c("observation_id = sample_001")),
    "no [block] marker",
    fixed = TRUE
  )
})

test_that("A duplicated block aborts immediately", {
  duplicated_lines <- c(
    "[observation_meta]",
    "observation_id = sample_001; label = 0",
    "[observation_meta]",
    "observation_id = sample_002; label = 1",
    "[time_meta]",
    "time_identifier = time_001; elapsed_time = 0; is_intervention_timepoint = FALSE",
    "[file_index]",
    "time_001.csv"
  )
  expect_error(panelio:::parse_manifest(duplicated_lines), "is defined more than once")
})

test_that("Text produced by write_manifest can be parsed again", {
  observation_meta <- data.frame(
    observation_id = c("sample_001", "sample_002"),
    label = factor(c("ctrl", "treat"), levels = c("ctrl", "treat")),
    stringsAsFactors = FALSE
  )
  time_meta <- data.frame(
    time_identifier = c("time_001", "time_002"),
    elapsed_time = c(0, 7),
    delta_t = c(0, 7),
    is_intervention_timepoint = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  outfile <- tempfile(fileext = ".txt")
  write_manifest(
    observation_meta = observation_meta,
    time_meta = time_meta,
    file_index = c("time_001.csv", "time_002.csv"),
    comments = "unit test",
    outfile = outfile
  )
  written_lines <- readLines(outfile)
  expect_true(any(startsWith(written_lines, "# unit test")))
  expect_true(any(written_lines == "[observation_meta]"))
  expect_true(any(written_lines == "[time_meta]"))
  expect_true(any(written_lines == "[file_index]"))
  expect_true(any(written_lines == "observation_id = sample_001; label = ctrl"))
  expect_true(any(written_lines == paste(
    "time_identifier = time_002; elapsed_time = 7; delta_t = 7;",
    "is_intervention_timepoint = TRUE"
  )))
  parsed <- panelio:::parse_manifest(written_lines)
  expect_identical(
    as.character(parsed$observation_meta$label),
    as.character(observation_meta$label)
  )
  expect_identical(parsed$time_meta, time_meta)
  expect_identical(parsed$file_index, c("time_001.csv", "time_002.csv"))
  unlink(outfile)
})
