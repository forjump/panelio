# =============================================================================
# R/utils_panel.R
# Internal validation helpers, manifest grammar helpers and structural
# conversions shared by the reader, the writer and the converters.
# Every function defined here is internal: none of them is exported.
# All tunable values come from R/constants.R; there is no magic number here.
# =============================================================================

# ---- Column name constants ---------------------------------------------------
OBSERVATION_COLUMN <- "observation"
TIME_IDENTIFIER_COLUMN <- "time_identifier"
ELAPSED_TIME_COLUMN <- "elapsed_time"
FEATURE_COLUMN <- "feature"
ABUNDANCE_VALUE_COLUMN <- "abundance_value"
LABEL_COLUMN <- "label"
IS_INTERVENTION_COLUMN <- "is_intervention"
OBSERVATION_ID_COLUMN <- "observation_id"
IS_INTERVENTION_META_COLUMN <- "is_intervention_timepoint"
DELTA_T_COLUMN <- "delta_t"

LONG_TABLE_COLUMNS <- c(
  OBSERVATION_COLUMN,
  TIME_IDENTIFIER_COLUMN,
  ELAPSED_TIME_COLUMN,
  FEATURE_COLUMN,
  ABUNDANCE_VALUE_COLUMN,
  LABEL_COLUMN,
  IS_INTERVENTION_COLUMN
)

OBSERVATION_META_COLUMNS <- c(OBSERVATION_ID_COLUMN, LABEL_COLUMN)

TIME_META_COLUMNS <- c(
  TIME_IDENTIFIER_COLUMN,
  ELAPSED_TIME_COLUMN,
  DELTA_T_COLUMN,
  IS_INTERVENTION_META_COLUMN
)

# ---- Error helpers -----------------------------------------------------------
# Message parts are concatenated so that long diagnostics can be split across
# source lines without embedding escape sequences.
panel_abort <- \(..., class_name = "panelio_error") {
  message_text <- paste0(c(...), collapse = "")
  stop(
    structure(
      class = c(class_name, "error", "condition"),
      list(message = message_text, call = NULL)
    )
  )
}

assert_scalar_text <- \(value, argument_name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    panel_abort(sprintf("Argument %s must be a non-empty string.", argument_name))
  }

  invisible(TRUE)
}

truncate_report <- \(values, report_limit) {
  value_text <- as.character(values)

  if (length(value_text) > report_limit) {
    return(paste0(
      paste(value_text[seq_len(report_limit)], collapse = ", "),
      sprintf(" ... (%d total)", length(value_text))
    ))
  }

  paste(value_text, collapse = ", ")
}

# ---- Value formatting and parsing -------------------------------------------
format_numeric_value <- \(value) {
  digit_template <- sprintf("%%.%dg", CSV_WRITE_DIGITS)
  short_form <- as.character(value)
  round_trip_value <- suppressWarnings(as.numeric(short_form))
  exact_match <- !is.na(value) & !is.na(round_trip_value) &
    (round_trip_value == value)
  needs_long_form <- !exact_match

  if (!any(needs_long_form)) {
    return(short_form)
  }

  long_form <- sprintf(digit_template, value)
  short_form[needs_long_form] <- long_form[needs_long_form]
  short_form
}

format_manifest_value <- \(value) {
  if (length(value) != 1L) {
    panel_abort("A manifest field must be a scalar value.")
  }

  if (is.numeric(value)) {
    return(format_numeric_value(value))
  }

  if (is.logical(value)) {
    return(as.character(value))
  }

  as.character(value)
}

coerce_manifest_values <- \(raw_values) {
  if (length(raw_values) == 0L) {
    return(character(0))
  }

  logical_tokens <- c(LOGICAL_TRUE_TOKENS, LOGICAL_FALSE_TOKENS)

  if (all(raw_values %in% logical_tokens)) {
    return(as.logical(raw_values))
  }

  if (all(grepl(NUMERIC_VALUE_PATTERN, raw_values))) {
    return(as.numeric(raw_values))
  }

  as.character(raw_values)
}

# ---- Long table structure ---------------------------------------------------
select_panel_columns <- \(long_table) {
  selected <- long_table[, LONG_TABLE_COLUMNS, drop = FALSE]
  rownames(selected) <- NULL
  selected
}

# Put a long table into a deterministic, input order independent ordering:
# observation, then time_identifier, then feature. Used to compare two long
# tables that are expected to hold the same panel content.
sort_panel_table <- \(long_table) {
  observation_levels <- sort(unique(as.character(long_table[[OBSERVATION_COLUMN]])))
  time_levels <- sort(unique(as.character(long_table[[TIME_IDENTIFIER_COLUMN]])))
  feature_levels <- sort(unique(as.character(long_table[[FEATURE_COLUMN]])))
  order_index <- order(
    match(long_table[[OBSERVATION_COLUMN]], observation_levels),
    match(long_table[[TIME_IDENTIFIER_COLUMN]], time_levels),
    match(long_table[[FEATURE_COLUMN]], feature_levels)
  )
  sorted_table <- long_table[order_index, , drop = FALSE]
  rownames(sorted_table) <- NULL
  sorted_table
}

assert_long_table <- \(x, context_name) {
  if (!is.data.frame(x)) {
    panel_abort(sprintf("%s requires a long-format data.frame.", context_name))
  }

  missing_columns <- setdiff(LONG_TABLE_COLUMNS, names(x))

  if (length(missing_columns) > 0) {
    panel_abort(
      sprintf(
        "%s is missing required columns: %s (expected %s).",
        context_name,
        paste(missing_columns, collapse = ", "),
        paste(LONG_TABLE_COLUMNS, collapse = ", ")
      )
    )
  }

  numeric_candidate_columns <- c(ELAPSED_TIME_COLUMN, ABUNDANCE_VALUE_COLUMN)
  non_numeric_columns <- numeric_candidate_columns[
    !vapply(x[numeric_candidate_columns], is.numeric, logical(1L))
  ]

  if (length(non_numeric_columns) > 0) {
    panel_abort(
      sprintf(
        "The following columns of %s must be numeric: %s.",
        context_name,
        paste(non_numeric_columns, collapse = ", ")
      )
    )
  }

  empty_columns <- LONG_TABLE_COLUMNS[
    vapply(x[LONG_TABLE_COLUMNS], \(column) any(is.na(column)), logical(1L))
  ]

  if (length(empty_columns) > 0) {
    panel_abort(
      sprintf(
        "The following columns of %s must not contain missing values: %s.",
        context_name,
        paste(empty_columns, collapse = ", ")
      )
    )
  }

  assert_label_grouping(x[[LABEL_COLUMN]])
  invisible(TRUE)
}

assert_unique_panel_keys <- \(long_table) {
  key_frame <- long_table[, c(
    OBSERVATION_COLUMN,
    TIME_IDENTIFIER_COLUMN,
    FEATURE_COLUMN
  ), drop = FALSE]
  duplicate_flag <- duplicated(key_frame) | duplicated(key_frame, fromLast = TRUE)

  if (!any(duplicate_flag)) {
    return(invisible(TRUE))
  }

  duplicate_rows <- key_frame[duplicate_flag, , drop = FALSE]
  duplicate_keys <- unique(
    do.call(paste, c(unname(as.list(duplicate_rows)), sep = " | "))
  )

  panel_abort(sprintf(
    "Duplicate composite key (observation, time_identifier, feature) in the long table: %s.",
    truncate_report(duplicate_keys, DUPLICATE_KEY_REPORT_COUNT)
  ))
}

assert_balanced_panel <- \(long_table, observation_levels, time_levels, feature_levels) {
  observation_count <- length(observation_levels)
  time_count <- length(time_levels)
  feature_count <- length(feature_levels)
  expected_cell_count <- observation_count * time_count * feature_count

  if (nrow(long_table) != expected_cell_count) {
    panel_abort(sprintf(
      paste0(
        "Unbalanced panel detected: expected %d x %d x %d = %d cells but found %d rows. ",
        "Only balanced longitudinal panels are supported."
      ),
      observation_count,
      time_count,
      feature_count,
      expected_cell_count,
      nrow(long_table)
    ))
  }

  observation_coverage <- table(
    as.character(long_table[[TIME_IDENTIFIER_COLUMN]]),
    as.character(long_table[[OBSERVATION_COLUMN]])
  )
  uncovered_time_layer <- rowSums(observation_coverage == 0L) > 0L
  uncovered_observation <- colSums(observation_coverage == 0L) > 0L

  if (any(uncovered_time_layer) || any(uncovered_observation)) {
    panel_abort(sprintf(
      paste0(
        "Unbalanced panel detected: every time layer must share the same observation_id set; ",
        "%d time layers miss observations and %d observations miss time layers."
      ),
      sum(uncovered_time_layer),
      sum(uncovered_observation)
    ))
  }

  feature_coverage <- table(
    as.character(long_table[[OBSERVATION_COLUMN]]),
    as.character(long_table[[FEATURE_COLUMN]])
  )

  if (any(feature_coverage == 0L)) {
    panel_abort(sprintf(
      "Unbalanced panel detected: %d observation x feature combinations are missing.",
      sum(feature_coverage == 0L)
    ))
  }

  # Every observation in every time layer carries one sparse composition.
  composition_sum <- tapply(
    as.numeric(long_table[[ABUNDANCE_VALUE_COLUMN]]),
    list(
      as.character(long_table[[OBSERVATION_COLUMN]]),
      as.character(long_table[[TIME_IDENTIFIER_COLUMN]])
    ),
    sum
  )
  assert_compositional_sum(composition_sum, "observation x time_identifier")
  assert_single_breakpoint(
    unique_time_breakpoint(long_table),
    unique_time_elapsed(long_table)
  )

  invisible(TRUE)
}

unique_time_breakpoint <- \(long_table) {
  first_flag <- !duplicated(long_table[[TIME_IDENTIFIER_COLUMN]])
  as.logical(long_table[[IS_INTERVENTION_COLUMN]][first_flag])
}

unique_time_elapsed <- \(long_table) {
  first_flag <- !duplicated(long_table[[TIME_IDENTIFIER_COLUMN]])
  as.numeric(long_table[[ELAPSED_TIME_COLUMN]][first_flag])
}

# ---- Metadata table validation ----------------------------------------------
assert_meta_table <- \(meta_table, required_columns, meta_name) {
  if (!is.data.frame(meta_table)) {
    panel_abort(sprintf("%s must be a data.frame.", meta_name))
  }

  missing_columns <- setdiff(required_columns, names(meta_table))

  if (length(missing_columns) > 0) {
    panel_abort(sprintf(
      "%s is missing required columns: %s.",
      meta_name,
      paste(missing_columns, collapse = ", ")
    ))
  }

  if (nrow(meta_table) == 0L) {
    panel_abort(sprintf("%s must not be an empty table.", meta_name))
  }

  has_missing_value <- any(vapply(
    meta_table[required_columns],
    \(column) any(is.na(column)),
    logical(1L)
  ))

  if (has_missing_value) {
    panel_abort(sprintf("%s must not contain missing values.", meta_name))
  }

  invisible(TRUE)
}

assert_observation_meta <- \(observation_meta) {
  assert_meta_table(observation_meta, OBSERVATION_META_COLUMNS, "observation_meta")

  if (anyDuplicated(observation_meta[[OBSERVATION_ID_COLUMN]])) {
    panel_abort("observation_id must be unique in observation_meta.")
  }

  assert_label_grouping(observation_meta[[LABEL_COLUMN]])
  invisible(TRUE)
}

# label is the grouping variable of the panel, so it travels as an R factor and
# the number of distinct groups stays small.
assert_label_grouping <- \(label_value) {
  if (!is.factor(label_value)) {
    panel_abort(
      "label is the grouping variable and must be an R factor, but a ",
      paste(class(label_value), collapse = "/"),
      " was supplied."
    )
  }

  if (anyNA(label_value)) {
    panel_abort("label must not contain missing values.")
  }

  group_count <- nlevels(droplevels(label_value))

  if (group_count > MAX_LABEL_GROUP_COUNT) {
    panel_abort(sprintf(
      "label is the grouping variable, so at most %d groups are allowed but %d were found: %s.",
      MAX_LABEL_GROUP_COUNT,
      group_count,
      truncate_report(levels(droplevels(label_value)), UNKNOWN_ID_REPORT_COUNT)
    ))
  }

  invisible(TRUE)
}

# A panel carries at most one intervention breakpoint. Along the elapsed_time
# order the flag must be FALSE up front and, once it turns TRUE, it may never
# return to FALSE.
assert_single_breakpoint <- \(is_intervention_value, elapsed_time_value) {
  time_order <- order(elapsed_time_value)
  ordered_flag <- as.logical(is_intervention_value)[time_order]
  transition_count <- sum(ordered_flag[-1L] != ordered_flag[-length(ordered_flag)])

  if (ordered_flag[1L]) {
    panel_abort(
      "is_intervention must be FALSE for the earliest time layer; the panel is ",
      "expected to start before the intervention."
    )
  }

  if (any(ordered_flag[-1L] < ordered_flag[-length(ordered_flag)])) {
    panel_abort(
      "is_intervention must stay TRUE once it turns TRUE; a panel may hold at ",
      "most one intervention breakpoint, so FALSE may never follow TRUE."
    )
  }

  if (transition_count > MAX_INTERVENTION_BREAKPOINTS) {
    panel_abort(sprintf(
      paste0(
        "A panel may hold at most %d intervention breakpoint but %d FALSE to ",
        "TRUE transitions were found."
      ),
      MAX_INTERVENTION_BREAKPOINTS,
      transition_count
    ))
  }

  invisible(TRUE)
}

assert_time_meta <- \(time_meta) {
  assert_meta_table(time_meta, TIME_META_COLUMNS, "time_meta")

  if (anyDuplicated(time_meta[[TIME_IDENTIFIER_COLUMN]])) {
    panel_abort("time_identifier must be unique in time_meta.")
  }

  if (!is.numeric(time_meta[[ELAPSED_TIME_COLUMN]])) {
    panel_abort("elapsed_time must be numeric in time_meta.")
  }

  if (!is.numeric(time_meta[[DELTA_T_COLUMN]])) {
    panel_abort("delta_t must be numeric in time_meta.")
  }

  if (!is.logical(time_meta[[IS_INTERVENTION_META_COLUMN]])) {
    panel_abort("is_intervention_timepoint must be logical in time_meta.")
  }

  assert_delta_t_consistency(
    time_meta[[ELAPSED_TIME_COLUMN]],
    time_meta[[DELTA_T_COLUMN]]
  )
  assert_single_breakpoint(
    time_meta[[IS_INTERVENTION_META_COLUMN]],
    time_meta[[ELAPSED_TIME_COLUMN]]
  )

  invisible(TRUE)
}

# delta_t is the interval since the previous time layer and may be non uniform.
# The earliest layer has no predecessor, so its delta_t must be 0; every later
# layer must match the elapsed_time difference to its predecessor.
assert_delta_t_consistency <- \(elapsed_time_value, delta_t_value) {
  time_order <- order(elapsed_time_value)
  ordered_elapsed <- as.numeric(elapsed_time_value)[time_order]
  ordered_delta <- as.numeric(delta_t_value)[time_order]
  expected_delta <- c(0, diff(ordered_elapsed))

  if (any(abs(ordered_delta - expected_delta) > DELTA_T_TOLERANCE)) {
    panel_abort(sprintf(
      paste0(
        "delta_t must equal the elapsed_time difference to the previous time ",
        "layer, with 0 for the earliest layer. Expected %s but found %s."
      ),
      paste(expected_delta, collapse = ", "),
      paste(ordered_delta, collapse = ", ")
    ))
  }

  invisible(TRUE)
}

assert_array_meta_alignment <- \(abundance_array, observation_meta, time_meta) {
  if (!is.array(abundance_array) || length(dim(abundance_array)) != 3L) {
    panel_abort(
      "abundance_array must be a three-dimensional array ",
      "[observation, feature, time_identifier]."
    )
  }

  assert_observation_meta(observation_meta)
  assert_time_meta(time_meta)

  array_dimnames <- dimnames(abundance_array)

  if (is.null(array_dimnames) || length(array_dimnames) != 3L) {
    panel_abort("abundance_array must carry complete dimnames.")
  }

  observation_dimnames <- array_dimnames[[1L]]
  feature_dimnames <- array_dimnames[[2L]]
  time_dimnames <- array_dimnames[[3L]]

  if (is.null(observation_dimnames) || is.null(feature_dimnames) || is.null(time_dimnames)) {
    panel_abort("None of the three abundance_array dimnames may be NULL.")
  }

  duplicated_dimnames <- c(
    observation_dimnames,
    feature_dimnames,
    time_dimnames
  )

  if (anyDuplicated(duplicated_dimnames)) {
    panel_abort("abundance_array dimnames must not contain duplicates.")
  }

  expected_dim <- c(
    length(observation_meta[[OBSERVATION_ID_COLUMN]]),
    length(feature_dimnames),
    length(time_meta[[TIME_IDENTIFIER_COLUMN]])
  )

  if (!identical(as.integer(dim(abundance_array)), as.integer(expected_dim))) {
    panel_abort(sprintf(
      "Array dimensions disagree with the metadata key lengths: array %s, metadata %s.",
      paste(dim(abundance_array), collapse = " x "),
      paste(expected_dim, collapse = " x ")
    ))
  }

  if (!identical(
    as.character(observation_dimnames),
    as.character(observation_meta[[OBSERVATION_ID_COLUMN]])
  )) {
    panel_abort("The first array dimension must align exactly with the observation_meta keys.")
  }

  if (!identical(
    as.character(time_dimnames),
    as.character(time_meta[[TIME_IDENTIFIER_COLUMN]])
  )) {
    panel_abort("The third array dimension must align exactly with the time_meta keys.")
  }

  identifier_pattern <- paste0(
    "^",
    TIME_FILE_PREFIX,
    TIME_IDENTIFIER_DIGIT_CLASS,
    "{",
    TIME_IDENTIFIER_DIGITS,
    "}$"
  )
  malformed_identifier <- grep(
    identifier_pattern,
    time_dimnames,
    invert = TRUE,
    value = TRUE
  )

  if (length(malformed_identifier) > 0L) {
    panel_abort(sprintf(
      "time_identifier must match %s followed by %d digits; offending values: %s.",
      TIME_FILE_PREFIX,
      TIME_IDENTIFIER_DIGITS,
      truncate_report(malformed_identifier, UNKNOWN_ID_REPORT_COUNT)
    ))
  }

  if (length(time_dimnames) > MAX_TIME_LAYER_COUNT) {
    panel_abort(sprintf(
      "Time layer count %d exceeds the configured maximum %d.",
      length(time_dimnames),
      MAX_TIME_LAYER_COUNT
    ))
  }

  if (anyNA(abundance_array)) {
    panel_abort(
      "abundance_array must not contain missing cells; an unbalanced panel was detected."
    )
  }

  # Every observation in every time layer holds one sparse composition, so the
  # abundance values must sum to the target along the feature dimension.
  assert_compositional_sum(
    as.vector(apply(abundance_array, c(1L, 3L), sum)),
    "observation x time_identifier"
  )

  invisible(TRUE)
}

# Each abundance matrix holds sparse compositional data, which means the values
# of one observation in one time layer must sum to the configured target.
assert_compositional_sum <- \(group_sum_value, group_label) {
  deviation <- abs(as.numeric(group_sum_value) - COMPOSITION_SUM_TARGET)
  offending <- deviation > COMPOSITION_TOLERANCE

  if (any(offending)) {
    panel_abort(sprintf(
      paste0(
        "Abundance values must be compositional: every %s group must sum to %s ",
        "within %s, but %d group(s) deviate, for example %s."
      ),
      group_label,
      format(COMPOSITION_SUM_TARGET, scientific = FALSE),
      format(COMPOSITION_TOLERANCE, scientific = FALSE),
      sum(offending),
      truncate_report(
        format(group_sum_value[offending], digits = 12L),
        UNKNOWN_ID_REPORT_COUNT
      )
    ))
  }

  invisible(TRUE)
}

# ---- Manifest grammar -------------------------------------------------------
strip_manifest_comments <- \(raw_lines) {
  comment_flag <- startsWith(raw_lines, COMMENT_PREFIX)
  trimmed_lines <- trimws(raw_lines)
  content_lines <- trimmed_lines[!comment_flag & nzchar(trimmed_lines)]
  content_lines
}

parse_manifest_pairs <- \(raw_lines) {
  assignment_lines <- grep(PAIR_PATTERN, raw_lines, fixed = TRUE, value = TRUE)

  if (length(assignment_lines) == 0L) {
    return(data.frame())
  }

  field_chunks <- strsplit(assignment_lines, FIELD_SEPARATOR, fixed = TRUE)
  flat_fields <- unlist(field_chunks, use.names = FALSE)
  line_index <- rep(seq_along(field_chunks), times = lengths(field_chunks))
  split_position <- regexpr(PAIR_PATTERN, flat_fields, fixed = TRUE)
  pair_keys <- trimws(substr(flat_fields, 1L, split_position - 1L))
  pair_values <- trimws(substr(flat_fields, split_position + 1L, nchar(flat_fields)))
  key_names <- unique(pair_keys)

  pair_matrix <- matrix(
    NA_character_,
    nrow = length(key_names),
    ncol = length(assignment_lines),
    dimnames = list(key_names, NULL)
  )
  pair_matrix[cbind(match(pair_keys, key_names), line_index)] <- pair_values
  parsed_columns <- lapply(seq_len(nrow(pair_matrix)), \(row_index) {
    coerce_manifest_values(pair_matrix[row_index, ])
  })
  names(parsed_columns) <- key_names

  as.data.frame(parsed_columns, stringsAsFactors = FALSE, optional = TRUE, check.names = FALSE)
}

# The manifest carries label as plain text, so the grouping column is turned
# back into an R factor as soon as the observation metadata is materialised. A
# missing label column is left alone so that the column validation downstream
# can report the actual problem.
coerce_observation_meta_frame <- \(observation_meta) {
  if (nrow(observation_meta) == 0L || !(LABEL_COLUMN %in% names(observation_meta))) {
    return(observation_meta)
  }

  observation_meta[[LABEL_COLUMN]] <- factor(as.character(observation_meta[[LABEL_COLUMN]]))
  observation_meta
}

parse_manifest_blocks <- \(raw_lines) {
  content_lines <- strip_manifest_comments(raw_lines)
  block_flag <- grepl(BLOCK_LINE_PATTERN, content_lines)
  block_names <- sub(BLOCK_LINE_PATTERN, "\\1", content_lines[block_flag])
  block_start <- which(block_flag)
  block_end <- c(block_start[-1L] - 1L, length(content_lines))

  if (length(block_start) == 0L) {
    panel_abort("manifest.txt parsing failed: no [block] marker was found.")
  }

  block_lines <- Map(
    \(start_line, end_line) {
      if (end_line < start_line + 1L) {
        return(character(0))
      }

      content_lines[seq.int(start_line + 1L, end_line)]
    },
    block_start,
    block_end
  )
  names(block_lines) <- block_names
  block_lines
}

parse_manifest <- \(raw_lines) {
  required_blocks <- c(
    OBSERVATION_META_BLOCK,
    TIME_META_BLOCK,
    FILE_INDEX_BLOCK
  )
  block_lines <- parse_manifest_blocks(raw_lines)

  if (anyNA(match(required_blocks, names(block_lines)))) {
    missing_blocks <- required_blocks[is.na(match(required_blocks, names(block_lines)))]
    panel_abort(sprintf(
      "manifest.txt parsing failed: missing block %s.",
      paste(missing_blocks, collapse = ", ")
    ))
  }

  duplicated_blocks <- unique(names(block_lines)[duplicated(names(block_lines))])

  if (length(duplicated_blocks) > 0L) {
    panel_abort(sprintf(
      "manifest.txt parsing failed: block %s is defined more than once.",
      paste(duplicated_blocks, collapse = ", ")
    ))
  }

  observation_meta <- coerce_observation_meta_frame(
    parse_manifest_pairs(block_lines[[OBSERVATION_META_BLOCK]])
  )
  time_meta <- parse_manifest_pairs(block_lines[[TIME_META_BLOCK]])
  file_index <- trimws(block_lines[[FILE_INDEX_BLOCK]])

  if (nrow(observation_meta) == 0L) {
    panel_abort("manifest.txt parsing failed: the observation_meta block is empty.")
  }

  if (nrow(time_meta) == 0L) {
    panel_abort("manifest.txt parsing failed: the time_meta block is empty.")
  }

  if (length(file_index) == 0L) {
    panel_abort("manifest.txt parsing failed: the file_index block is empty.")
  }

  list(
    observation_meta = observation_meta,
    time_meta = time_meta,
    file_index = file_index
  )
}

# ---- Metadata attachment ----------------------------------------------------
attach_metadata <- \(long_table, observation_meta, time_meta) {
  observation_match <- match(
    long_table[[OBSERVATION_COLUMN]],
    observation_meta[[OBSERVATION_ID_COLUMN]]
  )
  time_match <- match(
    long_table[[TIME_IDENTIFIER_COLUMN]],
    time_meta[[TIME_IDENTIFIER_COLUMN]]
  )

  if (anyNA(observation_match)) {
    unknown_observation <- unique(
      long_table[[OBSERVATION_COLUMN]][is.na(observation_match)]
    )
    panel_abort(sprintf(
      "Observations that are not registered in observation_meta: %s.",
      truncate_report(unknown_observation, UNKNOWN_ID_REPORT_COUNT)
    ))
  }

  if (anyNA(time_match)) {
    unknown_time <- unique(long_table[[TIME_IDENTIFIER_COLUMN]][is.na(time_match)])
    panel_abort(sprintf(
      "Time layers that are not registered in time_meta: %s.",
      truncate_report(unknown_time, UNKNOWN_ID_REPORT_COUNT)
    ))
  }

  long_table[[ELAPSED_TIME_COLUMN]] <- time_meta[[ELAPSED_TIME_COLUMN]][time_match]
  long_table[[LABEL_COLUMN]] <- observation_meta[[LABEL_COLUMN]][observation_match]
  long_table[[IS_INTERVENTION_COLUMN]] <- time_meta[[IS_INTERVENTION_META_COLUMN]][time_match]
  select_panel_columns(long_table)
}

# ---- Structural conversions -------------------------------------------------
panel_wide_to_long <- \(wide_table, time_identifier) {
  feature_names <- names(wide_table)[-1L]
  abundance_matrix <- as.matrix(wide_table[, feature_names, drop = FALSE])
  storage.mode(abundance_matrix) <- "double"
  observation_count <- nrow(abundance_matrix)
  feature_count <- length(feature_names)
  observation_values <- as.character(wide_table[[1L]])

  # Transposing once makes the base R column major layout emit the canonical
  # order inside one time layer: observation slowest, feature fastest.
  canonical_slice <- t(abundance_matrix)

  data.frame(
    observation = rep(observation_values, each = feature_count),
    time_identifier = rep(
      time_identifier,
      times = observation_count * feature_count
    ),
    feature = rep(feature_names, times = observation_count),
    abundance_value = as.vector(canonical_slice),
    stringsAsFactors = FALSE
  )
}

panel_array_to_long <- \(abundance_array) {
  array_dim <- dim(abundance_array)
  observation_count <- array_dim[1L]
  feature_count <- array_dim[2L]
  time_count <- array_dim[3L]

  # The stored array is [observation, feature, time_identifier] in base R
  # column major order, so it is permuted to [feature, observation,
  # time_identifier] first. The emitted rows are then in the canonical order
  # time_identifier x observation x feature with feature varying fastest,
  # which is exactly the order produced by concatenating the archive layers.
  canonical_array <- aperm(abundance_array, c(2L, 1L, 3L))
  canonical_dimnames <- dimnames(canonical_array)

  data.frame(
    observation = rep(
      canonical_dimnames[[2L]],
      each = feature_count,
      times = time_count
    ),
    time_identifier = rep(
      canonical_dimnames[[3L]],
      each = feature_count * observation_count
    ),
    feature = rep(
      canonical_dimnames[[1L]],
      times = observation_count * time_count
    ),
    abundance_value = as.vector(canonical_array),
    stringsAsFactors = FALSE
  )
}

long_to_panel_array <- \(long_table, observation_levels, time_levels, feature_levels) {
  abundance_array <- array(
    NA_real_,
    dim = c(
      length(observation_levels),
      length(feature_levels),
      length(time_levels)
    ),
    dimnames = list(observation_levels, feature_levels, time_levels)
  )
  cell_index <- cbind(
    match(long_table[[OBSERVATION_COLUMN]], observation_levels),
    match(long_table[[FEATURE_COLUMN]], feature_levels),
    match(long_table[[TIME_IDENTIFIER_COLUMN]], time_levels)
  )
  abundance_array[cell_index] <- as.double(long_table[[ABUNDANCE_VALUE_COLUMN]])
  abundance_array
}

# ---- Time layer file names --------------------------------------------------
# delta_t is the interval since the previous time layer and is allowed to be non
# uniform. The earliest layer has no predecessor, so its delta_t is 0.
# delta_t is derived from elapsed_time in elapsed_time order, not in the row
# order of the incoming table, so a shuffled long table still yields the same
# result as a sorted one.
derive_delta_t <- \(elapsed_time_value) {
  elapsed_numeric <- as.numeric(elapsed_time_value)
  time_order <- order(elapsed_numeric)
  ordered_delta <- c(0, diff(elapsed_numeric[time_order]))
  ordered_delta[order(time_order)]
}

build_time_identifier <- \(time_index) {
  sprintf(
    "%s%0*d",
    TIME_FILE_PREFIX,
    TIME_IDENTIFIER_DIGITS,
    as.integer(time_index)
  )
}

build_time_identifier_sequence <- \(time_count) {
  if (time_count > MAX_TIME_LAYER_COUNT) {
    panel_abort(sprintf(
      "Time layer count %d exceeds the configured maximum %d.",
      time_count,
      MAX_TIME_LAYER_COUNT
    ))
  }

  vapply(seq_len(time_count), build_time_identifier, character(1L), USE.NAMES = FALSE)
}

time_file_name <- \(time_identifier) {
  paste0(time_identifier, TIME_FILE_SUFFIX)
}

# ---- csv helpers ------------------------------------------------------------
read_csv_table <- \(file_path) {
  if (!file.exists(file_path)) {
    panel_abort(sprintf("Archive is missing the file: %s", basename(file_path)))
  }

  utils::read.csv(
    file = file_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = CSV_FILE_ENCODING
  )
}

# write.table formats doubles with 15 significant digits, which is not enough to
# reproduce a double exactly. 17 significant digits do, so every numeric cell is
# written explicitly. A whole number also gets a ".0" suffix, otherwise
# type.convert would read the column back as integer instead of double.
format_csv_numeric <- \(numeric_value) {
  formatted <- sprintf(paste0("%.", CSV_WRITE_DIGITS, "g"), numeric_value)
  plain_integer <- is.finite(numeric_value) & !grepl("[.eE]", formatted)
  formatted[plain_integer] <- paste0(formatted[plain_integer], ".0")
  formatted
}

write_csv_table <- \(table_data, file_path) {
  numeric_columns <- vapply(table_data, is.numeric, logical(1L))
  csv_data <- table_data
  csv_data[numeric_columns] <- lapply(
    csv_data[numeric_columns],
    format_csv_numeric
  )
  utils::write.csv(
    x = csv_data,
    file = file_path,
    row.names = FALSE,
    na = "",
    fileEncoding = CSV_FILE_ENCODING
  )
  invisible(file_path)
}

# ---- External tar helpers (never call setwd) --------------------------------
assert_tar_available <- \() {
  tar_is_available <- nzchar(Sys.which(TAR_COMMAND))

  if (!tar_is_available) {
    panel_abort(sprintf(
      "External command '%s' was not found; archive creation and extraction are unavailable.",
      TAR_COMMAND
    ))
  }

  invisible(TRUE)
}

create_archive <- \(archive_path, staging_dir) {
  assert_tar_available()
  unlink(archive_path)
  tar_arguments <- c(
    TAR_CREATE_FLAGS,
    archive_path,
    "-C",
    staging_dir,
    TAR_STAGING_ENTRY
  )
  tar_output <- suppressWarnings(system2(
    TAR_COMMAND,
    tar_arguments,
    stdout = TRUE,
    stderr = TRUE
  ))

  if (!file.exists(archive_path)) {
    panel_abort(sprintf(
      "Archive creation failed: %s\n%s",
      archive_path,
      paste(tar_output, collapse = "\n")
    ))
  }

  normalizePath(archive_path, mustWork = TRUE)
}

extract_archive <- \(archive_path, target_dir) {
  assert_tar_available()
  tar_arguments <- c(
    TAR_EXTRACT_FLAGS,
    archive_path,
    "-C",
    target_dir
  )
  tar_output <- suppressWarnings(system2(
    TAR_COMMAND,
    tar_arguments,
    stdout = TRUE,
    stderr = TRUE
  ))

  if (!file.exists(target_dir)) {
    panel_abort(sprintf("Archive extraction directory is missing: %s", target_dir))
  }

  invisible(tar_output)
}
