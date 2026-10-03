# =============================================================================
# R/panelio-package.R
# Package-level documentation and dynamic-library registration.
# =============================================================================

#' panelio: Microbiome Balanced Longitudinal Panel I/O
#'
#' Archive, read, write, preprocess and smooth balanced longitudinal microbiome
#' abundance panels. A dataset is stored as a gzipped tar archive holding one
#' wide table per time layer plus a machine parsable manifest carrying the
#' observation and time metadata.
#'
#' The abundance representation is the three-dimensional array
#' [observation, feature, time_identifier]. Preprocessing operates on the
#' feature axis; smoothing operates on the time axis via compiled Fortran
#' kernels.
#'
#' @useDynLib panelio, .registration = TRUE
#' @keywords internal
"_PACKAGE"
