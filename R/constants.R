# =============================================================================
# R/constants.R
# Every tunable value of the package lives here as a named, documented
# constant. There is no external config file and no yaml dependency: all
# "magic numbers" are replaced by these named constants, so the package
# depends only on base R (utils). Keeping every value in one file makes the
# tuning surface explicit and auditable.
# =============================================================================

# ---- archive layout ---------------------------------------------------------
TIME_FILE_PREFIX <- "time_"
TIME_FILE_SUFFIX <- ".csv"
MANIFEST_FILE_NAME <- "manifest.txt"
ARCHIVE_FILE_TEMPLATE <- "%s_%s.tar.gz"

# ---- csv / encoding ---------------------------------------------------------
# 17 significant digits uniquely determine an IEEE-754 double, so values such
# as 1/3, 1e-20, 1e300 or pi/4 survive a write/read round trip bit-for-bit.
CSV_FILE_ENCODING <- "UTF-8"
CSV_WRITE_DIGITS <- 17L

# ---- external tar invocation (no setwd ever used) --------------------------
TAR_COMMAND <- "tar"
TAR_CREATE_FLAGS <- "-czf"
TAR_EXTRACT_FLAGS <- "-xzf"
TAR_STAGING_ENTRY <- "."

# ---- manifest grammar -------------------------------------------------------
COMMENT_PREFIX <- "#"
BLOCK_LINE_PATTERN <- "^\\[([[:alnum:]_]+)\\]$"
KEY_VALUE_SEPARATOR <- "="
FIELD_SEPARATOR <- ";"
MANIFEST_FIELD_JOIN <- "; "
PAIR_PATTERN <- "="
NUMERIC_VALUE_PATTERN <- "^-?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
LOGICAL_TRUE_TOKENS <- "TRUE"
LOGICAL_FALSE_TOKENS <- "FALSE"
OBSERVATION_META_BLOCK <- "observation_meta"
TIME_META_BLOCK <- "time_meta"
FILE_INDEX_BLOCK <- "file_index"
MANIFEST_HEADER_COMMENT <- "manifest.txt machine-parsable"

# ---- time layer metadata ----------------------------------------------------
# elapsed_time is a generic cumulative time; its unit is defined by the data
# set itself and may be hours, days or weeks. delta_t is the interval since the
# previous time layer and may be non uniform; the earliest layer has delta_t 0.
DELTA_T_TOLERANCE <- 1e-8

# ---- intervention breakpoint ------------------------------------------------
# A panel may hold at most one intervention breakpoint. The flag must start
# FALSE and, once it turns TRUE, must never fall back to FALSE.
MAX_INTERVENTION_BREAKPOINTS <- 1L

# ---- grouping (observation label) -------------------------------------------
# label is the grouping variable, carried as an R factor; the number of distinct
# groups stays small (two groups being the common case).
MAX_LABEL_GROUP_COUNT <- 4L

# ---- compositional abundance ------------------------------------------------
# Every observation within every time layer holds a sparse composition, so the
# abundance values of one observation in one time layer must sum to the target.
COMPOSITION_SUM_TARGET <- 1
COMPOSITION_TOLERANCE <- 1e-6

# ---- limits and reporting ----------------------------------------------------
TIME_IDENTIFIER_DIGITS <- 3L
TIME_IDENTIFIER_DIGIT_CLASS <- "[0-9]"
MAX_TIME_LAYER_COUNT <- 1000L
DUPLICATE_KEY_REPORT_COUNT <- 10L
MISSING_CELL_REPORT_COUNT <- 10L
UNKNOWN_ID_REPORT_COUNT <- 10L

# ---- preprocessing: feature-axis operations ---------------------------------
# All preprocessing operates on the feature (species) axis, dimension 2 of the
# three-dimensional abundance array [observation, feature, time_identifier].
NORMALIZATION_DEFAULT_METHOD <- "clr"
NORMALIZE_NONE   <- "none"
NORMALIZE_CLR    <- "clr"
NORMALIZE_TSS    <- "tss"
NORMALIZE_RLE    <- "rle"
NORMALIZE_CSS    <- "css"
NORMALIZE_LOG1P  <- "log1p"
NORMALIZE_ILR    <- "ilr"
NORMALIZE_RCLR   <- "rclr"
NORMALIZE_Z_SCORE <- "z_score"
NORMALIZE_MIN_MAX <- "min_max"
NORMALIZE_METHODS <- c(NORMALIZE_NONE, NORMALIZE_CLR, NORMALIZE_TSS,
                       NORMALIZE_RLE, NORMALIZE_CSS, NORMALIZE_LOG1P,
                       NORMALIZE_ILR, NORMALIZE_RCLR,
                       NORMALIZE_Z_SCORE, NORMALIZE_MIN_MAX)
# Normalization methods that require strictly positive abundances and therefore
# trigger the automatic zero-fill step unless zero_method = "none". rclr, z_score
# and min_max handle zeros themselves and never trigger the automatic fill.
NORMALIZE_REQUIRES_POSITIVE <- c(NORMALIZE_CLR, NORMALIZE_RLE, NORMALIZE_ILR)

ZERO_FILL_DEFAULT_METHOD <- "half_min"
ZERO_NONE           <- "none"
ZERO_CONSTANT       <- "constant"
ZERO_MULTIPLICATIVE <- "multiplicative"
ZERO_HALF_MIN       <- "half_min"
ZERO_SQRT_MIN       <- "sqrt_min"
ZERO_BAYESIAN       <- "bayesian"
ZERO_FILL_METHODS <- c(ZERO_NONE, ZERO_CONSTANT, ZERO_MULTIPLICATIVE,
                       ZERO_HALF_MIN, ZERO_SQRT_MIN, ZERO_BAYESIAN)
ZERO_FILL_CONSTANT <- 0.5      # fixed pseudo-count for the "constant" method
ZERO_FILL_DELTA    <- 1e-4     # replacement value for the multiplicative method
ZERO_FILL_FRACTION <- 0.5      # replacement fraction for the half_min / sqrt_min methods

FILTER_PREVALENCE_THRESHOLD <- 0.1  # default prevalence floor when the filter is active
FILTER_MIN_ABUNDANCE        <- 0    # a feature is "present" in a sample above this value

CSS_DEFAULT_QUANTILE <- 0.5   # cumulative-sum scaling quantile

# ---- smoothing: time-axis operations ----------------------------------------
# All smoothing operates on the time axis, dimension 3 of the abundance array
# [observation, feature, time_identifier]. The iteration lives in the Fortran
# kernels under src/; the R layer only dispatches one .Fortran call per method.
SMOOTH_DEFAULT_METHOD <- "ekf"
SMOOTH_EKF    <- "ekf"
SMOOTH_RTS    <- "rts"
SMOOTH_GAUSSIAN <- "gaussian"
SMOOTH_LOESS  <- "loess"
SMOOTH_SAVGOL <- "savgol"
SMOOTH_MA     <- "moving_average"
SMOOTH_EWMA   <- "ewma"
SMOOTH_SPLINE <- "spline"
SMOOTH_METHODS <- c(SMOOTH_EKF, SMOOTH_RTS, SMOOTH_GAUSSIAN, SMOOTH_LOESS,
                    SMOOTH_SAVGOL, SMOOTH_MA, SMOOTH_EWMA, SMOOTH_SPLINE)
# The Kalman family (EKF / RTS) smooths the whole series and intentionally
# ignores breakpoints; the other methods honour breakpoints so an intervention
# effect at the breakpoint is never smoothed away.
SMOOTH_KALMAN_METHODS <- c(SMOOTH_EKF, SMOOTH_RTS)

SMOOTH_KALMAN_Q_RATIO <- 0.1      # Q = qratio * R (process vs observation noise)
SMOOTH_KALMAN_INITVAR_SCALE <- 100   # P0 = initvar * R
SMOOTH_GAUSSIAN_SIGMA <- 2.0      # Gaussian kernel bandwidth (time-index units)
SMOOTH_LOESS_SPAN <- 0.75         # tricube loess span (fraction of segment)
SMOOTH_LOESS_DEGREE <- 1L
SMOOTH_SAVGOL_HALFWINDOW <- 2L
SMOOTH_SAVGOL_DEGREE <- 2L
SMOOTH_MA_HALFWINDOW <- 2L
SMOOTH_EWMA_ALPHA <- 0.3
SMOOTH_SPLINE_LAMBDA <- 1.0       # second-difference penalty strength

# ---- smoothing: paper-depth hierarchical NB-EKF / RTS -------------------------
# The Kalman family can run either on a Gaussian observation (scalar local level,
# the historical default) or on a negative-binomial observation model over raw
# counts, matching the CLSC hierarchical EKF-RTS of the paper: the latent state
# is the log relative abundance, the expected count is library size times the
# softmax of the state, and the measurement variance follows the delta method
# (mean + mean^2 / dispersion). A hierarchical EM step (shared per-feature
# innovation variance and dispersion across all observations) and optional
# scaling of the process noise by the sampling interval complete the paper depth.
SMOOTH_OBSERVATION_GAUSSIAN <- "gaussian"
SMOOTH_OBSERVATION_NB       <- "nb"
SMOOTH_OBSERVATION_DEFAULT  <- "gaussian"
SMOOTH_NB_DEFAULT_DISPERSION <- 1.0   # per-feature NB dispersion when not estimated
SMOOTH_NB_EM_ITERATIONS     <- 80L    # paper uses ~80 EM iterations
SMOOTH_NB_LOG_FLOOR         <- -16.0  # log-relative-abundance floor for zero counts
SMOOTH_NB_HALF_COUNT        <- 0.5    # pseudo-count for the log-relative initial state
SMOOTH_NB_MIN_DISPERSION    <- 1e-4   # lower clip of the EM dispersion estimate
SMOOTH_NB_MIN_VARIANCE      <- 1e-6   # lower clip of the EM innovation variance
SMOOTH_NB_EM_TOLERANCE      <- 1e-6   # EM convergence tolerance on (q, phi) change
SMOOTH_NB_PSEUDO_COUNT      <- 0.5    # pseudo-count used in the fitted-mean variance

# ---- causal estimation --------------------------------------------------------
# Per-feature causal effect estimation from the panel. Treatment assignment is
# the label column; the intervention onset t0 is the first intervention
# timepoint. All estimators report per-feature point estimates, standard
# errors, t-statistics, p-values and Benjamini-Hochberg FDR; time-varying
# (event-time / per-time) effects are reported separately.
CAUSAL_DEFAULT_METHOD <- "did"
CAUSAL_DID      <- "did"
CAUSAL_POOLED   <- "pooled_ols"
CAUSAL_EVENT    <- "event_study"
CAUSAL_SCM      <- "scm"
CAUSAL_METHODS  <- c(CAUSAL_DID, CAUSAL_POOLED, CAUSAL_EVENT, CAUSAL_SCM)

# Event study: relative time is measured against t0 (reltime = t - t0); the
# reference (omitted) period is the last pre-intervention time, reltime = -1.
CAUSAL_EVENT_BASE_REL <- -1L
# Ridge added to the event-study design Gram matrix for numerical stability.
CAUSAL_GRAM_RIDGE <- 1e-8
# Two-way fixed-effect within-demeaning (Hoethen / Prais-Winsten) iterations.
CAUSAL_MAX_DEMEAN_ITER <- 500L
CAUSAL_DEMEAN_TOLERANCE <- 1e-10
# Lawson-Hanson NNLS internal tolerance for synthetic control weights.
CAUSAL_NNLS_TOLERANCE <- 1e-12
# Two-sided p-value tail alpha used to summarise the FDR share in print().
CAUSAL_TAIL_ALPHA <- 0.05

# ---- causal: scm constraint (exactly-convex) ---------------------------------
# The synthetic-control estimator supports two weight constraints: the classic
# non-negative least squares (w >= 0, the historical default) and the paper's
# exactly-convex simplex (w >= 0 and sum(w) = 1), which keeps the counterfactual
# inside the convex hull of the donor trajectories. A scalar level offset
# (centred log-ratio translation) can absorb baseline differences between the
# treated unit and the convex donor mixture.
CAUSAL_SCM_NNLS    <- "nnls"
CAUSAL_SCM_SIMPLEX <- "simplex"
CAUSAL_SCM_CONSTRAINT_DEFAULT <- CAUSAL_SCM_NNLS
CAUSAL_SCM_OFFSET_DEFAULT     <- FALSE
CAUSAL_SCM_MAX_ITER           <- 4000L   # active-set KKT iterations for the simplex solver
CAUSAL_SCM_TOLERANCE          <- 1e-12   # KKT convergence tolerance of the simplex solver

# ---- causal: placebo / bootstrap / pre-fit diagnostics ----------------------
# The paper assesses significance with a placebo (permutation) procedure (every
# donor acting as a pseudo-treated unit, the remaining donors as the pseudo-pool,
# followed by a two-sample permutation test per taxon and BH-FDR across taxa),
# quantifies uncertainty with subject-level bootstrap percentile intervals, and
# summarises counterfactual reliability with the pre-period RMSE and its
# signal-to-noise ratio (Abadie-style diagnostic).
CAUSAL_PLACEBO_DEFAULT   <- FALSE
CAUSAL_BOOTSTRAP_DEFAULT <- FALSE
CAUSAL_PRE_FIT_DEFAULT   <- FALSE
CAUSAL_BOOTSTRAP_SAMPLES <- 500L     # subject-level bootstrap resamples
CAUSAL_BOOTSTRAP_LEVEL   <- 0.95     # percentile confidence level
CAUSAL_PERMUTATION_SAMPLES <- 999L   # relabelings of the two-sample permutation test
CAUSAL_SNR_EPS           <- 1e-10    # floor under the pre-fit RMSE for the SNR ratio

# Feature-level effect columns always reported, in this order.
CAUSAL_EFFECT_COLUMNS <- c(
  "feature", "method", "estimate", "standard_error", "statistic",
  "p_value", "fdr", "n_treated", "n_control",
  "mean_pre_treated", "mean_post_treated", "mean_pre_control", "mean_post_control",
  "effect_range"
)
