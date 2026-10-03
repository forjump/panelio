# 09_diag_signal.R
# Does the CLR-preprocessed trajectory actually carry the +/-1.2 signal on A2?
# Inspect effect-taxon treated pre-post differences at each pipeline stage.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

tr <- truth[truth$dataset == "A2", , drop = FALSE]
t0 <- as.integer(tr$t0)
eff <- sprintf("taxa_%03d", as.integer(strsplit(tr$effect_taxa, ";")[[1L]]))
signs <- as.integer(strsplit(tr$effect_signs, ";")[[1L]])

arch <- file.path(bench_dir, "A2", "A2_rep01_v1.tar.gz")
pa <- as_panel_array(read_panel_table(arch))
ab <- pa$abundance_array
om <- pa$observation_meta
tm <- pa$time_meta
trt <- om$label == "treat"
pre <- seq_len(t0 - 1L)
post <- t0:nrow(tm)

cat("=== raw relab, effect taxa ===\n")
for (f in eff) {
  v_pre_tr <- mean(ab[trt, f, pre]); v_post_tr <- mean(ab[trt, f, post])
  v_pre_co <- mean(ab[!trt, f, pre]); v_post_co <- mean(ab[!trt, f, post])
  cat(sprintf("%s raw relab: tr pre->post %+.4f -> %+.4f | co pre->post %+.4f -> %+.4f | diff %+.4f\n",
              f, v_pre_tr, v_post_tr, v_pre_co, v_post_co,
              (v_post_tr - v_pre_tr) - (v_post_co - v_pre_co)))
}

cat("\n=== CLR (half_min) effect taxa ===\n")
pp <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                       zero = list(method = "half_min"), normalize = "clr")$abundance_array
for (f in eff) {
  d_tr <- mean(pp[trt, f, post]) - mean(pp[trt, f, pre])
  d_co <- mean(pp[!trt, f, post]) - mean(pp[!trt, f, pre])
  cat(sprintf("%s CLR: treat diff %+.3f | control diff %+.3f | DiD %+.3f (truth %+.1f)\n",
              f, d_tr, d_co, d_tr - d_co, 1.2 * signs[which(eff == f)]))
}

cat("\n=== rclr (zero-robust CLR) effect taxa ===\n")
pp_r <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = NULL, normalize = "rclr")$abundance_array
for (f in eff) {
  d_tr <- mean(pp_r[trt, f, post]) - mean(pp_r[trt, f, pre])
  d_co <- mean(pp_r[!trt, f, post]) - mean(pp_r[!trt, f, pre])
  cat(sprintf("%s rclr: treat diff %+.3f | control diff %+.3f | DiD %+.3f (truth %+.1f)\n",
              f, d_tr, d_co, d_tr - d_co, 1.2 * signs[which(eff == f)]))
}

cat("\n=== log1p (log relab) effect taxa ===\n")
pp_l <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                         zero = NULL, normalize = "log1p")$abundance_array
for (f in eff) {
  d_tr <- mean(pp_l[trt, f, post]) - mean(pp_l[trt, f, pre])
  d_co <- mean(pp_l[!trt, f, post]) - mean(pp_l[!trt, f, pre])
  cat(sprintf("%s log1p: treat diff %+.3f | control diff %+.3f | DiD %+.3f (truth %+.1f)\n",
              f, d_tr, d_co, d_tr - d_co, 1.2 * signs[which(eff == f)]))
}
