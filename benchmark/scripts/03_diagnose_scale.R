# 03_diagnose_scale.R
# Diagnose the scale issue: what does the NB-smoothed output look like, and
# what scale do estimates take under different input pipelines?
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth <- read.csv("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

arch <- file.path(bench_dir, "A1", "A1_rep01_v1.tar.gz")
lt <- read_panel_table(arch)
pa <- as_panel_array(lt)
ab <- pa$abundance_array
tr1 <- truth[truth$dataset == "A1", , drop = FALSE]
t0 <- as.integer(tr1$t0)
eff <- sprintf("taxa_%03d", as.integer(strsplit(tr1$effect_taxa, ";")[[1L]]))
signs <- as.integer(strsplit(tr1$effect_signs, ";")[[1L]])

cat("=== raw counts distribution ===\n")
cat("range:", range(ab), " zero prop:", mean(ab == 0), " libsize range:",
    range(apply(ab, c(1, 3), sum)), "\n")

cat("\n=== smoothed (rts, nb, hierarchical) output distribution ===\n")
pp <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                       zero = NULL, normalize = "none")
sp <- smooth_array(pp$abundance_array, method = "rts", observation_model = "nb",
                   hierarchical = TRUE, em_iterations = 80, time_scaling = TRUE,
                   time_meta = pa$time_meta)
cat("range:", range(sp), " mean:", mean(sp), " sd:", sd(sp), "\n")
cat("per-feature sd head:", head(round(apply(sp, 2, sd), 6)), "\n")
# treated vs control means per feature (post period)
trt <- pa$observation_meta$label == "treat"
post_idx <- t0:nrow(pa$time_meta)
m_tr_post <- apply(sp[trt, , post_idx, drop = FALSE], 2, mean)
m_co_post <- apply(sp[!trt, , post_idx, drop = FALSE], 2, mean)
cat("diff (treat-control) post mean for effect taxa:\n")
print(round(m_tr_post[eff] - m_co_post[eff], 6))
cat("effect taxa smoothed values (post, treat):\n")
print(round(sp[trt, eff, post_idx, drop = FALSE][1:6, , 1], 6))

cat("\n=== CLR-preprocessed (no smoothing) scm / did ===\n")
pp_clr <- preprocess_array(ab, filter = list(min_prevalence = 0.1, min_abundance = 0),
                           zero = list(method = "half_min"), normalize = "clr")
for (m in c("scm", "did", "pooled_ols", "event_study")) {
  est <- panel_causal_estimate(feature ~ label | observation + time,
                               data = pp_clr$abundance_array, method = m,
                               scm_constraint = "simplex", scm_offset = TRUE,
                               intervention_index = t0,
                               observation_meta = pa$observation_meta,
                               time_meta = pa$time_meta)
  e <- est$effects
  er <- e[e$feature %in% eff, , drop = FALSE]
  cat(sprintf("method %-11s: est on effect taxa: %s\n", m,
              paste(sprintf("%.3f", er$estimate), collapse = ", ")))
  cat(sprintf("  true signs: %s   est signs: %s\n",
              paste(signs, collapse = ","), paste(sign(er$estimate), collapse = ",")))
}

cat("\n=== CLR-preprocessed + gaussian-smooth scm ===\n")
sp_g <- smooth_array(pp_clr$abundance_array, method = "rts",
                     observation_model = "gaussian", time_meta = pa$time_meta)
est <- panel_causal_estimate(feature ~ label | observation + time,
                             data = sp_g, method = "scm",
                             scm_constraint = "simplex", scm_offset = TRUE,
                             intervention_index = t0,
                             observation_meta = pa$observation_meta,
                             time_meta = pa$time_meta)
er <- est$effects[est$effects$feature %in% eff, , drop = FALSE]
cat("scm(gauss-smooth) est:", paste(sprintf("%.3f", er$estimate), collapse = ", "), "\n")
cat("true signs:", paste(signs, collapse = ","), " est signs:",
    paste(sign(er$estimate), collapse = ","), "\n")
