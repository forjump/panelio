# 01_probe_metadata.R
# Probe one archive per dataset: dims, label levels, time_meta (t0 / flags),
# feature id naming, and cross-check against truth_table.csv.
# Vectorized only (lapply / vapply / apply); no for/while loops.
suppressMessages(pkgload::load_all("/mnt/c/Users/Zhang/Desktop/clsc_r", quiet = TRUE))

bench_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/panelio"
truth_path <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/truth/truth_table.csv"
out_path <- "/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results/dataset_metadata.csv"

truth <- read.csv(truth_path, stringsAsFactors = FALSE, check.names = FALSE)
datasets <- sprintf("A%d", 1:12)

probe_one <- function(ds) {
  dpath <- file.path(bench_dir, ds)
  archives <- sort(list.files(dpath, pattern = "\\.tar\\.gz$"))
  lt <- read_panel_table(file.path(dpath, archives[1L]))
  pa <- as_panel_array(lt)
  ab <- pa$abundance_array
  om <- pa$observation_meta
  tm <- pa$time_meta
  tr <- truth[truth$dataset == ds, , drop = FALSE]
  list(
    dataset = ds,
    n_archives = length(archives),
    first_archive = archives[1L],
    dims = paste(dim(ab), collapse = "x"),
    obs_meta_cols = paste(names(om), collapse = ","),
    n_obs = nrow(om),
    label_values = paste(unique(as.character(om$label)), collapse = ","),
    label_levels = if (is.factor(om$label)) paste(levels(om$label), collapse = ",") else "NA",
    time_meta_cols = paste(names(tm), collapse = ","),
    n_time = nrow(tm),
    time_ids = paste(as.character(tm$time_identifier), collapse = ","),
    intervention_flags = paste(as.character(tm$is_intervention_timepoint), collapse = ","),
    feature_ids_head = paste(head(dimnames(ab)[[2L]], 10L), collapse = ","),
    feature_ids_tail = paste(tail(dimnames(ab)[[2L]], 3L), collapse = ","),
    truth_t0 = tr$t0,
    truth_T = tr$T,
    truth_effect_size = tr$effect_size,
    truth_effect_taxa = tr$effect_taxa,
    truth_effect_signs = tr$effect_signs
  )
}

rows <- lapply(datasets, probe_one)
meta_df <- as.data.frame(do.call(rbind, lapply(rows, unlist)), stringsAsFactors = FALSE)
write.csv(meta_df, out_path, row.names = FALSE)
cat("WROTE:", out_path, "\n")
print(meta_df)
