# 00_check_env.R
# Check the R environment: pkgload availability, package load strategy, src dir.
cat("R version:", R.version.string, "\n")
cat("pkgload available:", requireNamespace("pkgload", quietly = TRUE), "\n")
src_dir <- "/mnt/c/Users/Zhang/Desktop/clsc_r/src"
cat("src files:\n")
print(list.files(src_dir))
so_files <- list.files(src_dir, pattern = "\\.so$")
cat("compiled .so present:", length(so_files) > 0L, "\n")
if (length(so_files) > 0L) print(so_files)
cat("benchmark/results dir:", dir.exists("/mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/results"), "\n")
