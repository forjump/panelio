## ---------------------------------------------------------------------------
## make_figures.R
## Panelio paper figures (ggprism style, discrete offset axes, no grid lines).
##
## Run from WSL R (Ubuntu) on the Windows-mounted project:
##   wsl -d Ubuntu-26.04 bash -lc "Rscript '/mnt/c/Users/Zhang/Desktop/clsc_r/paper/figures/make_figures.R'"
##
## Requires: ggplot2 (>=3.5), ggprism, dplyr, tidyr, readr
## Outputs:  fig_power.pdf/.png, fig_type1.pdf/.png, fig_sign.pdf/.png, fig_real.pdf/.png
## ---------------------------------------------------------------------------

suppressMessages({library(ggplot2); library(ggprism); library(dplyr); library(tidyr); library(readr)})

## ---- paths (this script lives in paper/figures; results live two levels up) ----
cmdfa <- commandArgs(FALSE)
script_self <- sub("^--file=", "", cmdfa[grepl("^--file=", cmdfa)])
figd <- if (length(script_self) && nzchar(script_self[1])) dirname(normalizePath(script_self[1])) else "."
rm(cmdfa, script_self)
res <- normalizePath(file.path(figd, "..", "..", "benchmark", "results"), mustWork = TRUE)
if (!dir.exists(figd)) dir.create(figd, recursive = TRUE)

## ---- colour-blind friendly palette; SCM highlighted ----
pal <- c(did = "#000000", pooled_ols = "#E69F00", event_study = "#56B4E9", scm = "#D55E00")
mlab <- c(did = "DiD", pooled_ols = "Pooled OLS", event_study = "Event study", scm = "SCM")
mlev <- c("did", "pooled_ols", "event_study", "scm")

base_theme <- theme_prism(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.line = element_line(colour = "black")
  )
offset_xd <- scale_x_discrete(guide = "prism_offset")
offset_y <- function(...) scale_y_continuous(guide = "prism_offset", ...)

save2 <- function(p, stem, w, h) {
  pdf(file.path(figd, paste0(stem, ".pdf")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(figd, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
}

## ---- load simulation summary ----
sim <- read_csv(file.path(res, "summary_simulation.csv"), show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = paste0("A", 1:12)),
         method = factor(method, levels = mlev))

## ---- Fig 1: Power after FDR control ----
pwr <- sim %>% filter(method != "NA", !is.na(power_fdr))
p1 <- ggplot(pwr, aes(dataset, power_fdr, colour = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal, labels = mlab) +
  labs(x = "Simulation design", y = "Power (BH-FDR, effect features)") +
  base_theme + offset_xd +
  offset_y(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05))
save2(p1, "fig_power", 7.5, 4.6)

## ---- Fig 2: Type I error ----
ty <- sim %>% filter(!is.na(typeI_fdr))
p2 <- ggplot(ty, aes(dataset, typeI_fdr, colour = method, group = method)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal, labels = mlab) +
  labs(x = "Simulation design", y = "Type I error (non-effect features)") +
  base_theme + offset_xd +
  offset_y(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.04))
save2(p2, "fig_type1", 7.5, 4.6)

## ---- Fig 3: sign recovery ----
sg <- sim %>% filter(method != "NA", !is.na(sign_hit))
p3 <- ggplot(sg, aes(dataset, sign_hit, colour = method, group = method)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal, labels = mlab) +
  labs(x = "Simulation design", y = "Sign-hit rate (effect features)") +
  base_theme + offset_xd +
  offset_y(limits = c(0.3, 0.7), breaks = seq(0.3, 0.7, 0.05))
save2(p3, "fig_sign", 7.5, 4.6)

## ---- Fig 4: real-data SCM forest ----
shortg <- function(feat) {
  g <- sub(".*\\|g__", "", feat)
  ifelse(g == feat, sub(".*\\|o__", "", feat), g)
}
rd_f <- function(fn) read_csv(file.path(res, fn), show_col_types = FALSE) %>%
  mutate(genus = shortg(feature),
         lo = estimate - 1.96 * standard_error,
         hi = estimate + 1.96 * standard_error,
         sig = fdr < 0.05) %>%
  filter(sig)

## ---- Fig 4: real-data SCM forest; panel A (DIABIMMUNE, top significant genera) and
##      panel B (iHMP-IBD, all significant genera) placed side by side ----
dh <- rd_f("real_diabimmune_scm.csv") %>%
  arrange(desc(abs(estimate))) %>% head(15) %>%   ## keep the best-plotted significant genera
  mutate(panel = "A  DIABIMMUNE (antibiotic)")
ih <- rd_f("real_ihmp_scm.csv") %>%
  mutate(panel = "B  iHMP-IBD (MGX)")
rd <- bind_rows(dh, ih) %>%
  mutate(panel = factor(panel, levels = c("A  DIABIMMUNE (antibiotic)", "B  iHMP-IBD (MGX)")),
         genus = reorder(genus, estimate))

p4 <- ggplot(rd, aes(estimate, genus, xmin = lo, xmax = hi, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2.2) +
  geom_errorbarh(width = 0.25, linewidth = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "#888888", `TRUE` = "#D55E00"),
                      labels = c(`FALSE` = "n.s.", `TRUE` = "FDR < 0.05"),
                      name = NULL) +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(x = "SCM effect estimate (CLR units, post minus counterfactual)", y = NULL) +
  theme_prism(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_y_discrete(guide = "prism_offset")
save2(p4, "fig_real", 7.6, 4.4)

## ---- Fig 5a: normalization choice on the null design (the pivotal sensitivity) ----
nm <- read_csv(file.path(res, "summary_sensitivity_normalize.csv"), show_col_types = FALSE)
n11 <- nm %>% filter(dataset == "A11", !is.na(typeI_fdr)) %>%
  mutate(norm = sub("^norm_", "", config),
         geome = norm %in% c("clr", "ilr", "rclr"),
         norm = reorder(norm, typeI_fdr))
p5 <- ggplot(n11, aes(typeI_fdr, norm, fill = geome)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#BBBBBB"),
                    labels = c(`TRUE` = "Geometric (CLR/ILR/rCLR)", `FALSE` = "Other"),
                    name = NULL) +
  labs(x = "Type I error on A11 (null design)", y = NULL) +
  theme_prism(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(guide = "prism_offset", limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1))
save2(p5, "fig_norm", 6.2, 3.4)

## ---- Fig 5b: ablation variants across designs ----
ab <- read_csv(file.path(res, "summary_ablation.csv"), show_col_types = FALSE) %>%
  filter(!is.na(power_fdr)) %>%
  mutate(design = factor(dataset, levels = c("A1", "A2", "A6", "A11", "A12")),
         variant = factor(variant, levels = c(
           "base", "ab_nnls", "ab_no_offset", "ab_smooth_spline", "ab_smooth_ewma",
           "ab_kalman_gauss", "ab_no_filter", "ab_rclr_nofill", "ab_no_normalize",
           "ab_nb_hier", "ab_nb_nohier")))
p6 <- ggplot(ab, aes(variant, power_fdr, colour = design, group = design)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(A1 = "#000000", A2 = "#E69F00", A6 = "#56B4E9",
                                 A11 = "#009E73", A12 = "#D55E00")) +
  labs(x = NULL, y = "Power (BH-FDR)") +
  base_theme + offset_xd +
  offset_y(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8,
                                    margin = margin(t = 2, b = 2, l = 6, r = 0)))
save2(p6, "fig_ablation", 8, 4.2)

cat("figures written to", figd, "\n")
