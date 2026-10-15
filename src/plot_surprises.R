## Surprise analysis: what does the landscape capture beyond age?
##
## 1. Age-residualized energy vs outcomes
## 2. PC loadings (top ICD-10 codes per PC)
## 3. Basin comorbidity profiles (within age strata)
## 4. Discordant case characterization
##
## Usage:
##   Rscript src/plot_surprises.R \
##     <patient_scores_tsv> <eigenvectors_csv> <col_ids_txt> \
##     <disease> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
scores_path <- args[1]
eigvec_path <- args[2]
col_ids_path <- args[3]
disease     <- args[4]
out_dir     <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== surprise analysis: %s ===\n", disease))

pt <- fread(scores_path)
cat(sprintf("  %d patients\n", nrow(pt)))

# =====================================================
# 1. Age-residualized energy
# =====================================================
cat("  1. Age-residualized energy ...\n")

pt_valid <- pt[is.finite(age) & is.finite(energy) & is.finite(as.numeric(obs_years))]
pt_valid[, time_years := as.numeric(obs_years)]
pt_valid[, death_bin := as.integer(death == "dead")]
pt_valid[, hosp_bin := as.integer(hospitalized == "yes")]
pt_valid[, sex_f := factor(sex)]

# Regress out age from energy
fit_resid <- lm(energy ~ age, data = pt_valid)
pt_valid[, E_resid := residuals(fit_resid)]
pt_valid[, E_resid_z := scale(E_resid)[, 1]]

# Cox with raw energy
fit_raw <- coxph(Surv(time_years, death_bin) ~ scale(energy)[,1] + age + sex_f,
                 data = pt_valid[time_years > 0])
hr_raw <- summary(fit_raw)$conf.int[1, ]
pv_raw <- coef(summary(fit_raw))[1, 5]

# Cox with age-residualized energy
fit_res <- coxph(Surv(time_years, death_bin) ~ E_resid_z + age + sex_f,
                 data = pt_valid[time_years > 0])
hr_res <- summary(fit_res)$conf.int[1, ]
pv_res <- coef(summary(fit_res))[1, 5]

cat(sprintf("  Raw energy:     HR=%.4f [%.4f, %.4f] p=%.2e\n",
            hr_raw[1], hr_raw[3], hr_raw[4], pv_raw))
cat(sprintf("  Age-residual:   HR=%.4f [%.4f, %.4f] p=%.2e\n",
            hr_res[1], hr_res[3], hr_res[4], pv_res))

# Partial R² of age in explaining energy
r2_age <- summary(fit_resid)$r.squared
cat(sprintf("  R²(energy ~ age) = %.4f (age explains %.1f%% of energy variance)\n",
            r2_age, 100 * r2_age))

# Save summary
resid_summary <- data.table(
  disease = disease,
  metric = c("HR_raw_energy", "HR_age_residual_energy", "R2_energy_age"),
  value = c(hr_raw[1], hr_res[1], r2_age),
  CI_lower = c(hr_raw[3], hr_res[3], NA),
  CI_upper = c(hr_raw[4], hr_res[4], NA),
  pvalue = c(pv_raw, pv_res, NA)
)
fwrite(resid_summary, file.path(out_dir, sprintf("resid_summary_%s.tsv", disease)),
       sep = "\t")

# Plot: age-residualized energy by death
set.seed(42)
n_s <- min(30000L, nrow(pt_valid))
ps <- pt_valid[sample.int(nrow(pt_valid), n_s)]

g1 <- ggplot(ps, aes(x = age, y = E_resid, color = death)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  scale_color_manual(values = c(alive = "#3a6dd9", dead = "#d94a3a")) +
  labs(title = sprintf("%s: Age-residualized energy by death\nR²(energy~age)=%.3f, HR(resid)=%.3f",
                       disease, r2_age, hr_res[1]),
       x = "Age", y = "Energy residual (age regressed out)", color = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("resid_energy_death_%s.png", disease)),
       g1, width = 9, height = 6, dpi = 150)

# By hospitalization
g2 <- ggplot(ps, aes(x = age, y = E_resid, color = hospitalized)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  scale_color_manual(values = c(no = "#3a6dd9", yes = "#d94a3a")) +
  labs(title = sprintf("%s: Age-residualized energy by hospitalization",
                       disease),
       x = "Age", y = "Energy residual", color = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("resid_energy_hosp_%s.png", disease)),
       g2, width = 9, height = 6, dpi = 150)

# =====================================================
# 2. PC loadings: top ICD-10 codes per PC
# =====================================================
cat("  2. PC loadings ...\n")

eigvec <- as.matrix(read.csv(eigvec_path, header = FALSE))

col_lines <- readLines(col_ids_path, encoding = "UTF-8")
code_names <- character(nrow(eigvec))
for (ln in col_lines) {
  parts <- strsplit(ln, ",")[[1]]
  if (length(parts) >= 2) {
    code <- trimws(parts[1])
    id <- as.integer(trimws(parts[2]))
    if (!is.na(id) && id >= 1 && id <= length(code_names)) {
      code_names[id] <- code
    }
  }
}

n_top <- 20
loading_plots <- list()
for (pc in 1:7) {
  v <- eigvec[, pc]
  top_pos <- order(-v)[1:n_top]
  top_neg <- order(v)[1:n_top]
  top_all <- unique(c(top_pos, top_neg))
  ld <- data.table(
    code = code_names[top_all],
    loading = v[top_all],
    pc = paste0("PC", pc)
  )
  ld <- ld[order(-abs(loading))][1:min(n_top, nrow(ld))]
  ld[, code := factor(code, levels = rev(code))]

  g <- ggplot(ld, aes(x = loading, y = code,
                       fill = ifelse(loading > 0, "pos", "neg"))) +
    geom_col() +
    scale_fill_manual(values = c(pos = "#d94a3a", neg = "#3a6dd9"),
                       guide = "none") +
    labs(title = sprintf("PC%d top loadings", pc), x = "Loading", y = NULL) +
    theme_minimal(base_size = 9)
  loading_plots[[pc]] <- g
}

# Save as 7 panels
png(file.path(out_dir, sprintf("pc_loadings_%s.png", disease)),
    width = 2400, height = 2400, res = 150)
grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 3)))
for (pc in 1:7) {
  row <- ((pc - 1) %/% 3) + 1
  col_pos <- ((pc - 1) %% 3) + 1
  print(loading_plots[[pc]],
        vp = grid::viewport(layout.pos.row = row, layout.pos.col = col_pos))
}
grid::grid.text(
  sprintf("%s: Top ICD-10 loadings per PC (dim=7)", disease),
  y = 0.995, gp = grid::gpar(fontsize = 14))
dev.off()
cat(sprintf("    -> pc_loadings_%s.png\n", disease))

# =====================================================
# 3. Discordant cases: young + low energy, old + high energy
# =====================================================
cat("  3. Discordant cases ...\n")

pt_valid[, age_group := cut(age, breaks = c(0, 50, 60, 70, 80, 120),
                             labels = c("<50", "50-59", "60-69", "70-79", "80+"),
                             include.lowest = TRUE)]
pt_valid[, E_tertile := cut(energy,
  breaks = quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("low", "mid", "high"), include.lowest = TRUE)]

# Cross-tab: age group x energy tertile x death rate
crosstab <- pt_valid[!is.na(age_group) & !is.na(E_tertile),
  .(n = .N,
    death_rate = mean(death_bin, na.rm = TRUE),
    hosp_rate = mean(hosp_bin, na.rm = TRUE)),
  by = .(age_group, E_tertile)]
fwrite(crosstab, file.path(out_dir, sprintf("crosstab_%s.tsv", disease)),
       sep = "\t")
cat("  Cross-tab (age x energy x outcomes):\n")
print(crosstab[order(age_group, E_tertile)])

# Heatmap: death rate by age group x energy tertile
g3 <- ggplot(crosstab, aes(x = E_tertile, y = age_group, fill = death_rate)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%\nn=%s",
                                 100 * death_rate,
                                 formatC(n, big.mark = ","))),
            size = 3.5) +
  scale_fill_gradient(low = "white", high = "#d94a3a") +
  labs(title = sprintf("%s: Death rate by age group x energy tertile", disease),
       x = "Energy tertile", y = "Age group", fill = "Death rate") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("heatmap_age_energy_death_%s.png", disease)),
       g3, width = 7, height = 6, dpi = 150)
cat(sprintf("    -> heatmap_age_energy_death_%s.png\n", disease))

# Same for hospitalization
g4 <- ggplot(crosstab, aes(x = E_tertile, y = age_group, fill = hosp_rate)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%\nn=%s",
                                 100 * hosp_rate,
                                 formatC(n, big.mark = ","))),
            size = 3.5) +
  scale_fill_gradient(low = "white", high = "#d94a3a") +
  labs(title = sprintf("%s: Hospitalization rate by age group x energy tertile",
                       disease),
       x = "Energy tertile", y = "Age group", fill = "Hosp rate") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("heatmap_age_energy_hosp_%s.png", disease)),
       g4, width = 7, height = 6, dpi = 150)

# =====================================================
# 4. Within-age-stratum basin comparison
# =====================================================
cat("  4. Basin comparison within age strata ...\n")

basin_by_age <- pt_valid[!is.na(age_group) & !is.na(basin),
  .(n = .N,
    death_rate = mean(death_bin, na.rm = TRUE),
    hosp_rate = mean(hosp_bin, na.rm = TRUE),
    mean_energy = mean(energy, na.rm = TRUE),
    mean_admissions = mean(n_admissions, na.rm = TRUE)),
  by = .(age_group, basin)]
fwrite(basin_by_age,
       file.path(out_dir, sprintf("basin_by_age_%s.tsv", disease)),
       sep = "\t")
cat("  Basin outcomes within age strata:\n")
print(basin_by_age[order(age_group, basin)])

g5 <- ggplot(basin_by_age, aes(x = age_group, y = death_rate,
                                fill = basin, group = basin)) +
  geom_col(position = "dodge") +
  labs(title = sprintf("%s: Death rate by basin within age group", disease),
       x = "Age group", y = "Death rate", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("basin_death_by_age_%s.png", disease)),
       g5, width = 9, height = 5, dpi = 150)

cat("Done\n")
