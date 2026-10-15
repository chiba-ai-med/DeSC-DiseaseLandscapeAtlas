## Variable-level pairwise analysis
##
## Scatter/box plots between key variables (age, energy, death, hosp, etc.)
## + correlation matrix heatmap. One set per disease family.
##
## Usage:
##   Rscript src/plot_variable_pairs.R <patient_scores_tsv> <disease> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
scores_path <- args[1]
disease     <- args[2]
out_dir     <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== variable pairs: %s ===\n", disease))

pt <- fread(scores_path)
cat(sprintf("  %d patients\n", nrow(pt)))

set.seed(42)
n_sample <- min(50000L, nrow(pt))
pt_sub <- pt[sample.int(nrow(pt), n_sample)]
rm(pt)

# Numeric columns for correlation
num_cols <- c("age", "energy", "n_admissions", "obs_years",
              paste0("PC", 1:7))
num_cols <- intersect(num_cols, names(pt_sub))

# Binary/factor columns
pt_sub[, death_bin := as.integer(death == "dead")]
pt_sub[, hosp_bin := as.integer(hospitalized == "yes")]
pt_sub[, basin_num := as.numeric(basin_id)]

all_num <- c(num_cols, "death_bin", "hosp_bin", "basin_num")
all_num <- intersect(all_num, names(pt_sub))

# =====================================================
# 1. Correlation matrix heatmap
# =====================================================
cat("  Computing correlation matrix ...\n")
mat <- as.matrix(pt_sub[, ..all_num])
mat <- mat[complete.cases(mat), ]
cor_mat <- cor(mat, use = "pairwise.complete.obs")

cor_long <- as.data.table(expand.grid(
  var1 = colnames(cor_mat), var2 = colnames(cor_mat),
  stringsAsFactors = FALSE))
cor_long[, r := as.vector(cor_mat)]

# Pretty labels
label_map <- c(age = "Age", energy = "Energy",
               n_admissions = "N admissions", obs_years = "Obs years",
               death_bin = "Death", hosp_bin = "Hospitalized",
               basin_num = "Basin",
               PC1 = "PC1", PC2 = "PC2", PC3 = "PC3", PC4 = "PC4",
               PC5 = "PC5", PC6 = "PC6", PC7 = "PC7")
cor_long[, v1 := ifelse(var1 %in% names(label_map), label_map[var1], var1)]
cor_long[, v2 := ifelse(var2 %in% names(label_map), label_map[var2], var2)]
ord <- c("Age", "Energy", "Death", "Hospitalized", "N admissions",
         "Obs years", "Basin", paste0("PC", 1:7))
ord <- intersect(ord, unique(cor_long$v1))
cor_long[, v1 := factor(v1, levels = ord)]
cor_long[, v2 := factor(v2, levels = rev(ord))]

g_cor <- ggplot(cor_long, aes(x = v1, y = v2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.5) +
  scale_fill_gradient2(low = "#3a6dd9", mid = "white", high = "#d94a3a",
                        limits = c(-1, 1)) +
  labs(title = sprintf("%s: Variable correlation matrix (n=%d)",
                       disease, nrow(mat)),
       fill = "r") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
ggsave(file.path(out_dir, sprintf("corr_%s.png", disease)), g_cor,
       width = 10, height = 8, dpi = 150)
cat(sprintf("    -> corr_%s.png\n", disease))

# =====================================================
# 2. Key scatter pairs (non-PC variables)
# =====================================================
key_pairs <- list(
  c("age", "energy"),
  c("age", "n_admissions"),
  c("energy", "n_admissions")
)

plots <- list()
for (pr in key_pairs) {
  d <- pt_sub[is.finite(get(pr[1])) & is.finite(get(pr[2]))]
  p <- ggplot(d, aes(x = .data[[pr[1]]], y = .data[[pr[2]]])) +
    geom_point(size = 0.1, alpha = 0.1, color = "gray30") +
    geom_smooth(method = "loess", color = "#d94a3a", se = FALSE) +
    labs(x = pr[1], y = pr[2]) +
    theme_minimal(base_size = 10)
  plots[[length(plots) + 1]] <- p
}

# Scatter colored by death
d <- pt_sub[is.finite(age) & is.finite(energy)]
p_ae_death <- ggplot(d, aes(x = age, y = energy, color = death)) +
  geom_point(size = 0.1, alpha = 0.15) +
  scale_color_manual(values = c(alive = "#3a6dd9", dead = "#d94a3a")) +
  labs(x = "Age", y = "Energy", color = "Death") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
plots[[length(plots) + 1]] <- p_ae_death

# Age-energy colored by hospitalization
p_ae_hosp <- ggplot(d, aes(x = age, y = energy, color = hospitalized)) +
  geom_point(size = 0.1, alpha = 0.15) +
  scale_color_manual(values = c(no = "#3a6dd9", yes = "#d94a3a")) +
  labs(x = "Age", y = "Energy", color = "Hospitalized") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
plots[[length(plots) + 1]] <- p_ae_hosp

# Age-energy colored by basin
p_ae_basin <- ggplot(d, aes(x = age, y = energy, color = basin)) +
  geom_point(size = 0.1, alpha = 0.15) +
  labs(x = "Age", y = "Energy", color = "Basin") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
plots[[length(plots) + 1]] <- p_ae_basin

png(file.path(out_dir, sprintf("var_pairs_%s.png", disease)),
    width = 1800, height = 1200, res = 150)
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
for (pi in seq_along(plots)) {
  row <- ((pi - 1) %/% 3) + 1
  col_pos <- ((pi - 1) %% 3) + 1
  print(plots[[pi]],
        vp = grid::viewport(layout.pos.row = row, layout.pos.col = col_pos))
}
grid::grid.text(
  sprintf("%s: Variable pairwise relationships (n=%d)", disease, n_sample),
  y = 0.99, gp = grid::gpar(fontsize = 13))
dev.off()
cat(sprintf("    -> var_pairs_%s.png\n", disease))

# =====================================================
# 3. Box plots: energy by age group x outcome
# =====================================================
pt_sub[, age_group := cut(age, breaks = c(0, 50, 60, 70, 80, 120),
                           labels = c("<50", "50-59", "60-69", "70-79", "80+"),
                           include.lowest = TRUE)]

# Energy by age group, colored by death
g_box1 <- ggplot(pt_sub[!is.na(age_group) & is.finite(energy)],
                 aes(x = age_group, y = energy, fill = death)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.1) +
  scale_fill_manual(values = c(alive = "#3a6dd9", dead = "#d94a3a")) +
  labs(title = sprintf("%s: Energy by age group and death", disease),
       x = "Age group", y = "Energy", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("box_energy_age_death_%s.png", disease)),
       g_box1, width = 9, height = 5, dpi = 150)
cat(sprintf("    -> box_energy_age_death_%s.png\n", disease))

# Energy by age group, colored by hospitalization
g_box2 <- ggplot(pt_sub[!is.na(age_group) & is.finite(energy)],
                 aes(x = age_group, y = energy, fill = hospitalized)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.1) +
  scale_fill_manual(values = c(no = "#3a6dd9", yes = "#d94a3a")) +
  labs(title = sprintf("%s: Energy by age group and hospitalization", disease),
       x = "Age group", y = "Energy", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("box_energy_age_hosp_%s.png", disease)),
       g_box2, width = 9, height = 5, dpi = 150)
cat(sprintf("    -> box_energy_age_hosp_%s.png\n", disease))

# Basin distribution by age group
g_basin_age <- ggplot(pt_sub[!is.na(age_group) & !is.na(basin)],
                       aes(x = age_group, fill = basin)) +
  geom_bar(position = "fill") +
  labs(title = sprintf("%s: Basin composition by age group", disease),
       x = "Age group", y = "Proportion", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, sprintf("bar_basin_age_%s.png", disease)),
       g_basin_age, width = 8, height = 5, dpi = 150)
cat(sprintf("    -> bar_basin_age_%s.png\n", disease))

cat("Done\n")
