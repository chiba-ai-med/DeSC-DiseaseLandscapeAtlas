## BasinNA Bidirectionality Forest Plot
## For each disease, shows the adjusted HR of the largest basin vs BasinNA
## Color: Type A (BasinNA=high risk, HR<1) vs Type B (BasinNA=low risk, HR>1)
##
## Usage: Rscript src/plot_basinna_bidirectionality.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cox <- fread("summary/basin_exploration/cox_basin_membership.tsv")
report <- fread("summary/basin_exploration/per_disease_basin_report.tsv")

# Keep only diseases where BasinNA (ref="0") is the reference
cox_na <- cox[ref == "0"]

# From the report, find the largest non-BasinNA basin per disease
largest <- report[basin_id != 0,
                  .SD[which.max(n_patients)],
                  by = disease][, .(disease, largest_basin = basin_id)]

# Join: get the HR of the largest basin vs BasinNA for each disease
plt <- merge(cox_na, largest, by = "disease")
plt <- plt[basin_id == as.character(largest_basin)]

# Classify: HR < 1 means basin is safer than BasinNA (Type A: BasinNA = high risk)
#           HR > 1 means basin is riskier than BasinNA (Type B: BasinNA = low risk)
plt[, type := ifelse(HR < 1, "Type A: BasinNA = high risk",
                              "Type B: BasinNA = low risk")]
plt[, sig := p < 0.05]

# Sort by HR
plt[, disease_label := gsub("_", " ", disease)]
plt[, disease_label := tools::toTitleCase(disease_label)]
setorder(plt, HR)
plt[, disease_label := factor(disease_label, levels = disease_label)]

# Colors: warm red for Type A (BasinNA dangerous), cool blue for Type B (BasinNA safe)
type_colors <- c("Type A: BasinNA = high risk" = "#C4463A",
                 "Type B: BasinNA = low risk"  = "#3A7CC4")

p <- ggplot(plt, aes(x = HR, y = disease_label, color = type)) +
  geom_vline(xintercept = 1.0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                 height = 0.3, linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.5) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = "") +
  scale_color_manual(values = type_colors, name = "") +
  scale_x_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5),
                labels = c("0.6","0.7","0.8","0.9","1.0","1.1","1.2","1.3","1.5")) +
  labs(
    title = "BasinNA Bidirectionality",
    subtitle = "Adjusted HR of largest basin vs BasinNA (age + sex adjusted Cox)",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.y = element_text(size = 7),
    plot.margin = margin(10, 15, 10, 5)
  ) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))

out_path <- "plot/basinna_bidirectionality_forest.png"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path, p, width = 8, height = 12, dpi = 200, bg = "white")
cat(sprintf("Saved: %s\n", out_path))

# Summary stats
cat(sprintf("\nType A (BasinNA=high risk, basin HR<1): %d diseases\n",
            sum(plt$HR < 1)))
cat(sprintf("Type B (BasinNA=low risk, basin HR>1): %d diseases\n",
            sum(plt$HR > 1)))
cat(sprintf("Significant (p<0.05): %d / %d\n",
            sum(plt$p < 0.05), nrow(plt)))
