library(data.table)
library(ggplot2)

outdir <- "plot/Figures"

plot_pair <- function(disease, pcx, pcy, color_by, outfile, n_sample = 200000) {
  f <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", disease)
  cols <- c(pcx, pcy, color_by)
  if (color_by == "sex") cols <- c(pcx, pcy, "sex")
  dt <- fread(f, select = cols)

  if (nrow(dt) > n_sample) dt <- dt[sample(.N, n_sample)]

  if (color_by == "age") {
    p <- ggplot(dt, aes(x = .data[[pcx]], y = .data[[pcy]], colour = age)) +
      geom_point(size = 0.3, alpha = 0.3) +
      scale_colour_viridis_c(option = "plasma", name = "Age") +
      labs(title = sprintf("%s (n=%s)", gsub("_", " ", disease),
                           format(fread(f, select = 1L)[, .N], big.mark = ",")),
           x = pcx, y = pcy) +
      theme_bw(base_size = 16) +
      theme(legend.position = "right",
            plot.title = element_text(size = 18, face = "bold"))
  } else {
    dt[, sex := factor(sex, levels = c("female", "male"))]
    p <- ggplot(dt, aes(x = .data[[pcx]], y = .data[[pcy]], colour = sex)) +
      geom_point(size = 0.3, alpha = 0.2) +
      scale_colour_manual(values = c("female" = "#E64B8B", "male" = "#4B8BE6"), name = "Sex") +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = sprintf("%s (n=%s)", gsub("_", " ", disease),
                           format(fread(f, select = 1L)[, .N], big.mark = ",")),
           x = pcx, y = pcy) +
      theme_bw(base_size = 16) +
      theme(legend.position = "right",
            plot.title = element_text(size = 18, face = "bold"))
  }

  ggsave(outfile, p, width = 7, height = 6, dpi = 300)
  cat(sprintf("Saved: %s\n", outfile))
}

# Age example: diabetes PC1 vs PC2
plot_pair("diabetes", "PC1", "PC2", "age",
          file.path(outdir, "Fig02a_diabetes_PC1_PC2_age.png"))

# Sex examples: try several candidates
# heart_failure PC3 (r=-0.226), stroke PC7 (r=-0.206), UTI PC7 (r=-0.448)
plot_pair("heart_failure", "PC1", "PC3", "sex",
          file.path(outdir, "Fig02_hf_PC1_PC3_sex.png"))
plot_pair("stroke", "PC1", "PC7", "sex",
          file.path(outdir, "Fig02_stroke_PC1_PC7_sex.png"))
plot_pair("urinary_tract_infection", "PC1", "PC7", "sex",
          file.path(outdir, "Fig02_uti_PC1_PC7_sex.png"))

cat("Done.\n")
