## Stacked area plot of ICD-10 code composition per basin before event
##
## Reads basin_code_trajectory_{endpoint}.tsv and produces stacked area plots.
##
## Usage:
##   Rscript src/plot_code_trajectory_stacked.R <disease> <out_dir> <endpoint>
##   endpoint: "death" or "hosp"

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease  <- args[1]
out_dir  <- args[2]
endpoint <- args[3]

cat(sprintf("=== Stacked area: %s (%s) ===\n", disease, endpoint))

tsv_path <- file.path(out_dir,
  sprintf("basin_code_trajectory_%s.tsv", endpoint))
dt <- fread(tsv_path)

event_label <- if (endpoint == "death") "death" else "hospitalization"
basins <- sort(unique(dt$basin_id[dt$basin_id > 0L]))

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

for (b in basins) {
  df <- dt[basin_id == b]
  if (nrow(df) < 10) next

  top_codes <- df[, .(total_freq = sum(freq)), by = .(code, en_label)]
  setorder(top_codes, -total_freq)
  top15 <- top_codes[1:min(15, nrow(top_codes))]

  df_top <- df[code %in% top15$code]
  df_top[, en_label := factor(en_label, levels = rev(top15$en_label))]

  # Stacked area
  g <- ggplot(df_top, aes(x = -periods, y = freq, fill = en_label)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    scale_x_continuous(
      breaks = seq(-19, 0, by = 2),
      labels = function(x) ifelse(x == 0, event_label, as.character(-x))) +
    scale_fill_manual(values = rep(
      c("#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
        "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
        "#af7aa1","#86bcb6","#d37295","#fabfd2","#b6992d"), 2)[1:nrow(top15)],
      name = "ICD-10") +
    labs(title = sprintf("%s / Basin %d: disease code composition before %s",
                          disease, b, event_label),
         subtitle = "Top 15 codes by total frequency",
         x = sprintf("Periods before %s", event_label),
         y = "Frequency (per patient)") +
    theme_pub +
    theme(legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))

  png_path <- file.path(out_dir,
    sprintf("basin_code_stacked_%s_basin%d.png", endpoint, b))
  ggsave(png_path, g, width = 14, height = 8, dpi = 200, bg = "white")
  cat(sprintf("  Written: %s\n", png_path))
}

cat("Done.\n")
