library(data.table)
library(ggplot2)
library(grid)

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
fig_dir <- file.path(base, "plot/Figures/Fig05")

# --- Panel a: Stable vs moved ---
a <- fread(file.path(fig_dir, "Fig05a_trajectory_overview.tsv"))
a_long <- melt(a, id.vars = c("n_basins", "n_diseases"),
               variable.name = "category", value.name = "pct")
a_long$category <- factor(a_long$category,
  levels = c("Moved", "Stable_NonBasin", "Stable_Basin"),
  labels = c("Moved", "Stable in\nNon-Basin", "Stable in\nBasin"))
a_long$n_basins_label <- paste(a_long$n_basins, "basins")
a_long$n_basins_label <- factor(a_long$n_basins_label,
  levels = c("1 basins", "2 basins", "3 basins", "4 basins"))
a_long$label <- paste0(a_long$pct, "%")
a_long[pct < 0.5, label := ""]

pa <- ggplot(a_long, aes(x = n_basins_label, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            colour = "black", fontface = "bold", size = 12) +
  scale_fill_manual(values = c("Moved" = "#5CB85C", "Stable in\nNon-Basin" = "#D9534F",
                                "Stable in\nBasin" = "#337AB7"),
                    name = "Trajectory pattern") +
  labs(title = "Stable vs. moved before death", x = NULL, y = "Patients (%)") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 40),
        axis.text = element_text(size = 32),
        axis.title = element_text(size = 36),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "Fig05a_trajectory_overview.png"), pa,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved Fig05a\n")

# --- Panel b: Starting in basin ---
b <- fread(file.path(fig_dir, "Fig05b_trajectory_basin_start.tsv"))
b_long <- melt(b, id.vars = c("n_basins", "n_diseases"),
               variable.name = "category", value.name = "pct")
b_long$category <- factor(b_long$category,
  levels = c("To_NonBasin", "Diff_Basin", "Same_Basin"),
  labels = c("To Non-Basin", "Different Basin", "Same Basin"))
b_long$n_basins_label <- factor(paste(b_long$n_basins, "basins"),
  levels = c("1 basins", "2 basins", "3 basins", "4 basins"))
b_long$label <- paste0(b_long$pct, "%")
b_long[pct < 0.5, label := ""]

pb <- ggplot(b_long, aes(x = n_basins_label, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            colour = "black", fontface = "bold", size = 12) +
  scale_fill_manual(values = c("To Non-Basin" = "#D9534F", "Different Basin" = "#F0AD4E",
                                "Same Basin" = "#337AB7"),
                    name = "End position") +
  labs(title = "Starting in basin", x = NULL, y = "Patients (%)") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 40),
        axis.text = element_text(size = 32),
        axis.title = element_text(size = 36),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "Fig05b_trajectory_basin_start.png"), pb,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved Fig05b\n")

# --- Panel c: Starting outside basin ---
cc <- fread(file.path(fig_dir, "Fig05c_trajectory_nonbasin_start.tsv"))
cc_long <- melt(cc, id.vars = c("n_basins", "n_diseases"),
                variable.name = "category", value.name = "pct")
cc_long$category <- factor(cc_long$category,
  levels = c("Stay_NonBasin", "To_Basin"),
  labels = c("Stay Non-Basin", "To Basin"))
cc_long$n_basins_label <- factor(paste(cc_long$n_basins, "basins"),
  levels = c("1 basins", "2 basins", "3 basins", "4 basins"))
cc_long$label <- paste0(cc_long$pct, "%")

pc <- ggplot(cc_long, aes(x = n_basins_label, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            colour = "black", fontface = "bold", size = 12) +
  scale_fill_manual(values = c("Stay Non-Basin" = "#D9534F", "To Basin" = "#337AB7"),
                    name = "End position") +
  labs(title = "Starting outside basin", x = NULL, y = "Patients (%)") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 40),
        axis.text = element_text(size = 32),
        axis.title = element_text(size = 36),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "Fig05c_trajectory_nonbasin_start.png"), pc,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved Fig05c\n")

# --- Composite (legend fix: smaller text + keys to prevent truncation) ---
legend_theme <- theme(legend.position = "bottom",
                      legend.text = element_text(size = 132),
                      legend.title = element_text(size = 144, face = "bold"),
                      legend.key.size = unit(5, "cm"),
                      legend.spacing.x = unit(1, "cm"),
                      legend.margin = margin(0, 0, 0, 0),
                      plot.margin = margin(10, 80, 10, 10))

ga <- ggplotGrob(pa + legend_theme)
gb <- ggplotGrob(pb + legend_theme)
gc <- ggplotGrob(pc + legend_theme)

png(file.path(fig_dir, "Fig05_trajectory.png"), width = 21600, height = 12000, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3, widths = unit(c(1, 1.3, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, clip = "off"))
grid.text("a", x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
grid.draw(ga)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, clip = "off"))
grid.text("b", x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
grid.draw(gb)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, clip = "off"))
grid.text("c", x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
grid.draw(gc)
popViewport()

dev.off()
cat("Saved composite Fig05\n")
