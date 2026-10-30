## Recompute all Fig05 data with CORRECT definitions
## Panel b: same_basin/diff_basin/to_nonbasin among ALL basin-start patients
library(data.table)

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
families <- fread(file.path(base, "config/disease_families.tsv"))
diseases <- families$family_id[families$family_id != "hf_cancer"]

classify_disease <- function(d) {
  traj_path <- file.path(base, "plot", d, "trajectories_dead.tsv")
  if (!file.exists(traj_path)) return(NULL)
  traj <- fread(traj_path, select = c("kojin_id", "basin_id", "is_basin", "periods_before_death"))
  if (nrow(traj) == 0) return(NULL)

  basin_path <- file.path(base, "output", d, "all_comorbid/landscaper/7/Basin.tsv")
  nb <- length(scan(basin_path, quiet = TRUE))

  pp <- traj[, .(
    n_steps = .N,
    all_basin = all(is_basin == TRUE),
    all_nonbasin = all(is_basin == FALSE),
    first_basin = is_basin[which.max(periods_before_death)],
    last_basin = is_basin[which.min(periods_before_death)],
    first_basin_id = basin_id[which.max(periods_before_death)],
    last_basin_id = basin_id[which.min(periods_before_death)],
    n_unique_basins = uniqueN(basin_id[is_basin == TRUE])
  ), by = kojin_id]

  n_dead <- nrow(pp)

  # Panel a: all patients
  stable_basin_n <- pp[all_basin == TRUE & n_unique_basins == 1, .N]
  stable_nonbasin_n <- pp[all_nonbasin == TRUE, .N]
  moved_n <- n_dead - stable_basin_n - stable_nonbasin_n

  stable_basin_pct <- stable_basin_n / n_dead * 100
  stable_nonbasin_pct <- stable_nonbasin_n / n_dead * 100
  moved_pct <- moved_n / n_dead * 100

  # Panel b: among ALL basin-start patients
  bs <- pp[first_basin == TRUE]
  n_bs <- nrow(bs)
  if (n_bs > 0) {
    same_basin_pct <- bs[first_basin_id == last_basin_id & last_basin == TRUE, .N] / n_bs * 100
    diff_basin_pct <- bs[first_basin_id != last_basin_id & last_basin == TRUE, .N] / n_bs * 100
    to_nonbasin_pct <- bs[last_basin == FALSE, .N] / n_bs * 100
  } else {
    same_basin_pct <- diff_basin_pct <- to_nonbasin_pct <- NA_real_
  }

  # Panel c: among nonbasin-start
  nbs <- pp[first_basin == FALSE]
  n_nbs <- nrow(nbs)
  if (n_nbs > 0) {
    to_basin_pct <- nbs[last_basin == TRUE, .N] / n_nbs * 100
    stay_nonbasin_pct <- nbs[last_basin == FALSE, .N] / n_nbs * 100
  } else {
    to_basin_pct <- stay_nonbasin_pct <- NA_real_
  }

  # Supplementary trajectory_classification categories
  always_nonbasin <- stable_nonbasin_pct
  basin_to_nonbasin <- pp[first_basin == TRUE & last_basin == FALSE, .N] / n_dead * 100
  basin_transition <- pp[first_basin == TRUE & last_basin == TRUE &
                          first_basin_id != last_basin_id, .N] / n_dead * 100
  nonbasin_to_basin <- pp[first_basin == FALSE & last_basin == TRUE, .N] / n_dead * 100
  stable_in_basin <- stable_basin_pct
  mixed <- 100 - always_nonbasin - basin_to_nonbasin - basin_transition -
           nonbasin_to_basin - stable_in_basin

  data.table(
    disease = d, n_basins = nb, n_dead = n_dead,
    n_basin_start = n_bs, n_nonbasin_start = n_nbs,
    stable_basin = stable_basin_pct, stable_nonbasin = stable_nonbasin_pct, moved = moved_pct,
    same_basin = same_basin_pct, diff_basin = diff_basin_pct, to_nonbasin = to_nonbasin_pct,
    to_basin = to_basin_pct, stay_nonbasin = stay_nonbasin_pct,
    always_nonbasin_s = always_nonbasin, basin_to_nonbasin_s = basin_to_nonbasin,
    basin_transition_s = basin_transition, nonbasin_to_basin_s = nonbasin_to_basin,
    stable_in_basin_s = stable_in_basin, mixed_s = mixed
  )
}

cat("Processing 51 diseases...\n")
dt <- rbindlist(lapply(diseases, function(d) {
  cat(sprintf("  %s\n", d))
  classify_disease(d)
}))

# Write detail TSV
out_dir <- file.path(base, "plot/Figures/Fig05")
fwrite(dt[, .(disease, n_basins, n_dead,
              stable_basin, stable_nonbasin, moved,
              same_basin, diff_basin, to_nonbasin,
              to_basin, stay_nonbasin)],
       file.path(out_dir, "Fig05b_trajectory_basin_start_detail.tsv"), sep = "\t")

# Summary TSVs for the figure panels
fig05a <- dt[, .(n_diseases = .N,
                  Stable_Basin = round(median(stable_basin), 1),
                  Stable_NonBasin = round(median(stable_nonbasin), 1),
                  Moved = round(median(moved), 1)),
             by = n_basins]
setorder(fig05a, n_basins)

fig05b <- dt[, .(n_diseases = .N,
                  Same_Basin = round(median(same_basin, na.rm = TRUE), 1),
                  Diff_Basin = round(median(diff_basin, na.rm = TRUE), 1),
                  To_NonBasin = round(median(to_nonbasin, na.rm = TRUE), 1)),
             by = n_basins]
setorder(fig05b, n_basins)

fig05c <- dt[!is.na(to_basin), .(n_diseases = .N,
                  To_Basin = round(median(to_basin, na.rm = TRUE), 1),
                  Stay_NonBasin = round(median(stay_nonbasin, na.rm = TRUE), 1)),
             by = n_basins]
setorder(fig05c, n_basins)

fwrite(fig05a, file.path(out_dir, "Fig05a_trajectory_overview.tsv"), sep = "\t")
fwrite(fig05b, file.path(out_dir, "Fig05b_trajectory_basin_start.tsv"), sep = "\t")
fwrite(fig05c, file.path(out_dir, "Fig05c_trajectory_nonbasin_start.tsv"), sep = "\t")

# Supplementary trajectory_classification.tsv
suppl <- dt[, .(disease, n_dead,
                always_nonbasin = always_nonbasin_s,
                basin_to_nonbasin = basin_to_nonbasin_s,
                basin_transition = basin_transition_s,
                mixed = mixed_s,
                nonbasin_to_basin = nonbasin_to_basin_s,
                stable_in_basin = stable_in_basin_s,
                n_basins)]
fwrite(suppl, file.path(base, "plot/Figures/Supplementary/trajectory_classification.tsv"), sep = "\t")

cat("\n=== Panel a summary ===\n"); print(fig05a)
cat("\n=== Panel b summary ===\n"); print(fig05b)
cat("\n=== Panel c summary ===\n"); print(fig05c)

# Key statistics for manuscript
cat("\n=== Manuscript statistics ===\n")
s1 <- dt[n_basins == 1]
s2 <- dt[n_basins == 2]
s34 <- dt[n_basins >= 3]

cat(sprintf("1-basin: median same_basin = %.1f%% (n=%d diseases)\n",
            median(s1$same_basin, na.rm=TRUE), nrow(s1)))
cat(sprintf("2-basin: median same_basin = %.1f%% (n=%d diseases)\n",
            median(s2$same_basin, na.rm=TRUE), nrow(s2)))
cat(sprintf("3-4-basin: median same_basin = %.1f%% (n=%d diseases)\n",
            median(s34$same_basin, na.rm=TRUE), nrow(s34)))

cat(sprintf("\n1-basin: median stable_basin (panel a) = %.1f%%\n",
            median(s1$stable_basin)))
cat(sprintf("All: median stable_nonbasin = %.1f%% (max=%.1f%% in %s)\n",
            median(dt$stable_nonbasin), max(dt$stable_nonbasin),
            dt$disease[which.max(dt$stable_nonbasin)]))

cat(sprintf("\nNon-basin-start → basin (panel c):\n"))
for (nb in sort(unique(dt$n_basins))) {
  sub <- dt[n_basins == nb & !is.na(to_basin)]
  if (nrow(sub) > 0) {
    cat(sprintf("  %d-basin: median to_basin = %.1f%% (n=%d diseases)\n",
                nb, median(sub$to_basin), nrow(sub)))
  }
}

# Diff basin: how many move to different basin
cat(sprintf("\n3-4-basin: median diff_basin = %.1f%%\n",
            median(s34$diff_basin, na.rm=TRUE)))

cat("\nDone. All TSVs written.\n")
