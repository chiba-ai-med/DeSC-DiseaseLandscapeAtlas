## Post-process: verify cohort table consistency with Cox results
## Reads landmark_cohort.tsv per disease and checks counts against summary TSVs
##
## Usage: Rscript src/post_process_landmark.R

suppressPackageStartupMessages(library(data.table))

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

cat("=== Verifying landmark cohort tables vs results ===\n\n")

overall <- fread("summary/landmark_all_diseases_overall.tsv")
per_basin <- fread("summary/landmark_all_diseases_per_basin.tsv")
atlas <- fread("summary/atlas_summary.tsv")

diseases <- overall[direction %in% c("PASS", "FAIL"), disease]
n_ok <- 0L
n_fail <- 0L

for (d in diseases) {
  cohort_path <- file.path("output", d, "all_comorbid", "landmark", "landmark_cohort.tsv")
  if (!file.exists(cohort_path)) {
    cat(sprintf("  MISSING: %s\n", cohort_path))
    n_fail <- n_fail + 1L
    next
  }

  dt <- fread(cohort_path)
  ov <- overall[disease == d]

  # Basin check: all > 0
  if (any(dt$basin_id <= 0, na.rm = TRUE)) {
    cat(sprintf("  ERROR %s: basin_id <= 0 found in cohort table\n", d))
    n_fail <- n_fail + 1L
  }

  # Mortality population (riyouka=1)
  mort <- dt[shibou_flg_riyouka == 1L]
  mort_n <- nrow(mort)
  mort_d <- sum(mort$event_death, na.rm = TRUE)
  mort_c <- mort_n - mort_d

  if (mort_n != ov$N || mort_d != ov$deaths || mort_c != ov$censored) {
    cat(sprintf("  MISMATCH %s: cohort N=%d/d=%d/c=%d vs result N=%d/d=%d/c=%d\n",
                d, mort_n, mort_d, mort_c, ov$N, ov$deaths, ov$censored))
    n_fail <- n_fail + 1L
  } else {
    n_ok <- n_ok + 1L
  }

  # Per-basin N/deaths check
  pb <- per_basin[disease == d]
  for (b in unique(pb$basin_id)) {
    sub_c <- mort[basin_id == b]
    sub_r <- pb[basin_id == b]
    if (nrow(sub_c) != sub_r$N || sum(sub_c$event_death) != sub_r$deaths) {
      cat(sprintf("  BASIN MISMATCH %s basin %d: cohort N=%d/d=%d vs result N=%d/d=%d\n",
                  d, b, nrow(sub_c), sum(sub_c$event_death), sub_r$N, sub_r$deaths))
      n_fail <- n_fail + 1L
    }
  }
}

# ==========================================================
# Hospitalization verification
# ==========================================================
hosp_overall <- fread("summary/landmark_all_diseases_hosp_overall.tsv")
hosp_per_basin <- fread("summary/landmark_all_diseases_hosp_per_basin.tsv")
hosp_diseases <- hosp_overall[direction %in% c("PASS", "FAIL"), disease]

cat("\n=== Hospitalization verification ===\n")
n_hosp_ok <- 0L
n_hosp_fail <- 0L
excl_counts <- list()

for (d in hosp_diseases) {
  cohort_path <- file.path("output", d, "all_comorbid", "landmark", "landmark_cohort.tsv")
  if (!file.exists(cohort_path)) {
    cat(sprintf("  MISSING: %s\n", cohort_path))
    n_hosp_fail <- n_hosp_fail + 1L
    next
  }

  dt <- fread(cohort_path)
  ov_h <- hosp_overall[disease == d]

  # Pre-landmark exclusion count
  n_excluded <- sum(dt$hosp_excluded, na.rm = TRUE)
  excl_counts[[d]] <- n_excluded

  # Hosp population = all eligible minus pre-landmark admissions
  hosp_pop <- dt[hosp_excluded == FALSE]
  hosp_n <- nrow(hosp_pop)
  hosp_ev <- sum(hosp_pop$event_hosp, na.rm = TRUE)
  hosp_cens <- hosp_n - hosp_ev

  if (hosp_n != ov_h$N || hosp_ev != ov_h$events || hosp_cens != ov_h$censored) {
    cat(sprintf("  HOSP MISMATCH %s: cohort N=%d/ev=%d/c=%d vs result N=%d/ev=%d/c=%d\n",
                d, hosp_n, hosp_ev, hosp_cens, ov_h$N, ov_h$events, ov_h$censored))
    n_hosp_fail <- n_hosp_fail + 1L
  } else {
    n_hosp_ok <- n_hosp_ok + 1L
  }

  # Per-basin hosp check
  pb_h <- hosp_per_basin[disease == d]
  for (b in unique(pb_h$basin_id)) {
    sub_c <- hosp_pop[basin_id == b]
    sub_r <- pb_h[basin_id == b]
    if (nrow(sub_c) != sub_r$N || sum(sub_c$event_hosp) != sub_r$events) {
      cat(sprintf("  HOSP BASIN MISMATCH %s basin %d: cohort N=%d/ev=%d vs result N=%d/ev=%d\n",
                  d, b, nrow(sub_c), sum(sub_c$event_hosp), sub_r$N, sub_r$events))
      n_hosp_fail <- n_hosp_fail + 1L
    }
  }
}

cat(sprintf("\nHospitalization verification: %d OK, %d failures\n", n_hosp_ok, n_hosp_fail))
n_fail <- n_fail + n_hosp_fail

# Exclusion summary
cat("\n=== Pre-landmark hospitalization exclusion counts ===\n")
for (d in names(excl_counts)) {
  if (excl_counts[[d]] > 0) {
    cat(sprintf("  %s: %d patients excluded\n", d, excl_counts[[d]]))
  }
}
cat(sprintf("  Total diseases with exclusions: %d/%d\n",
            sum(sapply(excl_counts, function(x) x > 0)), length(excl_counts)))

# Atlas n_basins check
for (d in atlas$disease_family) {
  a_nb <- atlas[disease_family == d, n_basins]
  pb_nb <- per_basin[disease == d, uniqueN(basin_id)]
  if (a_nb != pb_nb) {
    cat(sprintf("  ATLAS n_basins mismatch %s: atlas=%d vs per_basin distinct=%d\n",
                d, a_nb, pb_nb))
    n_fail <- n_fail + 1L
  }
}

cat(sprintf("\n=== Verification: %d OK, %d failures ===\n", n_ok, n_fail))
if (n_fail > 0) quit(status = 1)
