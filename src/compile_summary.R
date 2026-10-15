## Compile atlas summary from individual Cox sanity check results
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
out_path <- args[1]

fam_dt <- fread("config/disease_families.tsv")
families <- fam_dt$family_id
dims <- 7
preproc <- "all_comorbid"

rows <- list()
for (fam in families) {
  for (d in dims) {
    path <- sprintf("output/%s/%s/cox/%d/energy_sanity_check.tsv", fam, preproc, d)
    if (file.exists(path)) {
      dt <- fread(path)
      rows[[length(rows) + 1]] <- dt
    }
  }
}

if (length(rows) > 0) {
  result <- rbindlist(rows, fill = TRUE)
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  fwrite(result, out_path, sep = "\t")
  cat(sprintf("Summary written: %s (%d rows)\n", out_path, nrow(result)))
  print(result[, .(disease_family, dim, outcome, HR, pvalue, direction_ok,
                    n_patients, n_basins)])
} else {
  cat("No results found\n")
}
