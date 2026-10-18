## Compute Charlson (CCI) and Elixhauser comorbidity indices per patient
##
## Reads the binary COO matrix to extract patient-level ICD-10 codes,
## then applies Quan et al. 2005 ICD-10 mappings for CCI and
## van Walraven-weighted Elixhauser index.
##
## Usage:
##   Rscript src/compute_comorbidity_indices.R <disease> <base_dir> <out_tsv>

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
disease  <- args[1]
base_dir <- args[2]
out_tsv  <- args[3]

cat(sprintf("=== Comorbidity indices: %s ===\n", disease))

# --- Read COO matrix (binary zst) ---
coo_path <- file.path(base_dir, "cohort", "feature_matrix.coo.zst")
cat(sprintf("  Reading COO: %s\n", coo_path))

# --- Read col_ids (col_id -> ICD-10 code) ---
col_map <- fread(file.path(base_dir, "cohort", "col_ids.txt"),
                 header = FALSE, col.names = c("code", "col_id"))

# --- Read person_map (new_row_id -> kojin_id) ---
pmap <- fread(file.path(base_dir, "cohort", "person_map.tsv"))
setkey(pmap, new_row_id)
setkey(col_map, col_id)

# --- Read COO in chunks to handle large files ---
con <- pipe(sprintf("zstd -dcq '%s'", coo_path), "rb")
header <- readBin(con, "integer", n = 2, size = 4, endian = "little")
nrows <- header[1]; ncols <- header[2]
cat(sprintf("  header: nrows=%d, ncols=%d\n", nrows, ncols))

chunk_size <- 200000000L  # 100M pairs per chunk (200M ints)
person_codes_list <- list()
chunk_i <- 0L
total_pairs <- 0L
repeat {
  raw <- readBin(con, "integer", n = chunk_size, size = 4, endian = "little")
  if (length(raw) == 0) break
  chunk_i <- chunk_i + 1L
  npairs <- length(raw) %/% 2L
  total_pairs <- total_pairs + npairs
  chunk_dt <- data.table(
    new_row_id = raw[seq(1L, 2L * npairs, 2L)],
    col_id     = raw[seq(2L, 2L * npairs, 2L)]
  )
  rm(raw); gc(verbose = FALSE)
  chunk_dt <- merge(chunk_dt, pmap[, .(new_row_id, kojin_id)],
                    by = "new_row_id", all.x = FALSE)
  chunk_dt <- merge(chunk_dt, col_map, by = "col_id", all.x = FALSE)
  person_codes_list[[chunk_i]] <- unique(chunk_dt[, .(kojin_id, code)])
  rm(chunk_dt); gc(verbose = FALSE)
  cat(sprintf("  chunk %d: %d pairs (total %d)\n", chunk_i, npairs, total_pairs))
}
close(con)
cat(sprintf("  %d pairs read in %d chunks\n", total_pairs, chunk_i))

person_codes <- unique(rbindlist(person_codes_list))
rm(person_codes_list); gc(verbose = FALSE)
cat(sprintf("  %d unique person-code pairs\n", nrow(person_codes)))

# --- Simple comorbidity count ---
n_codes <- person_codes[, .(n_codes = uniqueN(code)), by = kojin_id]

# --- CCI: Quan et al. 2005 ICD-10 mapping ---
cci_patterns <- list(
  mi          = "^I21|^I22|^I252",
  chf         = "^I099|^I110|^I130|^I132|^I255|^I420|^I42[5-9]|^I43|^I50|^P290",
  pvd         = "^I70|^I71|^I731|^I73[89]|^I771|^I790|^I792|^K55[189]|^Z95[89]",
  cvd         = "^G4[56]|^H340|^I6",
  dementia    = "^F0[0-3]|^F051|^G30|^G311",
  copd        = "^I27[89]|^J4[0-7]|^J6[0-7]|^J684|^J70[13]",
  rheumatic   = "^M0[56]|^M315|^M3[2-4]|^M35[13]|^M360",
  pud         = "^K2[5-8]",
  mild_liver  = "^B18|^K70[0-39]|^K71[3-57]|^K7[34]|^K760|^K76[2-489]|^Z944",
  dm_uncomp   = "^E10[0169]|^E11[0169]|^E1[2-4][0169]",
  dm_comp     = "^E10[2-57]|^E11[2-57]|^E1[2-4][2-57]",
  hemiplegia  = "^G041|^G114|^G80[12]|^G8[12]|^G83[0-49]",
  renal       = "^I120|^I131|^N03[2-7]|^N05[2-7]|^N1[89]|^N250|^Z49[0-2]|^Z940|^Z992",
  malignancy  = "^C[0-2][0-6]|^C3[0-4]|^C3[7-9]|^C4[013]|^C4[5-8]|^C5[0-8]|^C6|^C7[0-6]|^C8[1-58]|^C9[0-7]",
  mod_liver   = "^I85[09]|^I864|^I982|^K704|^K711|^K72[19]|^K76[5-7]",
  metastatic  = "^C7[7-9]|^C80",
  aids        = "^B2[0-24]"
)
cci_weights <- c(mi=1, chf=1, pvd=1, cvd=1, dementia=1, copd=1, rheumatic=1,
                 pud=1, mild_liver=1, dm_uncomp=1, dm_comp=2, hemiplegia=2,
                 renal=2, malignancy=2, mod_liver=3, metastatic=6, aids=6)

# --- Elixhauser: van Walraven weights ---
elix_patterns <- list(
  chf            = "^I099|^I110|^I130|^I132|^I255|^I420|^I42[5-9]|^I43|^I50|^P290",
  arrhythmia     = "^I441|^I44[23]|^I456|^I459|^I47|^I48|^I49|^R00[018]|^T821|^Z450|^Z950",
  valvular       = "^A520|^I0[5-8]|^I091|^I098|^I3[4-9]|^Q23[0-3]|^Z95[2-4]",
  pulm_circ      = "^I26|^I27[0-289]",
  pvd            = "^I70|^I71|^I731|^I73[89]|^I771|^I790|^I792|^K55[189]|^Z95[89]",
  htn_uncomp     = "^I10",
  htn_comp       = "^I1[1-35]",
  paralysis      = "^G041|^G114|^G80[12]|^G8[12]|^G83[0-49]",
  neuro_other    = "^G1[0-3]|^G2[0-2]|^G254|^G255|^G31[289]|^G3[2-7]|^G4[01]|^G931|^G934|^R470|^R56",
  copd           = "^I27[89]|^J4[0-7]|^J6[0-7]|^J684|^J70[13]",
  dm_uncomp      = "^E10[0169]|^E11[0169]|^E1[2-4][0169]",
  dm_comp        = "^E10[2-57]|^E11[2-57]|^E1[2-4][2-57]",
  hypothyroid    = "^E0[0-3]|^E890",
  renal          = "^I120|^I131|^N03[2-7]|^N05[2-7]|^N1[89]|^N250|^Z49[0-2]|^Z940|^Z992",
  liver          = "^B18|^I85|^I864|^I982|^K70|^K71[13-57]|^K7[234]|^K76[02-59]|^Z944",
  pud            = "^K2[5-8]",
  aids           = "^B2[0-24]",
  lymphoma       = "^C8[1-58]|^C96|^C900|^C902",
  metastatic     = "^C7[7-9]|^C80",
  solid_tumor    = "^C[0-2][0-6]|^C3[0-4]|^C3[7-9]|^C4[013]|^C4[5-8]|^C5[0-8]|^C6|^C7[0-6]",
  rheumatic      = "^M0[56]|^M315|^M3[2-4]|^M35[13]|^M360",
  coagulopathy   = "^D6[5-8]|^D69[13-6]",
  obesity        = "^E66",
  weight_loss    = "^E4[0-6]|^R634|^R64",
  fluid_elec     = "^E222|^E86|^E87",
  blood_anemia   = "^D500",
  def_anemia     = "^D50[89]|^D5[1-3]",
  alcohol        = "^F10|^E52|^G621|^I426|^K292|^K70[019]|^T51|^Z502|^Z714|^Z721",
  drug_abuse     = "^F1[1-69]|^Z715|^Z722",
  psychoses      = "^F2[0-9]|^F302|^F312|^F315",
  depression     = "^F204|^F31[34]|^F3[23]|^F341|^F412|^F432"
)
elix_weights <- c(chf=7, arrhythmia=5, valvular=-1, pulm_circ=4, pvd=2,
                  htn_uncomp=0, htn_comp=0, paralysis=7, neuro_other=6, copd=3,
                  dm_uncomp=0, dm_comp=0, hypothyroid=0, renal=5, liver=11,
                  pud=0, aids=0, lymphoma=9, metastatic=12, solid_tumor=4,
                  rheumatic=0, coagulopathy=3, obesity=-4, weight_loss=6,
                  fluid_elec=5, blood_anemia=-2, def_anemia=-2,
                  alcohol=0, drug_abuse=-7, psychoses=0, depression=-3)

# --- Compute scores ---
all_codes <- unique(person_codes$code)
cat(sprintf("  %d unique ICD-10 codes in cohort\n", length(all_codes)))

# Pre-compute which codes match each category
compute_index <- function(patterns, weights) {
  code_cats <- lapply(names(patterns), function(cat) {
    pat <- patterns[[cat]]
    matched <- all_codes[grepl(pat, all_codes)]
    matched
  })
  names(code_cats) <- names(patterns)

  cat_hits <- rbindlist(lapply(names(code_cats), function(cat) {
    codes <- code_cats[[cat]]
    if (length(codes) == 0) return(NULL)
    person_codes[code %in% codes, .(category = cat), by = kojin_id]
  }))
  cat_hits <- unique(cat_hits)

  # Score per person: sum of weights for categories present
  person_score <- cat_hits[, .(score = sum(weights[category])), by = kojin_id]
  person_score
}

cat("  Computing CCI...\n")
cci_scores <- compute_index(cci_patterns, cci_weights)
setnames(cci_scores, "score", "cci")

cat("  Computing Elixhauser...\n")
elix_scores <- compute_index(elix_patterns, elix_weights)
setnames(elix_scores, "score", "elixhauser")

# --- Merge all scores ---
all_patients <- unique(pmap$kojin_id)
result <- data.table(kojin_id = all_patients)
result <- merge(result, n_codes, by = "kojin_id", all.x = TRUE)
result <- merge(result, cci_scores, by = "kojin_id", all.x = TRUE)
result <- merge(result, elix_scores, by = "kojin_id", all.x = TRUE)
result[is.na(n_codes), n_codes := 0L]
result[is.na(cci), cci := 0L]
result[is.na(elixhauser), elixhauser := 0L]

cat(sprintf("  CCI: median=%d, mean=%.1f, max=%d\n",
    median(result$cci), mean(result$cci), max(result$cci)))
cat(sprintf("  Elixhauser: median=%d, mean=%.1f, min=%d, max=%d\n",
    median(result$elixhauser), mean(result$elixhauser),
    min(result$elixhauser), max(result$elixhauser)))
cat(sprintf("  n_codes: median=%d, mean=%.0f, max=%d\n",
    median(result$n_codes), mean(result$n_codes), max(result$n_codes)))

fwrite(result, out_tsv, sep = "\t")
cat(sprintf("  Written: %s (%d patients)\n", out_tsv, nrow(result)))
cat("Done.\n")
