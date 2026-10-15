"""
Extract disease-family-specific hospitalizations from DeSC receipt data.

Joins receipt.csv (inpatient) + receipt_diseases.csv (primary diagnosis)
+ m_icd10.csv (ICD-10 mapping) to produce per-family hospitalization files.

Usage:
  python3 src/extract_disease_specific_hospitalization.py \
      /home/desc/utf8/all/receipt.csv \
      /home/desc/utf8/all/receipt_diseases.csv \
      /home/desc/utf8/all/m_icd10.csv \
      config/disease_families.tsv \
      output/
"""

import sys, os, csv

receipt_path = sys.argv[1]
diseases_path = sys.argv[2]
icd10_path = sys.argv[3]
families_path = sys.argv[4]
out_base = sys.argv[5]

# --- Step 1: Load m_icd10.csv (diseases_code -> icd10_code) ---
print("Step 1: Loading m_icd10.csv ...", file=sys.stderr, flush=True)
code2icd = {}
with open(icd10_path, encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)
    di = header.index('diseases_code')
    ii = header.index('icd10_code')
    for row in reader:
        if len(row) > max(di, ii):
            code2icd[row[di]] = row[ii]
print(f"  {len(code2icd):,} disease codes mapped", file=sys.stderr, flush=True)

# --- Step 2: Load disease families ---
print("Step 2: Loading disease families ...", file=sys.stderr, flush=True)
families = {}
with open(families_path, encoding='utf-8') as f:
    header = f.readline().strip().split('\t')
    ni = header.index('family_id')
    ii = header.index('icd10_prefixes')
    for line in f:
        fs = line.strip().split('\t')
        if len(fs) <= max(ni, ii):
            continue
        fid = fs[ni]
        prefixes = tuple(c.strip() for c in fs[ii].split(','))
        families[fid] = prefixes

print(f"  {len(families)} families", file=sys.stderr, flush=True)

# --- Step 3: Scan receipt.csv for inpatient receipts ---
print("Step 3: Scanning receipt.csv for inpatient receipts (nyugai_kbn_code=2) ...",
      file=sys.stderr, flush=True)
# Build: receipt_id -> (kojin_id, nyuin_ymd)
inpatient = {}
n_total = 0
with open(receipt_path, encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)
    ri = header.index('receipt_id')
    ki = header.index('kojin_id')
    ni = header.index('nyugai_kbn_code')
    yi = header.index('nyuin_ymd')
    for row in reader:
        n_total += 1
        if n_total % 100_000_000 == 0:
            print(f"  {n_total:,} rows, {len(inpatient):,} inpatient ...",
                  file=sys.stderr, flush=True)
        if len(row) > max(ri, ki, ni, yi) and row[ni] == '2':
            nyuin = row[yi] if yi < len(row) else ''
            inpatient[row[ri]] = (row[ki], nyuin)

print(f"  Total: {n_total:,} rows, {len(inpatient):,} inpatient receipts",
      file=sys.stderr, flush=True)

# --- Step 4: Scan receipt_diseases.csv ---
print("Step 4: Scanning receipt_diseases.csv (168GB) ...",
      file=sys.stderr, flush=True)

# Per family: kojin_id -> (earliest_admission_ymd, count)
fam_hosp = {fid: {} for fid in families}

n_total = 0
n_matched = 0
with open(diseases_path, encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)
    ri = header.index('receipt_id')
    di = header.index('diseases_code')
    si = header.index('shushoubyomei_flg_code')
    max_idx = max(ri, di, si)

    for row in reader:
        n_total += 1
        if n_total % 100_000_000 == 0:
            print(f"  {n_total:,} rows, {n_matched:,} matched ...",
                  file=sys.stderr, flush=True)

        if len(row) <= max_idx:
            continue
        if row[si] != '1':
            continue

        receipt_id = row[ri]
        if receipt_id not in inpatient:
            continue

        diseases_code = row[di]
        icd10 = code2icd.get(diseases_code, '')
        if not icd10:
            continue

        kojin_id, nyuin_ymd = inpatient[receipt_id]

        for fid, prefixes in families.items():
            if any(icd10.startswith(p) for p in prefixes):
                n_matched += 1
                cur = fam_hosp[fid].get(kojin_id)
                if cur is None:
                    fam_hosp[fid][kojin_id] = [nyuin_ymd, 1]
                else:
                    cur[1] += 1
                    if nyuin_ymd and (not cur[0] or nyuin_ymd < cur[0]):
                        cur[0] = nyuin_ymd

print(f"  Total: {n_total:,} rows, {n_matched:,} family matches",
      file=sys.stderr, flush=True)

# --- Step 5: Write per-family output ---
print("Step 5: Writing output ...", file=sys.stderr, flush=True)
for fid in sorted(families.keys()):
    hosp = fam_hosp[fid]
    odir = os.path.join(out_base, fid, "all_comorbid", "cohort")
    os.makedirs(odir, exist_ok=True)
    outpath = os.path.join(odir, "person_hospitalization_specific.tsv")
    with open(outpath, 'w') as out:
        out.write("kojin_id\tfirst_admission_ymd\tn_admissions\n")
        for kid in sorted(hosp.keys(), key=int):
            adm_ymd, n_adm = hosp[kid]
            out.write(f"{kid}\t{adm_ymd}\t{n_adm}\n")
    print(f"  {fid}: {len(hosp):,} patients with disease-specific admission",
          file=sys.stderr, flush=True)

print("Done.", file=sys.stderr, flush=True)
