"""
Extract disease-family-specific COO matrices from the full DeSC COO.
Handles arbitrary number of families from disease_families.tsv.

Pass 1: Identify rows for all families simultaneously (~17GB for 51 families)
Pass 2: Write COOs in batches of BATCH_SIZE families to limit memory

Usage:
  python3 src/extract_cohort_coo.py \
      data/coo.txt.zst data/col_id_disease_name_small.txt \
      config/disease_families.tsv output/
"""

import struct, sys, os, subprocess
import numpy as np

coo_zst_path = sys.argv[1]
col_ids_path = sys.argv[2]
families_path = sys.argv[3]
out_base = sys.argv[4]

BATCH_SIZE = 10  # families per pass-2 batch

# --- Parse families ---
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
        codes = [c.strip() for c in fs[ii].split(',')]
        families[fid] = codes

fam_list = sorted(families.keys())
n_fam = len(fam_list)
print(f"Families: {n_fam}", file=sys.stderr)

# --- Parse column IDs ---
col_codes = {}
with open(col_ids_path, encoding='utf-8') as f:
    for line in f:
        line = line.strip().replace('﻿', '')
        if not line:
            continue
        commas = [i for i, c in enumerate(line) if c == ',']
        if len(commas) < 2:
            continue
        code = line[:commas[0]].strip()
        cid = int(line[commas[0]+1:commas[1]].strip())
        col_codes[cid] = code

NCOLS = max(col_codes.keys())

# --- Build column sets per family ---
fam_index_cols = {}
fam_keep_cols = {}
fam_col_remap = {}
all_col_ids = set(col_codes.keys())

for fid, prefixes in families.items():
    idx = set()
    for cid, code in col_codes.items():
        if any(code.startswith(p) for p in prefixes):
            idx.add(cid)
    fam_index_cols[fid] = idx
    kept = sorted(all_col_ids - idx)
    fam_keep_cols[fid] = kept
    remap = np.zeros(NCOLS + 1, dtype=np.int32)
    for new_id, old_id in enumerate(kept, 1):
        remap[old_id] = new_id
    fam_col_remap[fid] = remap
    print(f"  {fid}: {len(idx)} index cols, {len(kept)} kept", file=sys.stderr)

# --- Column arrays for vectorized matching ---
fam_col_arrays = {}
for fid, idx_cols in fam_index_cols.items():
    fam_col_arrays[fid] = np.array(sorted(idx_cols), dtype=np.int32)

# --- Helper ---
def open_coo(path):
    proc = subprocess.Popen(['zstd', '-dc', path], stdout=subprocess.PIPE)
    hdr = proc.stdout.read(8)
    nrows, ncols = struct.unpack('<ii', hdr)
    return proc, proc.stdout, nrows, ncols

BUF = 8 * 1024 * 1024

# ===============================================================
# Pass 1: Identify rows per family (all families simultaneously)
# ===============================================================
print(f"\n=== Pass 1: identifying rows for {n_fam} families ===",
      file=sys.stderr)
proc1, stream1, nrows, ncols_hdr = open_coo(coo_zst_path)
print(f"COO: {nrows:,} rows x {ncols_hdr} cols", file=sys.stderr)

# Skip families whose COO already exists
skip_extract = set()
for fid in fam_list:
    odir = os.path.join(out_base, fid, "all_comorbid", "cohort")
    if os.path.exists(os.path.join(odir, "feature_matrix.coo.zst")):
        skip_extract.add(fid)
if skip_extract:
    print(f"  Skipping {len(skip_extract)} families with existing COOs",
          file=sys.stderr)

need_extract = [f for f in fam_list if f not in skip_extract]
if not need_extract:
    print("All COOs exist, nothing to extract", file=sys.stderr)
    proc1.kill()
    proc1.wait()
    sys.exit(0)

row_masks = {f: np.zeros(nrows + 1, dtype=np.bool_) for f in need_extract}

count = 0
left = b''
while True:
    chunk = stream1.read(BUF)
    if not chunk:
        break
    data = left + chunk
    np_ = len(data) // 8
    rem = len(data) % 8
    pairs = np.frombuffer(data[:np_ * 8], dtype=np.int32).reshape(-1, 2)
    left = data[np_ * 8:] if rem else b''
    rows, cols = pairs[:, 0], pairs[:, 1]
    for fid in need_extract:
        hits = np.isin(cols, fam_col_arrays[fid])
        if np.any(hits):
            row_masks[fid][rows[hits]] = True
    count += np_
    if count % 500_000_000 < np_:
        print(f"  {count:,} entries ...", file=sys.stderr)
proc1.wait()

for fid in need_extract:
    n = int(np.sum(row_masks[fid]))
    print(f"  {fid}: {n:,} rows", file=sys.stderr)

# ===============================================================
# Pass 2: Write COOs in batches
# ===============================================================
batches = [need_extract[i:i+BATCH_SIZE]
           for i in range(0, len(need_extract), BATCH_SIZE)]
print(f"\n=== Pass 2: {len(batches)} batches of ≤{BATCH_SIZE} ===",
      file=sys.stderr)

for bi, batch in enumerate(batches):
    print(f"\n--- Batch {bi+1}/{len(batches)}: {batch} ---", file=sys.stderr)

    # Compute rank arrays for this batch
    batch_rank = {}
    batch_nrows = {}
    for fid in batch:
        rank = np.cumsum(row_masks[fid]).astype(np.int32)
        batch_rank[fid] = rank
        batch_nrows[fid] = int(rank[nrows])
        print(f"  {fid}: {batch_nrows[fid]:,} rows", file=sys.stderr)

    # Open output files
    out_files = {}
    out_raw = {}
    for fid in batch:
        odir = os.path.join(out_base, fid, "all_comorbid", "cohort")
        os.makedirs(odir, exist_ok=True)
        path = os.path.join(odir, "feature_matrix.coo")
        out_raw[fid] = path
        f = open(path, 'wb')
        n_kept = len(fam_keep_cols[fid])
        f.write(struct.pack('<ii', batch_nrows[fid], n_kept))
        out_files[fid] = f

    # Stream COO
    proc2, stream2, _, _ = open_coo(coo_zst_path)
    nnz = {fid: 0 for fid in batch}
    count2 = 0
    left = b''
    while True:
        chunk = stream2.read(BUF)
        if not chunk:
            break
        data = left + chunk
        np_ = len(data) // 8
        rem = len(data) % 8
        pairs = np.frombuffer(data[:np_ * 8], dtype=np.int32).reshape(-1, 2)
        left = data[np_ * 8:] if rem else b''
        rows, cols = pairs[:, 0], pairs[:, 1]
        for fid in batch:
            rank = batch_rank[fid]
            col_remap = fam_col_remap[fid]
            cur = rank[rows]
            prev = rank[rows - 1]
            row_in = cur > prev
            new_cols = col_remap[cols]
            col_in = new_cols > 0
            mask = row_in & col_in
            n_hit = int(np.sum(mask))
            if n_hit > 0:
                out_data = np.column_stack(
                    [cur[mask], new_cols[mask]]).astype(np.int32)
                out_files[fid].write(out_data.tobytes())
                nnz[fid] += n_hit
        count2 += np_
        if count2 % 500_000_000 < np_:
            print(f"  {count2:,} entries ...", file=sys.stderr)
    proc2.wait()

    # Close, compress, write metadata
    for fid in batch:
        out_files[fid].close()
        raw = out_raw[fid]
        zst = raw + ".zst"
        print(f"  {fid}: {nnz[fid]:,} nnz, compressing ...", file=sys.stderr)
        subprocess.run(['zstd', '-f', '--rm', raw, '-o', zst], check=True)
        odir = os.path.dirname(raw)
        with open(os.path.join(odir, "matrix_stats.tsv"), 'w') as sf:
            sf.write(f"disease_family\t{fid}\n")
            sf.write(f"preproc\tall_comorbid\n")
            sf.write(f"n_patients\t{batch_nrows[fid]}\n")
            sf.write(f"n_features\t{len(fam_keep_cols[fid])}\n")
            sf.write(f"nnz\t{nnz[fid]}\n")
        with open(os.path.join(odir, "col_ids.txt"), 'w', encoding='utf-8') as cf:
            for old_id in fam_keep_cols[fid]:
                new_id = int(fam_col_remap[fid][old_id])
                cf.write(f"{col_codes[old_id]},{new_id}\n")
    del batch_rank

del row_masks

print(f"\n=== All extractions complete ===", file=sys.stderr)
for fid in need_extract:
    stats = os.path.join(out_base, fid, "all_comorbid", "cohort",
                         "matrix_stats.tsv")
    if os.path.exists(stats):
        with open(stats) as f:
            for line in f:
                if line.startswith("n_patients"):
                    n = line.strip().split('\t')[1]
                    print(f"  {fid}: {n} rows", file=sys.stderr)
