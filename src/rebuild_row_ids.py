"""
Build person_map (new_row_id -> old_row_id -> kojin_id) for all families.
Single pass through COO + scan of row_id_number_small.txt.

Usage:
  python3 src/rebuild_row_ids.py \
      data/coo.txt.zst data/col_id_disease_name_small.txt \
      config/disease_families.tsv data/row_id_number_small.txt output/
"""

import struct, sys, os, subprocess
import numpy as np

coo_zst = sys.argv[1]
col_ids_path = sys.argv[2]
families_path = sys.argv[3]
row_ids_orig_path = sys.argv[4]
out_base = sys.argv[5]

# --- Parse families ---
families = {}
with open(families_path, encoding='utf-8') as f:
    header = f.readline().strip().split('\t')
    ni, ii = header.index('family_id'), header.index('icd10_prefixes')
    for line in f:
        fs = line.strip().split('\t')
        if len(fs) <= max(ni, ii):
            continue
        families[fs[ni]] = [c.strip() for c in fs[ii].split(',')]

# Skip families that already have person_map
need = {}
for fid, prefixes in families.items():
    pmap = os.path.join(out_base, fid, "all_comorbid", "cohort", "person_map.tsv")
    if not os.path.exists(pmap):
        need[fid] = prefixes

if not need:
    print("All person_maps exist, nothing to do", file=sys.stderr)
    sys.exit(0)

print(f"Building person_map for {len(need)} families", file=sys.stderr)

# --- Parse col IDs ---
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

fam_col_arrays = {}
for fid, prefixes in need.items():
    s = set()
    for cid, code in col_codes.items():
        if any(code.startswith(p) for p in prefixes):
            s.add(cid)
    fam_col_arrays[fid] = np.array(sorted(s), dtype=np.int32)
    print(f"  {fid}: {len(s)} index cols", file=sys.stderr)

# --- Pass: identify rows ---
proc = subprocess.Popen(['zstd', '-dc', coo_zst], stdout=subprocess.PIPE)
hdr = proc.stdout.read(8)
nrows = struct.unpack('<i', hdr[:4])[0]
print(f"COO: {nrows:,} rows", file=sys.stderr)

masks = {f: np.zeros(nrows + 1, dtype=np.bool_) for f in need}
BUF = 8 * 1024 * 1024
left = b''
cnt = 0
while True:
    chunk = proc.stdout.read(BUF)
    if not chunk:
        break
    data = left + chunk
    np_ = len(data) // 8
    rem = len(data) % 8
    pairs = np.frombuffer(data[:np_ * 8], dtype=np.int32).reshape(-1, 2)
    left = data[np_ * 8:] if rem else b''
    rows, cols = pairs[:, 0], pairs[:, 1]
    for fid, arr in fam_col_arrays.items():
        hits = np.isin(cols, arr)
        if np.any(hits):
            masks[fid][rows[hits]] = True
    cnt += np_
    if cnt % 500_000_000 < np_:
        print(f"  {cnt:,}", file=sys.stderr)
proc.wait()

# Write row_ids.txt and collect needed old_ids
all_needed = set()
fam_old_ids = {}
for fid in need:
    old_ids = np.where(masks[fid])[0]
    old_ids = old_ids[old_ids > 0]
    odir = os.path.join(out_base, fid, "all_comorbid", "cohort")
    os.makedirs(odir, exist_ok=True)
    path = os.path.join(odir, "row_ids.txt")
    with open(path, 'w') as f:
        for new_id, oid in enumerate(old_ids, 1):
            f.write(f"{oid},{new_id}\n")
    fam_old_ids[fid] = old_ids
    all_needed.update(old_ids.tolist())
    print(f"  {fid}: {len(old_ids):,} rows", file=sys.stderr)

del masks
print(f"Total unique old_ids: {len(all_needed):,}", file=sys.stderr)

# --- Map old_id -> kojin_id ---
print("Reading row_id_number_small.txt ...", file=sys.stderr)
old2kojin = {}
with open(row_ids_orig_path) as f:
    for ln_num, line in enumerate(f, 1):
        if ln_num not in all_needed:
            continue
        parts = line.strip().split()
        pp = parts[0]
        kojin_id = pp[pp.rfind('_') + 1:]
        old2kojin[ln_num] = kojin_id
        if len(old2kojin) % 5_000_000 == 0:
            print(f"  mapped {len(old2kojin):,}/{len(all_needed):,}",
                  file=sys.stderr)
        if len(old2kojin) == len(all_needed):
            break

print(f"Mapped {len(old2kojin):,} old_ids", file=sys.stderr)

# Write person_map
for fid in need:
    odir = os.path.join(out_base, fid, "all_comorbid", "cohort")
    path = os.path.join(odir, "person_map.tsv")
    with open(path, 'w') as f:
        f.write("new_row_id\told_row_id\tkojin_id\n")
        for new_id, oid in enumerate(fam_old_ids[fid], 1):
            kojin = old2kojin.get(int(oid), "NA")
            f.write(f"{new_id}\t{oid}\t{kojin}\n")
    print(f"  {fid}: person_map -> {path}", file=sys.stderr)

print("Done", file=sys.stderr)
