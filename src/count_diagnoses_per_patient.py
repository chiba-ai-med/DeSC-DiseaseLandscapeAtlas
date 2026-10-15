"""Count unique diagnosis codes per patient from binary COO matrix.

COO format: 8-byte header + pairs of int32 (row_id, col_id).
For each row, count distinct columns (= ICD codes).
Then map to patient via person_map (latest period per kojin_id).

Usage:
    python src/count_diagnoses_per_patient.py <disease_family>
"""
import sys
import struct
import subprocess
from collections import defaultdict

fid = sys.argv[1]

coo_path = f"output/{fid}/all_comorbid/cohort/feature_matrix.coo.zst"
person_map_path = f"output/{fid}/all_comorbid/cohort/person_map.tsv"
out_path = f"output/{fid}/all_comorbid/cohort/n_diagnoses.tsv"

print(f"Counting diagnoses for {fid}...")

row_ncols = defaultdict(int)
proc = subprocess.Popen(["zstd", "-d", "-c", coo_path], stdout=subprocess.PIPE)

# Skip 8-byte header
proc.stdout.read(8)

n = 0
buf_size = 8 * 100000  # Read 100k entries at a time
while True:
    buf = proc.stdout.read(buf_size)
    if not buf:
        break
    n_entries = len(buf) // 8
    for k in range(n_entries):
        i, j = struct.unpack_from('<ii', buf, k * 8)
        row_ncols[i] += 1
    n += n_entries
    if n % 10_000_000 < 100000:
        print(f"  {n:,} entries processed, {len(row_ncols):,} unique rows...")

proc.wait()
print(f"  Total: {n:,} entries, {len(row_ncols):,} unique rows")

# Read person_map
kojin_rows = defaultdict(list)
with open(person_map_path) as f:
    f.readline()  # header
    for line in f:
        parts = line.strip().split('\t')
        new_row_id = int(parts[0])
        old_row_id = int(parts[1])
        kojin_id = parts[2]
        kojin_rows[kojin_id].append((old_row_id, new_row_id))

with open(out_path, 'w') as f:
    f.write("kojin_id\tn_diagnoses\n")
    for kojin_id, entries in kojin_rows.items():
        latest = max(entries, key=lambda x: x[0])
        new_row_id = latest[1]
        n_diag = row_ncols.get(new_row_id, 0)
        f.write(f"{kojin_id}\t{n_diag}\n")

print(f"Output: {out_path}")
print(f"  {len(kojin_rows):,} patients")
