"""Convert rowid_window_lookup.tsv to compact binary for fast R access.

Output format: for each old_row_id (1..N), two uint16 values:
  start_encoded = year * 12 + month  (e.g., 2014/04 = 2014*12+4 = 24172)
  end_encoded   = year * 12 + month

Total size: N * 4 bytes. R reads with readBin and indexes by old_row_id.
"""
import struct
import sys
import time

inpath = "data/rowid_window_lookup.tsv"
outpath = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad/window_lookup.bin"

def encode_ym(s):
    return int(s[:4]) * 12 + int(s[5:7])

t0 = time.time()
n = 0
with open(inpath) as f, open(outpath, 'wb') as out:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 2:
            out.write(struct.pack('<HH', 0, 0))
            n += 1
            continue
        window = parts[1]
        sep = window.index('..')
        start_s = window[:sep]
        end_s = window[sep+2:]
        out.write(struct.pack('<HH', encode_ym(start_s), encode_ym(end_s)))
        n += 1
        if n % 50000000 == 0:
            print(f"  {n/1e6:.0f}M rows, {time.time()-t0:.0f}s", flush=True)

t1 = time.time()
print(f"Done: {n} rows in {t1-t0:.0f}s")
print(f"Output: {outpath} ({n*4/1e9:.2f} GB)")
