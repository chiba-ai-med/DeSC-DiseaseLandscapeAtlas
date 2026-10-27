"""Convert rowid_window_lookup.tsv to compact binary for fast R access.

Output format: for each old_row_id (1..N), two uint16 values:
  start_encoded = year * 12 + (month - 1)  (0-indexed month)
  end_encoded   = year * 12 + (month - 1)

R decode: yr = enc %/% 12L, mo = enc %% 12L + 1L

Total size: N * 4 bytes. R reads with readBin and indexes by old_row_id.
"""
import struct
import sys
import time

inpath = sys.argv[1] if len(sys.argv) > 1 else "data/rowid_window_lookup.tsv"
outpath = sys.argv[2] if len(sys.argv) > 2 else "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad/window_lookup.bin"

def encode_ym(s):
    """Encode 'YYYY/MM' as uint16: year * 12 + (month - 1)."""
    return int(s[:4]) * 12 + (int(s[5:7]) - 1)

def decode_ym(enc):
    """Decode uint16 back to 'YYYY/MM'."""
    yr = enc // 12
    mo = enc % 12 + 1
    return f"{yr:04d}/{mo:02d}"

# Round-trip test for all 12 months
for m in range(1, 13):
    for y in [2014, 2020, 2022]:
        s = f"{y:04d}/{m:02d}"
        enc = encode_ym(s)
        dec = decode_ym(enc)
        assert dec == s, f"Round-trip failed: {s} -> {enc} -> {dec}"
print("Round-trip test passed for all 12 months x 3 years")

t0 = time.time()
n = 0
n_dec = 0
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
        se = encode_ym(start_s)
        ee = encode_ym(end_s)
        out.write(struct.pack('<HH', se, ee))
        # Spot-check December windows
        if end_s.endswith('/12'):
            n_dec += 1
            if n_dec <= 3:
                assert decode_ym(ee) == end_s, f"December decode failed: {end_s} -> {ee} -> {decode_ym(ee)}"
        n += 1
        if n % 50000000 == 0:
            print(f"  {n/1e6:.0f}M rows, {time.time()-t0:.0f}s", flush=True)

t1 = time.time()
print(f"Done: {n} rows in {t1-t0:.0f}s ({n_dec} December windows)")
print(f"Output: {outpath} ({n*4/1e9:.2f} GB)")

# Verify: read back first and last few entries and check against source
print("\nVerification: reading back binary and comparing to source...")
import os
expected_size = n * 4
actual_size = os.path.getsize(outpath)
assert actual_size == expected_size, f"Size mismatch: {actual_size} != {expected_size}"

with open(inpath) as f, open(outpath, 'rb') as binf:
    for i, line in enumerate(f):
        if i >= 10 and i < n - 10:
            continue
        parts = line.rstrip('\n').split('\t')
        binf.seek(i * 4)
        buf = struct.unpack('<HH', binf.read(4))
        expected_row_id = i + 1
        actual_row_id = int(parts[0])
        assert actual_row_id == expected_row_id, \
            f"Line {i}: old_row_id={actual_row_id}, expected {expected_row_id}"
        if len(parts) < 2:
            assert buf == (0, 0), f"Row {i+1}: expected (0,0), got {buf}"
        else:
            window = parts[1]
            sep = window.index('..')
            es = encode_ym(window[:sep])
            ee = encode_ym(window[sep+2:])
            assert buf == (es, ee), f"Row {i+1}: expected ({es},{ee}), got {buf}"
    print(f"Verified first 10 and last 10 rows (old_row_id alignment + values)")

print("All verification passed")
