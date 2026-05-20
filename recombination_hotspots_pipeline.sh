#!/usr/bin/env bash
# =============================================================================
# Germline Recombination Hotspot Pipeline
# Pratto et al. 2014 (GSE59836) × deCODE hg38
#
# GSE59836_Peak_data_Supplementary_File_1.txt.gz format (confirmed):
#   ## comment lines at top describing columns
#   Col 1: chrom   Col 2: start   Col 3: end
#   Col 4: AA1_Strength   Col 5: AA2_Strength
#   Col 6: AB1_Strength   ...
#
# Requirements: wget, python3, bedtools
# =============================================================================
set -euo pipefail

WORKDIR="germline_hotspots"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

SUPPL_BASE="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59836/suppl"
PEAK_FILE="GSE59836_Peak_data_Supplementary_File_1.txt.gz"

echo ""
echo "=================================================="
echo " Germline Recombination Hotspot Pipeline"
echo " Pratto 2014 (GSE59836) × deCODE hg38"
echo "=================================================="

# ── Step 1: liftOver + chain ───────────────────────────────────────────────
echo ""
echo "[1/5] Setting up liftOver..."
if [ ! -f liftOver ]; then
    wget -q --show-progress \
        https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O liftOver
fi
chmod +x liftOver
echo "  OK: $(./liftOver 2>&1 | head -1)"

if [ ! -f hg19ToHg38.over.chain.gz ]; then
    wget -q --show-progress \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi
echo "  OK: hg19ToHg38.over.chain.gz"

# ── Step 2: Download peak data ─────────────────────────────────────────────
echo ""
echo "[2/5] Downloading Pratto 2014 peak data..."
if [ ! -f "$PEAK_FILE" ]; then
    wget -q --show-progress "${SUPPL_BASE}/${PEAK_FILE}" -O "$PEAK_FILE"
fi
echo "  OK: $PEAK_FILE ($(du -h "$PEAK_FILE" | cut -f1))"

# Peek at the full comment block so we know all columns
echo ""
echo "  File header:"
gunzip -c "$PEAK_FILE" | grep '^##' | sed 's/^/    /'
echo ""

# ── Step 3: Parse → BED (hg19) ────────────────────────────────────────────
echo "[3/5] Parsing peak data → BED..."

python3 - <<'PYEOF'
import gzip, sys

INFILE  = "GSE59836_Peak_data_Supplementary_File_1.txt.gz"
OUTFILE = "pratto_hotspots_hg19.bed"

# Per the ## comments:
#   col 0 = chrom, col 1 = start, col 2 = end (0-based indexing)
#   col 3 = AA1_Strength, col 4 = AA2_Strength  (PRDM9-A homozygous individuals)
# We keep ALL hotspots (any individual); filter by AA strength > 0 if desired.

CHR_COL   = 0
START_COL = 1
END_COL   = 2
AA1_COL   = 3   # AA1_Strength
AA2_COL   = 4   # AA2_Strength

n_written  = 0
n_skipped  = 0
n_comments = 0

with gzip.open(INFILE, "rt") as fh, open(OUTFILE, "w") as out:
    for line in fh:
        line = line.rstrip("\n")

        # Skip comment/header lines
        if line.startswith("#") or line.startswith("Chr") or line.startswith("chrom"):
            n_comments += 1
            continue

        parts = line.split("\t")
        if len(parts) < 3:
            n_skipped += 1
            continue

        chrom = parts[CHR_COL].strip()
        if not chrom:
            n_skipped += 1
            continue
        # Ensure chr prefix
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        try:
            start = int(parts[START_COL])
            end   = int(parts[END_COL])
        except ValueError:
            n_skipped += 1
            continue

        if start < 0 or end <= start:
            n_skipped += 1
            continue

        # Use max of AA1+AA2 strength as score (best PRDM9-A signal)
        score = 0
        for col in [AA1_COL, AA2_COL]:
            if len(parts) > col:
                try:
                    v = float(parts[col])
                    if v > score:
                        score = v
                except ValueError:
                    pass

        ucsc_score = min(1000, int(score))
        name = f"hotspot_{chrom}_{start}"

        out.write(f"{chrom}\t{start}\t{end}\t{name}\t{ucsc_score}\t.\t{score:.2f}\n")
        n_written += 1

print(f"  Comment/header lines skipped : {n_comments}")
print(f"  Malformed lines skipped      : {n_skipped}")
print(f"  Hotspot records written      : {n_written:,}")
print(f"  Output                       : {OUTFILE}")

if n_written == 0:
    print("\n  [ERROR] No records written. Run:")
    print("  gunzip -c GSE59836_Peak_data_Supplementary_File_1.txt.gz | head -20")
    sys.exit(1)
PYEOF

sort -k1,1 -k2,2n pratto_hotspots_hg19.bed -o pratto_hotspots_hg19.bed

echo ""
echo "  First 5 records:"
head -5 pratto_hotspots_hg19.bed | \
    awk '{printf "    %-6s  %-10s  %-10s  score=%-6s  strength=%s\n",$1,$2,$3,$5,$7}'

# ── Step 4: LiftOver hg19 → hg38 ──────────────────────────────────────────
echo ""
echo "[4/5] Lifting over hg19 → hg38..."

./liftOver \
    pratto_hotspots_hg19.bed \
    hg19ToHg38.over.chain.gz \
    pratto_hotspots_hg38.bed \
    pratto_unmapped.bed

MAPPED=$(wc -l < pratto_hotspots_hg38.bed)
UNMAPPED=$(grep -vc '^#' pratto_unmapped.bed 2>/dev/null || echo 0)
echo "  Mapped   : $MAPPED  → pratto_hotspots_hg38.bed"
echo "  Unmapped : $UNMAPPED → pratto_unmapped.bed"

# ── Step 5: Fetch deCODE hotspots ─────────────────────────────────────────
echo ""
echo "[5/5] Fetching deCODE hg38 hotspots (>=50 cM/Mb) from UCSC..."

if [ ! -f decode_recombAvg_hotspots_hg38.bed ]; then
python3 - <<'PYEOF'
import urllib.request, json, time, sys

UCSC_API  = "https://api.genome.ucsc.edu"
GENOME    = "hg38"
TRACK     = "recombAvg"
THRESHOLD = 50.0
CHROMS    = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
PAUSE     = 0.3
OUT       = "decode_recombAvg_hotspots_hg38.bed"

def pick(row, candidates):
    for c in candidates:
        if c in row: return c
    return None

all_hot = []
for chrom in CHROMS:
    url = f"{UCSC_API}/getData/track?genome={GENOME}&track={TRACK}&chrom={chrom}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = json.loads(r.read())
    except Exception as e:
        print(f"  [WARN] {chrom}: {e}", file=sys.stderr)
        time.sleep(PAUSE); continue

    rows = None
    for key in [TRACK, "data", GENOME]:
        if key in data and isinstance(data[key], list):
            rows = data[key]; break
    if not rows:
        print(f"  [WARN] No data for {chrom}", file=sys.stderr)
        time.sleep(PAUSE); continue

    sf = pick(rows[0], ["chromStart", "start", "tStart"])
    ef = pick(rows[0], ["chromEnd",   "end",   "tEnd"])
    rf = pick(rows[0], ["value", "score", "decodeAvg"])

    if not all([sf, ef, rf]):
        print(f"  [ERROR] Unknown fields in {chrom}: {list(rows[0].keys())}", file=sys.stderr)
        sys.exit(1)

    hot = []
    for r in rows:
        rate = float(r.get(rf, 0) or 0)
        if rate >= THRESHOLD:
            hot.append([str(r["chrom"]), int(r[sf]), int(r[ef]), rate])

    # merge adjacent bins
    merged = []
    for entry in hot:
        if merged and entry[0] == merged[-1][0] and entry[1] <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], entry[2])
            merged[-1][3] = max(merged[-1][3], entry[3])
        else:
            merged.append(list(entry))

    all_hot.extend(merged)
    print(f"  {chrom}: {len(rows):>7,} bins → {len(merged):>4,} hotspots")
    time.sleep(PAUSE)

with open(OUT, "w") as fh:
    fh.write('track name="deCODE_hotspots_50cMpMb" '
             'description="deCODE hg38 >=50cM/Mb (Halldorsson 2019)" useScore=1\n')
    for chrom, s, e, r in all_hot:
        fh.write(f"{chrom}\t{s}\t{e}\thot_{chrom}_{s}"
                 f"\t{min(1000,int(r*10))}\t.\t{r:.4f}\n")

print(f"\n  Written {len(all_hot):,} intervals → {OUT}")
PYEOF
else
    echo "  Already present — skipping."
fi

# ── Step 6: Intersect ──────────────────────────────────────────────────────
echo ""
echo "[6/6] Intersecting Pratto DSBs ∩ deCODE crossover hotspots..."

if ! command -v bedtools &>/dev/null; then
    echo "  [ERROR] bedtools not found. Install: sudo apt install bedtools"
    exit 1
fi

grep -v '^track\|^browser' pratto_hotspots_hg38.bed \
    | awk 'NF>=3' | sort -k1,1 -k2,2n > pratto_clean.bed

grep -v '^track\|^browser' decode_recombAvg_hotspots_hg38.bed \
    | awk 'NF>=3' | sort -k1,1 -k2,2n > decode_clean.bed

bedtools intersect \
    -a pratto_clean.bed \
    -b decode_clean.bed \
    > final_germline_hotspots_hg38.bed

FINAL=$(wc -l < final_germline_hotspots_hg38.bed)

echo ""
echo "=================================================="
echo " DONE"
echo "=================================================="
echo ""
echo "  Working directory: $(pwd)"
echo ""
printf "  %-46s  %s\n" "pratto_hotspots_hg19.bed"            "All Pratto DSB hotspots (hg19)"
printf "  %-46s  %s\n" "pratto_hotspots_hg38.bed"            "All Pratto DSB hotspots (hg38)"
printf "  %-46s  %s\n" "decode_recombAvg_hotspots_hg38.bed"  "deCODE crossover hotspots (>=50 cM/Mb)"
printf "  %-46s  %s\n" "final_germline_hotspots_hg38.bed"    "★ High-confidence intersection"
echo ""
echo "  Final high-confidence hotspot count: $FINAL"
echo ""
echo "  Columns in final BED:"
echo "    chrom  start  end  name  ucsc_score  strand  AA_strength"
echo ""
echo "  These sites have BOTH:"
echo "    ✓ Active meiotic DSBs in real testis (Pratto 2014, DMC1-SSDS)"
echo "    ✓ Elevated crossover rate in the population (Halldorsson 2019)"
echo ""