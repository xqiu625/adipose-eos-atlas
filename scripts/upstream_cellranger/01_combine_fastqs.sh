#!/usr/bin/env bash
set -euo pipefail

# 01_combine_fastqs.sh — Upstream raw-data processing, step 1 of 4
# Concatenates the per-lane (L007/L008) FASTQs for each of the 4 GEM-reaction
# pools (RR0104, RR0508, RR0912, RR1316) into single per-pool R1/R2 files,
# ready for FastQC (02) and `cellranger multi` (04).
# Paths/filenames are HPCC-specific to the GSE311941 run; edit to reproduce.
#
# Previously this script ran with no `set -euo pipefail`, no check that inputs
# existed, no explicit working directory (relied on the caller's CWD, and
# wrote outputs into the same directory as the raw inputs), and no checksum
# to confirm what was written. A missing lane file would make `cat` print an
# error to stderr and continue, silently producing a truncated "combined"
# FASTQ that downstream steps would process without any indication anything
# was wrong. This version fails loudly on any missing/empty input, writes
# outputs to a separate directory, and records a sha256 manifest of every
# output file.

RAW_DIR="${RAW_DIR:-$(pwd)}"       # directory containing the raw per-lane FASTQs
OUT_DIR="${OUT_DIR:-$RAW_DIR/combined}"  # combined FASTQs are written here, not into RAW_DIR
mkdir -p "$OUT_DIR"

MANIFEST="$OUT_DIR/combine_fastqs_manifest.sha256"
: > "$MANIFEST"

# combine <out_name> <in1> <in2> [in3 ...]
combine() {
  local out="$1"; shift
  local inputs=("$@")
  for f in "${inputs[@]}"; do
    if [ ! -s "$RAW_DIR/$f" ]; then
      echo "ERROR: missing or empty input file: $RAW_DIR/$f" >&2
      exit 1
    fi
  done
  local out_path="$OUT_DIR/$out"
  cat "${inputs[@]/#/$RAW_DIR/}" > "$out_path"
  sha256sum "$out_path" >> "$MANIFEST"
  echo "Wrote $out_path (from ${#inputs[@]} input file(s))"
}

combine RR1316_S3_L001_R1_001.fastq.gz \
  RR13-16_S3_L007_R1_001.fastq.gz RR13-16_S3_L008_R1_001.fastq.gz
combine RR1316_S3_L001_R2_001.fastq.gz \
  RR13-16_S3_L007_R2_001.fastq.gz RR13-16_S3_L008_R2_001.fastq.gz

combine RR0104_combined_R1_001.fastq.gz \
  RR1-4_S1_L007_R1_001.fastq.gz RR1-4_S1_L008_R1_001.fastq.gz
combine RR0104_combined_R2_001.fastq.gz \
  RR1-4_S1_L007_R2_001.fastq.gz RR1-4_S1_L008_R2_001.fastq.gz

combine RR0508_S4_L001_R1_001.fastq.gz \
  RR5-8_S4_L007_R1_001.fastq.gz RR5-8_S4_L008_R1_001.fastq.gz
combine RR0508_S4_L001_R2_001.fastq.gz \
  RR5-8_S4_L007_R2_001.fastq.gz RR5-8_S4_L008_R2_001.fastq.gz

combine RR0912_S2_L001_R1_001.fastq.gz \
  RR9-12_S2_L007_R1_001.fastq.gz RR9-12_S2_L008_R1_001.fastq.gz
combine RR0912_S2_L001_R2_001.fastq.gz \
  RR9-12_S2_L007_R2_001.fastq.gz RR9-12_S2_L008_R2_001.fastq.gz

# RR0104 received an additional top-up sequencing run (RR001_S1_L001_*): merge
# it into the already-lane-combined RR0104 files above. The RR0104_combined_*
# inputs here are this script's own OUT_DIR outputs from above, not RAW_DIR.
if [ ! -s "$OUT_DIR/RR0104_combined_R1_001.fastq.gz" ] || [ ! -s "$OUT_DIR/RR0104_combined_R2_001.fastq.gz" ]; then
  echo "ERROR: RR0104_combined_* outputs missing from $OUT_DIR (earlier step must have failed)" >&2
  exit 1
fi
for f in RR001_S1_L001_R1_001.fastq.gz RR001_S1_L001_R2_001.fastq.gz; do
  if [ ! -s "$RAW_DIR/$f" ]; then
    echo "ERROR: missing or empty input file: $RAW_DIR/$f" >&2
    exit 1
  fi
done

cat "$OUT_DIR/RR0104_combined_R1_001.fastq.gz" "$RAW_DIR/RR001_S1_L001_R1_001.fastq.gz" \
  > "$OUT_DIR/RR0104_S1_L001_R1_001.fastq.gz"
sha256sum "$OUT_DIR/RR0104_S1_L001_R1_001.fastq.gz" >> "$MANIFEST"
cat "$OUT_DIR/RR0104_combined_R2_001.fastq.gz" "$RAW_DIR/RR001_S1_L001_R2_001.fastq.gz" \
  > "$OUT_DIR/RR0104_S1_L001_R2_001.fastq.gz"
sha256sum "$OUT_DIR/RR0104_S1_L001_R2_001.fastq.gz" >> "$MANIFEST"
echo "Wrote $OUT_DIR/RR0104_S1_L001_R{1,2}_001.fastq.gz (RR0104_combined + top-up run)"

echo "Done. Manifest: $MANIFEST"
