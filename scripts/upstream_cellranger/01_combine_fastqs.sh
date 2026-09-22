#!/usr/bin/env bash
set -euo pipefail

# 01_combine_fastqs.sh — Upstream raw-data processing, step 1 of 4
# Concatenates the per-lane (L007/L008) FASTQs for each of the 4 GEM-reaction
# pools (RR0104, RR0508, RR0912, RR1316) into single per-pool R1/R2 files,
# ready for FastQC (02) and `cellranger multi` (04).
# Paths/filenames are HPCC-specific to the GSE311941 run; edit to reproduce.
#
# Output directories/filenames are chosen to exactly match what
# 03_create_multiconfig_csv.sh's `[libraries]` `fastqs` field points at for
# each pool, and follow the `<fastq_id>_S#_L001_R#_001.fastq.gz` naming
# `cellranger multi` auto-discovers via the `fastq_id` column -- so 04's
# `cellranger multi` runs can find these files without any extra config.
# Combined outputs are written into the SAME directory as the raw per-lane
# inputs (not a separate "combined" subfolder): this is safe because the raw
# per-lane filenames (e.g. RR13-16_S3_L007_R1_001.fastq.gz, hyphenated) never
# share a `<fastq_id>_S` prefix with the combined outputs (e.g.
# RR1316_S3_L001_R1_001.fastq.gz, unhyphenated), so cellranger's fastq_id
# prefix match cannot pick up the raw lane files by mistake.
# 02_fastqc.sh must be kept in sync with the exact filenames/directories
# written here.
#
# Previously this script ran with no `set -euo pipefail`, no check that
# inputs existed, and no checksum to confirm what was written. A missing lane
# file would make `cat` print an error to stderr and continue, silently
# producing a truncated "combined" FASTQ that downstream steps would process
# without any indication anything was wrong. This version fails loudly on any
# missing/empty input and records a sha256 manifest of every output file.

# Same directory serves as both raw-lane input and cellranger fastqs target
RR1316_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/raw_data"
SHARED_DIR="/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data"  # RR0104, RR0508, RR0912
mkdir -p "$RR1316_DIR" "$SHARED_DIR"

MANIFEST="$SHARED_DIR/combine_fastqs_manifest.sha256"
: > "$MANIFEST"

# combine <dir> <out_name> <in1> <in2> [in3 ...]
combine() {
  local dir="$1" out="$2"; shift 2
  local inputs=("$@")
  for f in "${inputs[@]}"; do
    if [ ! -s "$dir/$f" ]; then
      echo "ERROR: missing or empty input file: $dir/$f" >&2
      exit 1
    fi
  done
  local out_path="$dir/$out"
  cat "${inputs[@]/#/$dir/}" > "$out_path"
  sha256sum "$out_path" >> "$MANIFEST"
  echo "Wrote $out_path (from ${#inputs[@]} input file(s))"
}

combine "$RR1316_DIR" RR1316_S3_L001_R1_001.fastq.gz \
  RR13-16_S3_L007_R1_001.fastq.gz RR13-16_S3_L008_R1_001.fastq.gz
combine "$RR1316_DIR" RR1316_S3_L001_R2_001.fastq.gz \
  RR13-16_S3_L007_R2_001.fastq.gz RR13-16_S3_L008_R2_001.fastq.gz

combine "$SHARED_DIR" RR0104_combined_R1_001.fastq.gz \
  RR1-4_S1_L007_R1_001.fastq.gz RR1-4_S1_L008_R1_001.fastq.gz
combine "$SHARED_DIR" RR0104_combined_R2_001.fastq.gz \
  RR1-4_S1_L007_R2_001.fastq.gz RR1-4_S1_L008_R2_001.fastq.gz

combine "$SHARED_DIR" RR0508_S4_L001_R1_001.fastq.gz \
  RR5-8_S4_L007_R1_001.fastq.gz RR5-8_S4_L008_R1_001.fastq.gz
combine "$SHARED_DIR" RR0508_S4_L001_R2_001.fastq.gz \
  RR5-8_S4_L007_R2_001.fastq.gz RR5-8_S4_L008_R2_001.fastq.gz

combine "$SHARED_DIR" RR0912_S2_L001_R1_001.fastq.gz \
  RR9-12_S2_L007_R1_001.fastq.gz RR9-12_S2_L008_R1_001.fastq.gz
combine "$SHARED_DIR" RR0912_S2_L001_R2_001.fastq.gz \
  RR9-12_S2_L007_R2_001.fastq.gz RR9-12_S2_L008_R2_001.fastq.gz

# RR0104 received an additional top-up sequencing run (RR001_S1_L001_*): merge
# it into the already-lane-combined RR0104 files above.
if [ ! -s "$SHARED_DIR/RR0104_combined_R1_001.fastq.gz" ] || [ ! -s "$SHARED_DIR/RR0104_combined_R2_001.fastq.gz" ]; then
  echo "ERROR: RR0104_combined_* outputs missing from $SHARED_DIR (earlier step must have failed)" >&2
  exit 1
fi
for f in RR001_S1_L001_R1_001.fastq.gz RR001_S1_L001_R2_001.fastq.gz; do
  if [ ! -s "$SHARED_DIR/$f" ]; then
    echo "ERROR: missing or empty input file: $SHARED_DIR/$f" >&2
    exit 1
  fi
done

cat "$SHARED_DIR/RR0104_combined_R1_001.fastq.gz" "$SHARED_DIR/RR001_S1_L001_R1_001.fastq.gz" \
  > "$SHARED_DIR/RR0104_S1_L001_R1_001.fastq.gz"
sha256sum "$SHARED_DIR/RR0104_S1_L001_R1_001.fastq.gz" >> "$MANIFEST"
cat "$SHARED_DIR/RR0104_combined_R2_001.fastq.gz" "$SHARED_DIR/RR001_S1_L001_R2_001.fastq.gz" \
  > "$SHARED_DIR/RR0104_S1_L001_R2_001.fastq.gz"
sha256sum "$SHARED_DIR/RR0104_S1_L001_R2_001.fastq.gz" >> "$MANIFEST"
echo "Wrote $SHARED_DIR/RR0104_S1_L001_R{1,2}_001.fastq.gz (RR0104_combined + top-up run)"

echo "Done. Manifest: $MANIFEST"
