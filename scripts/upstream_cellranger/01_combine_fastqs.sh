#!/usr/bin/env bash
set -euo pipefail

# Check existing pool FASTQs and create the configured RR1316 links.
# Input paths and source-file identities must match the intended dataset.

RR1316_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/raw_data"
SHARED_DIR="/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data"  # RR0104, RR0508, RR0912, RR001

MANIFEST="$SHARED_DIR/combine_fastqs_manifest.sha256"
: > "$MANIFEST"

# verify_exists <path> -- fails loudly rather than letting a missing/empty
# input pass silently into FastQC or cellranger.
verify_exists() {
  local f="$1"
  if [ ! -s "$f" ]; then
    echo "ERROR: missing or empty expected raw FASTQ: $f" >&2
    exit 1
  fi
  sha256sum "$f" >> "$MANIFEST"
  echo "Verified: $f"
}

# RR0104, RR0508, RR0912: already-combined pool-level files.
verify_exists "$SHARED_DIR/RR0104_S1_L001_R1_001.fastq.gz"
verify_exists "$SHARED_DIR/RR0104_S1_L001_R2_001.fastq.gz"
verify_exists "$SHARED_DIR/RR0508_S4_L001_R1_001.fastq.gz"
verify_exists "$SHARED_DIR/RR0508_S4_L001_R2_001.fastq.gz"
verify_exists "$SHARED_DIR/RR0912_S2_L001_R1_001.fastq.gz"
verify_exists "$SHARED_DIR/RR0912_S2_L001_R2_001.fastq.gz"

# RR1316: raw data was delivered as RR001; ensure it exists under the
# RR1316-named path cellranger's fastq_id discovery expects.
mkdir -p "$RR1316_DIR"
for suffix in R1 R2; do
  src="$SHARED_DIR/RR001_S1_L001_${suffix}_001.fastq.gz"
  dst="$RR1316_DIR/RR1316_S1_L001_${suffix}_001.fastq.gz"
  if [ ! -e "$dst" ]; then
    verify_exists "$src"
    ln -s "$src" "$dst"
    echo "Linked $dst -> $src"
  fi
  verify_exists "$dst"
done

echo "Done. Manifest: $MANIFEST"
