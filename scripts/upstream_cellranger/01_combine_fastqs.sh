#!/usr/bin/env bash
set -euo pipefail

# 01_combine_fastqs.sh — Upstream raw-data processing, step 1 of 4
# Establishes and verifies the per-pool raw FASTQ inputs to FastQC (02) and
# `cellranger multi` (04).
# Paths/filenames are HPCC-specific to the GSE311941 run; edit to reproduce.
#
# ==== RR1316 provenance correction (2026-09-22) =============================
# Per the reviewed, GEO-accepted submission package
# (Desktop/302025/submission/GEO_submission/: README_raw_data.txt,
# FILE_RENAMING_PLAN.txt, CORRECTION_ALL_FILES_PRESENT.txt,
# FILE_COPYING_COMMANDS.txt -- all consistent), RR1316's entire raw FASTQ data
# was delivered by the sequencing core under the filename "RR001", not as
# separate RR13-16_S3_L007/L008 lane files to be combined here. Those lane
# files do not exist anywhere on the current filesystem; only historical
# FastQC reports for an old "RR13-16_S3_combined" name survive (no live
# .fastq.gz for it). RR001 is NOT part of RR0104 -- no GEO document ties the
# two together, and RR0104's own combined file is self-contained. A previous
# version of this script incorrectly treated RR001 as an "RR0104 top-up run"
# concatenated onto RR0104's lanes; that step has been removed. The corrected,
# documented handling is simply: RR001_S1_L001_{R1,R2}_001.fastq.gz IS
# RR1316's raw data, made available at RR1316's expected fastq_id path via a
# symlink (matches FILE_RENAMING_PLAN.txt's own "Alternative: Symbolic Links"
# section, and avoids duplicating ~64GB).
#
# Because the raw per-lane (L007/L008) FASTQs for RR0104/RR0508/RR0912 no
# longer exist on this filesystem either (only their already-combined
# S#_L001 pool-level outputs remain, dated 2024), the "combine" step for
# those three pools is a verification, not a regeneration: it checks the
# expected combined file already exists and is nonempty. If reproducing this
# pipeline elsewhere with the true raw per-lane FASTQs still present, see
# this file's git history for the original cat-based combine logic.

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
