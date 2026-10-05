#!/usr/bin/env bash
set -euo pipefail

# Write Cell Ranger multi configs for four FLEX pools (BC001-BC004).
# Reference: mm10-2020-A; mouse transcriptome probe set v1.0.1.

CONFIG_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config"
mkdir -p "$CONFIG_DIR"


## RR1316
cat > "$CONFIG_DIR/RR1316_config.csv" <<'EOF'
[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR1316,/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
13,BC001,KO Male HFD
14,BC002,KO Male HFD
15,BC003,KO Female HFD
16,BC004,KO Female HFD
EOF
echo "Wrote $CONFIG_DIR/RR1316_config.csv"

## RR0104
cat > "$CONFIG_DIR/RR0104_config.csv" <<'EOF'
[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0104,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
1,BC001,WT Male HFD
2,BC002,WT Female HFD
3,BC003,KO Male HFD
4,BC004,KO Female HFD
EOF
echo "Wrote $CONFIG_DIR/RR0104_config.csv"

## RR0508
cat > "$CONFIG_DIR/RR0508_config.csv" <<'EOF'
[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0508,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
5,BC001,WT Female HFD Infected
6,BC002,WT Female HFD Infected
7,BC003,KO Female HFD Infected
8,BC004,KO Female HFD Infected
EOF
echo "Wrote $CONFIG_DIR/RR0508_config.csv"

## RR0912
cat > "$CONFIG_DIR/RR0912_config.csv" <<'EOF'
[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0912,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
9,BC001,WT Male HFD
10,BC002,WT Male HFD
11,BC003,WT Female HFD
12,BC004,WT Female HFD
EOF
echo "Wrote $CONFIG_DIR/RR0912_config.csv"

echo "Done: 4 multiconfig CSVs written to $CONFIG_DIR"
