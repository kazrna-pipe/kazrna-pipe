#!/usr/bin/env bash
# bin/fetch_tcga_esca.sh
#
# Downloads TCGA-ESCA gene expression counts (STAR, harmonised by GDC) and the
# sample metadata needed to tell tumour from adjacent-normal.
#
# These are open-access, pre-computed counts. The raw BAMs are controlled and
# require a dbGaP application (phs000178), which is why the TCGA arm of the
# analysis uses GDC's harmonised counts rather than passing through this
# pipeline. The Kazakh cohort is processed end to end by KazRNA-Pipe.
#
#   bash bin/fetch_tcga_esca.sh data/tcga_esca

set -euo pipefail

OUT="${1:-data/tcga_esca}"
API="https://api.gdc.cancer.gov"

mkdir -p "${OUT}"
cd "${OUT}"

command -v gdc-client >/dev/null || {
    echo "gdc-client not found. Install from:"
    echo "  https://gdc.cancer.gov/access-data/gdc-data-transfer-tool"
    exit 1
}

# ---- 1. Query the API directly ---------------------------------------------
# Filtering here rather than through the portal UI keeps the selection
# reproducible: this file records exactly which cohort was requested.
echo ">> Requesting the file manifest"
cat > filters.json <<'EOF'
{
  "filters": {
    "op": "and",
    "content": [
      {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-ESCA"]}},
      {"op": "in", "content": {"field": "data_category",   "value": ["Transcriptome Profiling"]}},
      {"op": "in", "content": {"field": "data_type",       "value": ["Gene Expression Quantification"]}},
      {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
      {"op": "in", "content": {"field": "access",          "value": ["open"]}}
    ]
  },
  "format": "TSV",
  "size": "2000",
  "return_type": "manifest"
}
EOF

curl -s --header "Content-Type: application/json" \
     --request POST --data @filters.json \
     "${API}/files" > gdc_manifest.tsv

n=$(( $(wc -l < gdc_manifest.tsv) - 1 ))
echo ">> Manifest lists ${n} files"
[ "${n}" -gt 100 ] || { echo "Unexpectedly few files; check filters.json"; exit 1; }

# ---- 2. Sample metadata ----------------------------------------------------
# The file UUID alone does not say whether a sample is tumour or normal. The
# sample barcode does: field 4 is 01 for primary tumour, 11 for solid tissue
# normal. Without this mapping the counts cannot be assigned to a condition.
echo ">> Requesting sample metadata"
cat > meta_filters.json <<'EOF'
{
  "filters": {
    "op": "and",
    "content": [
      {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-ESCA"]}},
      {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
      {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}
    ]
  },
  "fields": "file_id,file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.demographic.gender,cases.demographic.age_at_index,cases.diagnoses.primary_diagnosis,cases.diagnoses.tumor_stage",
  "format": "TSV",
  "size": "2000"
}
EOF

curl -s --header "Content-Type: application/json" \
     --request POST --data @meta_filters.json \
     "${API}/files" > sample_metadata.tsv

echo ">> Sample types present:"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /sample_type/) c=i; next} c{print $c}' \
    sample_metadata.tsv | sort | uniq -c

# ---- 3. Download -----------------------------------------------------------
echo ">> Downloading (about 850 MB)"
gdc-client download -m gdc_manifest.tsv -d .

echo
echo ">> Done. Files in $(pwd)"
echo "   Build the samplesheet from sample_metadata.tsv: sample type 'Primary Tumor'"
echo "   and 'Solid Tissue Normal' give the within-study contrast."
