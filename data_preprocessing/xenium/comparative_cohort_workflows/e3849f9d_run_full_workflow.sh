#!/usr/bin/env bash
set -euo pipefail

ROOT="${METABOLIC_COMPARATIVE_ROOT:-$(pwd)/comparative_seven_donor}"
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${ROOT}/inputs"
MODEL_DIR="${ROOT}/model_artifacts/resolvi_run"
OUTPUT_DIR="${ROOT}/outputs"

PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"
SEED="${SEED:-20260729}"
MAX_EPOCHS="${MAX_EPOCHS:-100}"
BATCH_SIZE="${BATCH_SIZE:-512}"

SFE="${INPUT_DIR}/sfe_seven_donors_annotated_reference_mapped_20260729.rds"
H5AD_INPUT="${INPUT_DIR}/seven_donor_resolvi_input_20260729.h5ad"
IMC_SPE="${INPUT_DIR}/spe2025-12-17 umap.RData"
MATRISOME="${INPUT_DIR}/Hs_Matrisome_Masterlist_Naba.csv"
TRAINING_LOG="${MODEL_DIR}/training.log"

mkdir -p "${MODEL_DIR}" "${OUTPUT_DIR}"

for required in "${SFE}" "${IMC_SPE}" "${MATRISOME}"; do
  if [[ ! -f "${required}" ]]; then
    printf 'Missing required input: %s\n' "${required}" >&2
    exit 1
  fi
done

if [[ ! -f "${H5AD_INPUT}" ]]; then
  "${R_BIN}" "${WORKFLOW_DIR}/dafeb4cf_export_seven_donor_sfe_to_h5ad.R" \
    "${SFE}" \
    "${H5AD_INPUT}"
fi

train_command=(
  "${PYTHON_BIN}"
  "${WORKFLOW_DIR}/a1989afe_train_seven_donor_resolvi.py"
  --input-h5ad "${H5AD_INPUT}"
  --output-dir "${MODEL_DIR}"
  --n-latent 30
  --max-epochs "${MAX_EPOCHS}"
  --batch-size "${BATCH_SIZE}"
  --seed "${SEED}"
  --accelerator cpu
  --device 1
  --num-workers 0
  --include-corrected-px-rate
  --posterior-samples 3
)

if [[ -s "${MODEL_DIR}/run_manifest.json" &&
      -s "${MODEL_DIR}/seven_donor_resolvi.h5ad" &&
      -s "${TRAINING_LOG}" ]]; then
  "${train_command[@]}"
else
  "${train_command[@]}" 2>&1 | tee "${TRAINING_LOG}"
fi

"${PYTHON_BIN}" "${WORKFLOW_DIR}/c249ddcb_validate_seven_donor_resolvi.py" \
  --input-h5ad "${H5AD_INPUT}" \
  --run-dir "${MODEL_DIR}" \
  --training-log "${TRAINING_LOG}" \
  --output-dir "${MODEL_DIR}/validation"

"${PYTHON_BIN}" "${WORKFLOW_DIR}/59ec9a4c_annotate_joint_resolvi.py" \
  --model-h5ad "${MODEL_DIR}/seven_donor_resolvi.h5ad" \
  --output-dir "${OUTPUT_DIR}" \
  --seed "${SEED}"

cp "${MODEL_DIR}/run_manifest.json" \
  "${OUTPUT_DIR}/resolvi_training_manifest.json"

"${R_BIN}" "${WORKFLOW_DIR}/7b088d88_03_attach_joint_resolvi_and_marker_qc.R" \
  "${SFE}" \
  "${OUTPUT_DIR}" \
  "${OUTPUT_DIR}"

"${R_BIN}" "${WORKFLOW_DIR}/096adafb_04_compute_seven_donor_scores.R" \
  "${MATRISOME}" \
  "${ROOT}"

"${R_BIN}" "${WORKFLOW_DIR}/b24c6f0c_05_generate_seven_donor_manuscript_panels.R" \
  "${IMC_SPE}" \
  "${ROOT}"


test -s "${OUTPUT_DIR}/sfe_seven_donors_annotated_final.rds"
"${R_BIN}" "${WORKFLOW_DIR}/210f098d_validate_report_outputs.R" "${ROOT}"
printf 'Workflow complete: %s\n' "${ROOT}"
