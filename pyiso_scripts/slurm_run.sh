#!/bin/bash
#SBATCH --job-name=pyiso_workflow
#SBATCH --account=naiss2024-22-1378
#SBATCH --array=1-100
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=32
#SBATCH --output=/proj/proteoforma_nsc/smooth_q_to_pep/pyiso_test/logs/pyiso_%A_%a.out
#SBATCH --error=/proj/proteoforma_nsc/smooth_q_to_pep/pyiso_test/logs/pyiso_%A_%a.err

echo "[INFO]START: $(date)"

module load buildtool-easybuild/4.9.4-hpc71cbb0050 GCCcore/12.3.0
module load Python/3.11.3
source /proj/proteoforma_nsc/11analysis/envs/annotate_peptide_env/bin/activate

# Set the working directory to the scripts folder.
cd /proj/proteoforma_nsc/smooth_q_to_pep/pyiso_scripts || exit 1

# Define the list of dataset codes.
DATASETS=(PXD003868 PXD004325 PXD004424 PXD004467 PXD004536 PXD004565 PXD004947 PXD004948 PXD005025 PXD013274)

# Determine the dataset and seed based on the job array index.
TASK_ID=$SLURM_ARRAY_TASK_ID
DATASET_INDEX=$(( (TASK_ID - 1) / 10 ))
SEED=$(( (TASK_ID - 1) % 10 + 1 ))
DATASET=${DATASETS[$DATASET_INDEX]}

# Parameters for the workflow.
KNOCKING_OUT_NUM_LIST="1,2,3,4"
OUT_DIR="pyiso_test"

echo "[INFO] Running dataset $DATASET with seed $SEED (TASK_ID=$TASK_ID)"

# Call the main workflow script.
# (Adjust the path to main.sh if necessary.)
bash /proj/proteoforma_nsc/smooth_q_to_pep/pyiso_scripts/main.sh "$KNOCKING_OUT_NUM_LIST" "$OUT_DIR" "$DATASET" "$SEED"

echo "[INFO]END: $(date)"