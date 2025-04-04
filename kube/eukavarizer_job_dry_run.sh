#!/bin/bash
#PBS -N eukavarizer_job_dry_run
#PBS -l select=1:ncpus=64:mem=512gb:scratch_local=400gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_job_dry_run.log

# Exit on errors, undefined vars, and failed pipes
set -euo pipefail

# Define paths
DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/dry_run_job/logs/eukavarizer_job_dry_run_sad.log"
mkdir -p "$(dirname "$LOGFILE")"

echo "=== Job EUKAVARIZER_JOB_DRY_RUN started on $(hostname) at $(date) ===" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

# Load required modules
module add openjdk/17
module add mambaforge

# Configure mamba channels
echo ">>> Configuring mamba channels..." | tee -a "$LOGFILE"
conda config --add channels luppo
conda config --add channels bioconda
conda config --add channels conda-forge

# Move to scratch space
cd "$SCRATCH"

# Clone the eukavarizer repo
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/eukavarizer.git | tee -a "$LOGFILE"

# Download Nextflow
echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

# Prepare Conda envs in scratch
mkdir -p ./.conda_pkgs ./.conda_envs ./.conda_dir
export CONDA_PKGS_DIRS=$SCRATCH/.conda_pkgs
export NXF_CONDA_CACHEDIR="$(pwd)/.conda_next"
export CONDA_PKGS_DIRS="$(pwd)/.conda_dir"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK=$DATADIR/dry_run_job/work
export NXF_LOG_FILE=$DATADIR/dry_run_job/logs/.nextflow_dry.log

# Enter pipeline directory
cd eukavarizer

# Actual pipeline run with inputs
echo ">>> Running main Nextflow pipeline" | tee -a "$LOGFILE"
../nextflow run main.nf -profile mamba,mix_medium,qc_off \
    --taxonomy_id 9606 \
    --reference_genome "$DATADIR/eukavarizer/data/9606/ref/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
    --sequence_dir "$DATADIR/eukavarizer/data/9606/gib" \
    --outdir "$DATADIR/dry_run_job/out" --seqtk_size 1.0 --seqtk_flag false | tee -a "$LOGFILE"

echo ">>> Cleaning up any broken conda environments..." | tee -a "$LOGFILE"
find "$NXF_CONDA_CACHEDIR" -type d -name "envs" -exec rm -rf {} + || true

# Clean scratch
echo "Cleaning up scratch..." | tee -a "$LOGFILE"
clean_scratch

echo "=== Job $PBS_JOBID completed at $(date) ===" | tee -a "$LOGFILE"
