#!/bin/bash
#PBS -N eukavarizer_job
#PBS -l select=1:ncpus=32:mem=512gb:scratch_local=400gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_job.log

# Exit on errors, undefined vars, and failed pipes
set -euo pipefail

# Define paths
DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/logs/eukavarizer_job_sad.log"

echo "=== Job $PBS_JOBID started on $(hostname) at $(date) ===" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

# Load required modules
module add openjdk/17
module add mambaforge

# Move to scratch space
cd "$SCRATCH"

# Clone the eukavarizer repo
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/eukavarizer.git | tee -a "$LOGFILE"

# Download Nextflow
echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

# Prepare Conda envs in scratch
mkdir -p ./.conda_pkgs ./.conda_envs
export CONDA_PKGS_DIRS=$SCRATCH/.conda_pkgs
export NXF_CONDA_CACHEDIR="$(pwd)/.conda_next"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK=/storage/brno2/home/luppo/work
export NXF_LOG_FILE=/storage/brno2/home/luppo/logs/.nextflow.log

# Enter pipeline directory
cd eukavarizer

# First dry-run to build cache
echo ">>> Running first Nextflow command (cache warm-up)" | tee -a "$LOGFILE"
../nextflow run main.nf -profile mamba,mix_medium,qc_off -resume | tee -a "$LOGFILE"

# Actual pipeline run with inputs
echo ">>> Running main Nextflow pipeline" | tee -a "$LOGFILE"
../nextflow run main.nf -profile mamba,mix_medium,qc_off \
  --taxonomy_id 9606 \
  --reference_genome "$DATADIR/data/9606/ref/hg38.fa.gz" \
  --sequence_dir "$DATADIR/eukavarizer/src/out/" \
  --outdir "$DATADIR/out_big" | tee -a "$LOGFILE"

echo ">>> Cleaning up any broken conda environments..." | tee -a "$LOGFILE"
find "$NXF_CONDA_CACHEDIR" -type d -name "envs" -exec rm -rf {} + || true

# Clean scratch
echo "Cleaning up scratch..." | tee -a "$LOGFILE"
clean_scratch

echo "=== Job $PBS_JOBID completed at $(date) ===" | tee -a "$LOGFILE"
