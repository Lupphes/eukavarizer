#!/bin/bash
#PBS -N eukavarizer_job
#PBS -l select=1:ncpus=32:mem=256gb:scratch_local=400gb
#PBS -l walltime=24:00:00
#PBS -m abe                           # Send mail on abort (a), begin (b), and end (e)
#PBS -M ondrej.sloup@protonmail.com         # Replace with your email address
#PBS -j oe                           # Join stdout and stderr into one file
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_${PBS_JOBID}.log

# Define working directories
DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR

echo "=== Job $PBS_JOBID started on `hostname` at `date` ==="
echo "Working in scratch: $SCRATCH"

# Load necessary modules
module add openjdk/17
module add mambaforge

# Move to scratch
cd $SCRATCH || { echo "Failed to enter scratch dir"; exit 1; }

# Copy project files
cp -r $DATADIR/eukavarizer . || { echo "Failed to copy eukavarizer project"; exit 2; }

# Download Nextflow
curl -s https://get.nextflow.io | bash || { echo "Failed to download Nextflow"; exit 3; }

# Create local conda directories
mkdir -p ./.conda_pkgs ./.conda_envs
export CONDA_PKGS_DIRS=$SCRATCH/.conda_pkgs
export CONDA_ENVS_PATH=$SCRATCH/.conda_envs

# Enter the pipeline directory
cd eukavarizer || { echo "Failed to enter eukavarizer dir"; exit 4; }

# First Nextflow run (resume cache)
echo ">>> Running first Nextflow command"
../nextflow run main.nf -profile mamba,mix_medium,qc_off -resume || { echo "Initial NF run failed"; exit 5; }

# Second Nextflow run with actual inputs
echo ">>> Running second Nextflow command"
../nextflow run main.nf -profile mamba,mix_medium,qc_off \
  --taxonomy_id 9606 \
  --reference_genome "$DATADIR/data/9606/ref/hg38.fa.gz" \
  --sequence_dir "$DATADIR/data/9606/gib" \
  --outdir "$DATADIR/out" || { echo "Second NF run failed"; exit 6; }

# Clean up scratch
echo "Cleaning up scratch..."
clean_scratch

echo "=== Job $PBS_JOBID completed at `date` ==="
