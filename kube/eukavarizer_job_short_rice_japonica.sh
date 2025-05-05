#!/bin/bash
#PBS -N eukavarizer_job_short_rice_japonica
#PBS -l select=1:ncpus=16:mem=256gb:scratch_local=2000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_job_short_rice_japonica.log

# Define paths
DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/short_job_rice_japonica/logs/eukavarizer_job_short_rice_japonica_sad.log"
mkdir -p "$(dirname "$LOGFILE")"

echo "=== Job EUKAVARIZER_JOB_SHORT started on $(hostname) at $(date) ===" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

# Load required modules
module add openjdk/17
module add mambaforge

# Configure mamba channels
# echo ">>> Configuring mamba channels..." | tee -a "$LOGFILE"
# conda config --add channels luppo
# conda config --add channels bioconda
# conda config --add channels conda-forge

echo ">>> Move to scratch..." | tee -a "$LOGFILE"
# Move to scratch space
cd "$SCRATCH"

# Clone the eukavarizer repo
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/eukavarizer.git | tee -a "$LOGFILE"

# Download Nextflow
echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

# Prepare Conda envs in scratch
mkdir -p "$SCRATCH/.conda_pkgs" "$SCRATCH/.conda_envs" "$SCRATCH/.conda_next" "$SCRATCH/.nextflow"
export CONDA_PKGS_DIRS="$SCRATCH/.conda_pkgs"
export NXF_CONDA_CACHEDIR="$SCRATCH/.conda_next"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK=$DATADIR/short_job_rice_japonica/work
export NXF_LOG_FILE=$DATADIR/short_job_rice_japonica/logs/.nextflow_short.log
export NXF_HOME="$SCRATCH/.nextflow"

# Enter pipeline directory
cd eukavarizer
chmod +x bin/svaba_annotate.py
chmod +x bin/simple-event-annotation.R

sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/conf/samplesheets/samplesheet_rice_japonica.csv" > "$SCRATCH/samplesheet_formatted.csv"

# Actual pipeline run with inputs
echo ">>> Running main Nextflow pipeline" | tee -a "$LOGFILE"
../nextflow run main.nf -profile mamba,rice_japonica,qc_off \
    --reference_genome "$DATADIR/eukavarizer/data/39947/ref/GCF_034140825.1_ASM3414082v1_genomic.fna.gz" \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/short_job_rice_japonica/out" | tee -a "$LOGFILE"

# Clean scratch
echo "Cleaning up scratch..." | tee -a "$LOGFILE"
clean_scratch

echo "=== Job $PBS_JOBID completed at $(date) ===" | tee -a "$LOGFILE"

# GCA_001623345.3_ZS97RS3_genomic.fna.gz # INDICA
# GCF_034140825.1_ASM3414082v1_genomic.fna.gz # JAPONICA
