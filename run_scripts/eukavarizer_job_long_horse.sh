#!/bin/bash
#PBS -N eukavarizer_job_long_horse
#PBS -l select=1:ncpus=16:mem=256gb:scratch_local=2000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_job_long_horse.log

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/long_job_horse/logs/eukavarizer_job_long_horse_sad.log"
mkdir -p "$(dirname "$LOGFILE")"

echo "=== Job EUKAVARIZER_JOB_SHORT started on $(hostname) at $(date) ===" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

module add openjdk/17
module add mambaforge

echo ">>> Installing newer micromamba in scratch..." | tee -a "$LOGFILE"
echo ">>> Move to scratch..." | tee -a "$LOGFILE"
# Move to scratch space
cd "$SCRATCH"

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba 2>&1 | tee -a "$LOGFILE"

mkdir -p "$SCRATCH/.mamba" "$SCRATCH/.mamba/pkgs"
export MAMBA_ROOT_PREFIX="$SCRATCH/.mamba"
export MAMBA_PKGS_DIRS="$SCRATCH/.mamba/pkgs"
export PATH="$SCRATCH/bin:$PATH"

ln -s "$SCRATCH/bin/micromamba" "$SCRATCH/bin/mamba"
ln -s "$SCRATCH/bin/micromamba" "$SCRATCH/bin/conda"

eval "$(./bin/micromamba shell hook -s bash)" 2>&1 | tee -a "$LOGFILE"
echo ">>> Micromamba version: $(micromamba --version)" | tee -a "$LOGFILE"
echo ">>> Mamba version: $(mamba --version)" | tee -a "$LOGFILE"
echo ">>> Which mamba: $(which mamba)" | tee -a "$LOGFILE"

# Clone the eukavarizer repo
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/eukavarizer.git | tee -a "$LOGFILE"

# Download Nextflow
echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

# Prepare Nextflow and Conda cache directories in scratch
mkdir -p "$SCRATCH/.conda_next" "$SCRATCH/.nextflow"
export CONDA_PKGS_DIRS="$SCRATCH/.mamba/pkgs"
export NXF_CONDA_CACHEDIR="$SCRATCH/.conda_next"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK="$DATADIR/long_job_horse/work"
export NXF_LOG_FILE="$DATADIR/long_job_horse/logs/.nextflow_short.log"
export NXF_HOME="$SCRATCH/.nextflow"
export MAMBA_ALWAYS_YES=true
export MAMBA_NO_BANNER=1

cd eukavarizer
chmod +x bin/svaba_annotate.py
chmod +x bin/simple-event-annotation.R

# Configure mamba channels
# echo ">>> Configuring mamba channels..." | tee -a "$LOGFILE"
# conda config --add channels luppo
# conda config --add channels bioconda
# conda config --add channels conda-forge

sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/assets/samplesheets/samplesheet_long_horse.csv" > "$SCRATCH/samplesheet_formatted.csv"

echo ">>> Running main Nextflow pipeline" | tee -a "$LOGFILE"
../nextflow run main.nf -profile mamba,horse,qc_off \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/long_job_horse/out" | tee -a "$LOGFILE"


echo "Cleaning up scratch..." | tee -a "$LOGFILE"
clean_scratch

echo "=== Job $PBS_JOBID completed at $(date) ===" | tee -a "$LOGFILE"
