#!/bin/bash
#========================================================================================
# PBS Job: Eukavarizer Pipeline - Equus caballus (Horse, Taxonomy 9796)
#========================================================================================
#PBS -N eukavarizer_horse_9796_long
#PBS -l select=1:ncpus=16:mem=256gb:scratch_local=2000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_horse_9796_long.log

#========================================================================================
# Configuration
#========================================================================================
# Sample: Equus caballus (Horse)
# Taxonomy ID: 9796
# Profile: horse (from conf/samples/horse.config)
# Run type: long reads

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/horse_9796_long/logs/eukavarizer.log"
mkdir -p "$(dirname "$LOGFILE")"

# Trap to ensure cleanup always happens
cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== Eukavarizer Pipeline: Horse (9796) - Long Reads ===" | tee -a "$LOGFILE"
echo "Started on $(hostname) at $(date)" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

#========================================================================================
# Environment Setup
#========================================================================================
echo ">>> Loading modules..." | tee -a "$LOGFILE"
module add openjdk/17

echo ">>> Move to scratch..." | tee -a "$LOGFILE"
cd "$SCRATCH"

#========================================================================================
# Repository and Tool Setup
#========================================================================================
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/eukavarizer.git | tee -a "$LOGFILE"

echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

#========================================================================================
# Conda and Nextflow Configuration
#========================================================================================
echo ">>> Configuring Nextflow and Conda..." | tee -a "$LOGFILE"
mkdir -p "$SCRATCH/.conda_pkgs" "$SCRATCH/.conda_next" "$SCRATCH/.nextflow"
export CONDA_PKGS_DIRS="$SCRATCH/.conda_pkgs"
export NXF_CONDA_CACHEDIR="$SCRATCH/.conda_next"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK="$DATADIR/horse_9796_long/work"
export NXF_LOG_FILE="$DATADIR/horse_9796_long/logs/.nextflow.log"
export NXF_HOME="$SCRATCH/.nextflow"
export MAMBA_ALWAYS_YES=true
export MAMBA_NO_BANNER=1

#========================================================================================
# Pipeline Preparation
#========================================================================================
cd eukavarizer
chmod +x bin/svaba_annotate.py
chmod +x bin/simple-event-annotation.R

echo ">>> Preparing samplesheet..." | tee -a "$LOGFILE"
sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/assets/samplesheets/samplesheet_long_horse.csv" > "$SCRATCH/samplesheet_formatted.csv"

#========================================================================================
# Pipeline Execution
#========================================================================================
echo ">>> Running Eukavarizer pipeline..." | tee -a "$LOGFILE"
"$SCRATCH/nextflow" run main.nf -profile mamba,horse,qc_off \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/horse_9796_long/out" | tee -a "$LOGFILE"

echo "=== Pipeline completed at $(date) ===" | tee -a "$LOGFILE"
