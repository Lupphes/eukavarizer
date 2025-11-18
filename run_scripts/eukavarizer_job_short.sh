#!/bin/bash
#========================================================================================
# PBS Job: Eukavarizer Pipeline - Homo sapiens (Human, Taxonomy 9606) - Short Reads
#========================================================================================
#PBS -N eukavarizer_human_9606_short
#PBS -l select=1:ncpus=32:mem=728gb:scratch_local=2000gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_human_9606_short.log

#========================================================================================
# Configuration
#========================================================================================
# Sample: Homo sapiens (Human)
# Taxonomy ID: 9606
# Reference: T2T-CHM13v2.0 (GCF_009914755.1)
# Profile: short_full (Delly, Manta, SVABA, GRIDSS, Tiddit, Dysgu)
# Run type: Short reads (Illumina)

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/human_9606_short/logs/eukavarizer.log"
mkdir -p "$(dirname "$LOGFILE")"
export COLUMNS=200 

# Trap to ensure cleanup always happens
cleanup() {
    echo "Not cleaning up scratch..." | tee -a "$LOGFILE"
    # clean_scratch
}
trap cleanup EXIT

echo "=== Eukavarizer Pipeline: Human (9606) - Short Reads ===" | tee -a "$LOGFILE"
echo "Started on $(hostname) at $(date)" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"

#========================================================================================
# Environment Setup
#========================================================================================
echo ">>> Loading modules..." | tee -a "$LOGFILE"
module add openjdk/17
module add micromamba/2.3.3

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
mkdir -p "$SCRATCH/.conda_pkgs" "$SCRATCH/.conda_next" "$SCRATCH/.nextflow" "$SCRATCH/.micromamba"
export CONDA_PKGS_DIRS="$SCRATCH/.conda_pkgs"
export NXF_CONDA_CACHEDIR="$SCRATCH/.conda_next"
export MAMBA_EXE=micromamba
export MAMBA_ROOT_PREFIX="$SCRATCH/.micromamba"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK="$SCRATCH/work"
export NXF_LOG_FILE="$DATADIR/human_9606_short/logs/.nextflow.log"
export NXF_HOME="$SCRATCH/.nextflow"
export MAMBA_ALWAYS_YES=true

#========================================================================================
# Pipeline Preparation
#========================================================================================
cd eukavarizer
chmod +x bin/svaba_annotate.py
chmod +x bin/simple-event-annotation.R

echo ">>> Preparing samplesheet..." | tee -a "$LOGFILE"
sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/assets/samplesheets/samplesheet_human_short.csv" > "$SCRATCH/samplesheet_formatted.csv"

echo ">>> Staging reference genome to scratch..." | tee -a "$LOGFILE"
mkdir -p "$SCRATCH/reference"
cp "$DATADIR/eukavarizer/data/9606/ref/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" "$SCRATCH/reference/"
echo "Reference ready at $SCRATCH/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" | tee -a "$LOGFILE"

#========================================================================================
# Pipeline Execution
#========================================================================================
echo ">>> Running Eukavarizer pipeline..." | tee -a "$LOGFILE"
"$SCRATCH/nextflow" run main.nf -profile mamba,short_full,qc_off \
    --taxonomy_id 9606 \
    --reference_genome "$SCRATCH/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/human_9606_short/out" \
    --seqtk_size 1.0 --seqtk_flag false --svaba_flag true --bwamem2 false --minimap2_flag false | tee -a "$LOGFILE"

echo "=== Pipeline completed at $(date) ===" | tee -a "$LOGFILE"
