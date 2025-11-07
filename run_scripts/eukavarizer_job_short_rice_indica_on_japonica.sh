#!/bin/bash
#========================================================================================
# PBS Job: Eukavarizer Pipeline - Oryza sativa (Rice Indica, Taxonomy 39946) on Japonica reference
#========================================================================================
#PBS -N eukavarizer_rice_39946_indica_on_japonica
#PBS -l select=1:ncpus=16:mem=256gb:scratch_local=2000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_rice_39946_indica_on_japonica.log

#========================================================================================
# Configuration
#========================================================================================
# Sample: Oryza sativa Indica Group (Rice Indica) on Japonica reference
# Taxonomy ID: 39946
# Reference: ZS97RS3 (GCA_001623345.3) - Japonica reference
# Profile: rice_indica (Delly, Manta, Dysgu, Tiddit, SVABA, Sniffles, CuteSV)
# Run type: Short reads (Illumina)

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/rice_39946_indica_on_japonica/logs/eukavarizer.log"
mkdir -p "$(dirname "$LOGFILE")"
export COLUMNS=200 

# Trap to ensure cleanup always happens
cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== Eukavarizer Pipeline: Rice Indica (39946) on Japonica reference ===" | tee -a "$LOGFILE"
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
export NXF_LOG_FILE="$DATADIR/rice_39946_indica_on_japonica/logs/.nextflow.log"
export NXF_HOME="$SCRATCH/.nextflow"
export MAMBA_ALWAYS_YES=true

#========================================================================================
# Pipeline Preparation
#========================================================================================
cd eukavarizer
chmod +x bin/svaba_annotate.py
chmod +x bin/simple-event-annotation.R

echo ">>> Preparing samplesheet..." | tee -a "$LOGFILE"
sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/assets/samplesheets/samplesheet_rice_indica.csv" > "$SCRATCH/samplesheet_formatted.csv"

echo ">>> Staging reference genome to scratch..." | tee -a "$LOGFILE"
mkdir -p "$SCRATCH/reference"
cp "$DATADIR/eukavarizer/data/39947/ref/GCF_034140825.1_ASM3414082v1_genomic.fna.gz" "$SCRATCH/reference/"
echo "Reference ready at $SCRATCH/reference/GCF_034140825.1_ASM3414082v1_genomic.fna.gz" | tee -a "$LOGFILE"

#========================================================================================
# Pipeline Execution
#========================================================================================

# GCA_001623345.3_ZS97RS3_genomic.fna.gz # INDICA
# GCF_034140825.1_ASM3414082v1_genomic.fna.gz # JAPONICA
# "$DATADIR/eukavarizer/data/39947/ref/GCF_034140825.1_ASM3414082v1_genomic.fna.gz" \

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034140825.1/
echo ">>> Running Eukavarizer pipeline..." | tee -a "$LOGFILE"
"$SCRATCH/nextflow" run main.nf -profile mamba,rice_indica,qc_off \
    --reference_genome "$SCRATCH/reference/GCF_034140825.1_ASM3414082v1_genomic.fna.gz" \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/rice_39946_indica_on_japonica/out" | tee -a "$LOGFILE"

echo "=== Pipeline completed at $(date) ===" | tee -a "$LOGFILE"
