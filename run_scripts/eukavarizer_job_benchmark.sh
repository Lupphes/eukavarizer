#!/bin/bash
#========================================================================================
# PBS Job: Eukavarizer Pipeline - Homo sapiens (Human, Taxonomy 9606) - Benchmark GRCh38
#========================================================================================
#PBS -N eukavarizer_human_9606_benchmark_grch38
#PBS -l select=1:ncpus=64:mem=768gb:scratch_local=1000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/eukavarizer_human_9606_benchmark_grch38.log

#========================================================================================
# Configuration
#========================================================================================
# Sample: Homo sapiens (Human)
# Taxonomy ID: 9606
# Reference: GRCh38 (GCA_000001405.15)
# Profile: long_full (Delly, Manta, Sniffles, CuteSV)
# Run type: Benchmark with GRCh38 reference

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/human_9606_benchmark_grch38/logs/eukavarizer.log"
mkdir -p "$(dirname "$LOGFILE")"

# Trap to ensure cleanup always happens
cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== Eukavarizer Pipeline: Human (9606) - Benchmark GRCh38 ===" | tee -a "$LOGFILE"
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
export NXF_WORK="$DATADIR/human_9606_benchmark_grch38/work"
export NXF_LOG_FILE="$DATADIR/human_9606_benchmark_grch38/logs/.nextflow.log"
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
sed "s|\$DATADIR|$DATADIR|g" "$DATADIR/eukavarizer/assets/samplesheets/samplesheet_human_pacbio.csv" > "$SCRATCH/samplesheet_formatted.csv"

#========================================================================================
# Pipeline Execution
#========================================================================================
echo ">>> Running Eukavarizer pipeline..." | tee -a "$LOGFILE"
"$SCRATCH/nextflow" run main.nf -profile mamba,long_full,qc_off \
    --taxonomy_id 9606 \
    --reference_genome "$DATADIR/eukavarizer/data/9606/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    --input "$SCRATCH/samplesheet_formatted.csv" \
    --outdir "$DATADIR/human_9606_benchmark_grch38/out" \
    --seqtk_size 1.0 --seqtk_flag false --minimap2_flag true | tee -a "$LOGFILE"

echo "=== Pipeline completed at $(date) ===" | tee -a "$LOGFILE"
