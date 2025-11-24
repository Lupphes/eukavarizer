#!/bin/bash
#========================================================================================
# PBS Job: SV Benchmark - CMRG GIAB
#========================================================================================
#PBS -N sv_benchmark_cmrg_giab
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/sv_benchmark_cmrg_giab.log

#========================================================================================
# Configuration
#========================================================================================
# SV Benchmark using CMRG and GIAB datasets
# Platforms: ONT (~40X) and PacBio Vega HiFi (~30X)
# SV Callers: Sniffles, CuteSV, Severus, Delly, Sawfish, Dysgu
# Requirements: 200GB storage, 32GB RAM, 4 cores

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/sv_benchmark_cmrg_giab/logs/sv_benchmark.log"
mkdir -p "$(dirname "$LOGFILE")"

# Trap to ensure cleanup always happens
cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== SV Benchmark: CMRG GIAB ===" | tee -a "$LOGFILE"
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
echo ">>> Cloning SV_Benchmark_CMRG_GIAB repository to data directory..." | tee -a "$LOGFILE"
if [ ! -d "$DATADIR/SV_Benchmark_CMRG_GIAB" ]; then
    cd "$DATADIR"
    git clone https://github.com/kcleal/SV_Benchmark_CMRG_GIAB.git | tee -a "$LOGFILE"
    cd SV_Benchmark_CMRG_GIAB
    echo ">>> Fetching benchmark data to data directory..." | tee -a "$LOGFILE"
    bash fetch_data.sh | tee -a "$LOGFILE"
else
    echo ">>> Repository already exists at $DATADIR/SV_Benchmark_CMRG_GIAB" | tee -a "$LOGFILE"
fi

echo ">>> Copying repository to scratch..." | tee -a "$LOGFILE"
cd "$SCRATCH"
cp -r "$DATADIR/SV_Benchmark_CMRG_GIAB" . | tee -a "$LOGFILE"

echo ">>> Downloading Nextflow..." | tee -a "$LOGFILE"
curl -s https://get.nextflow.io | bash | tee -a "$LOGFILE"

#========================================================================================
# Nextflow Configuration
#========================================================================================
echo ">>> Configuring Nextflow..." | tee -a "$LOGFILE"
mkdir -p "$SCRATCH/.nextflow"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK="$SCRATCH/work"
export NXF_LOG_FILE="$DATADIR/sv_benchmark_cmrg_giab/logs/.nextflow.log"
export NXF_HOME="$SCRATCH/.nextflow"

cd SV_Benchmark_CMRG_GIAB

#========================================================================================
# Pipeline Execution
#========================================================================================
echo ">>> Running SV Benchmark pipeline..." | tee -a "$LOGFILE"
"$SCRATCH/nextflow" pipeline.nf | tee -a "$LOGFILE"

#========================================================================================
# Copy Results
#========================================================================================
echo ">>> Copying results to storage..." | tee -a "$LOGFILE"
mkdir -p "$DATADIR/sv_benchmark_cmrg_giab/results"
cp -r results/* "$DATADIR/sv_benchmark_cmrg_giab/results/" | tee -a "$LOGFILE"

echo "=== Benchmark completed at $(date) ===" | tee -a "$LOGFILE"
