#!/bin/bash
#PBS -N sv_benchmark_cmrg_giab
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=1000gb
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/sv_benchmark_cmrg_giab.log

DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR
LOGFILE="$DATADIR/sv_benchmark_cmrg_giab/logs/sv_benchmark.log"
mkdir -p "$(dirname "$LOGFILE")"

# Pipeline Parameters
DATASET="giab"  # Options: "giab", "cmrg", "both"
BENCH_PARAMS="--passonly -r 1000 --dup-to-ins -p 0 --sizemax 50000"
TEST_MODE="false"
OUTDIR="$DATADIR/sv_benchmark_cmrg_giab/results"
REF_GENOME="$DATADIR/eukavarizer/data/9606/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ONT_DATA=""
# GIAB HG002 PacBio BAM (aligned to GRCh38) - downloaded by sv_benchmark_download_data.sh
PACBIO_DATA="$DATADIR/SV_Benchmark_CMRG_GIAB/giab_data/HG002.15kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam"

cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== SV Benchmark: CMRG GIAB ===" | tee -a "$LOGFILE"
echo "Started on $(hostname) at $(date)" | tee -a "$LOGFILE"

module add micromamba/2.3.3
module add openjdk/17
cd "$SCRATCH"

# Check data exists
if [ ! -d "$DATADIR/SV_Benchmark_CMRG_GIAB" ]; then
    echo "ERROR: Repository not found. Run sv_benchmark_download_data.sh first" | tee -a "$LOGFILE"
    exit 1
fi

#========================================================================================
# Repository and Tool Setup
#========================================================================================
echo ">>> Cloning repository..." | tee -a "$LOGFILE"
git clone https://github.com/Lupphes/SV_Benchmark_CMRG_GIAB.git | tee -a "$LOGFILE"

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
export NXF_LOG_FILE="$DATADIR/sv_benchmark_cmrg_giab/logs/.nextflow.log"
export NXF_HOME="$SCRATCH/.nextflow"
export MAMBA_ALWAYS_YES=true

cd SV_Benchmark_CMRG_GIAB

#========================================================================================
# Copy/Link Input Files to Current Directory
#========================================================================================
echo ">>> Copying input files to scratch..." | tee -a "$LOGFILE"

# Decompress and copy reference genome
if [ -n "$REF_GENOME" ]; then
    echo "  Copying reference genome..." | tee -a "$LOGFILE"
    if [[ "$REF_GENOME" == *.gz ]]; then
        gunzip -c "$REF_GENOME" > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        REF_GENOME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    else
        cp "$REF_GENOME" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        REF_GENOME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    fi
fi

# Copy PacBio data if provided
if [ -n "$PACBIO_DATA" ]; then
    echo "  Copying PacBio data..." | tee -a "$LOGFILE"
    PACBIO_BASENAME=$(basename "$PACBIO_DATA")
    cp "$PACBIO_DATA" "$PACBIO_BASENAME"
    # Copy index files (try both .bai and .bam.bai patterns)
    [ -f "${PACBIO_DATA}.bai" ] && cp "${PACBIO_DATA}.bai" "${PACBIO_BASENAME}.bai"
    [ -f "${PACBIO_DATA%.bam}.bai" ] && cp "${PACBIO_DATA%.bam}.bai" "${PACBIO_BASENAME%.bam}.bai"
    PACBIO_DATA="$PACBIO_BASENAME"
fi
# Copy ONT data if provided
if [ -n "$ONT_DATA" ]; then
    echo "  Copying ONT data..." | tee -a "$LOGFILE"
    ONT_BASENAME=$(basename "$ONT_DATA")
    cp "$ONT_DATA" "$ONT_BASENAME"
    # Copy index files (try .crai for CRAM or .bai for BAM)
    [ -f "${ONT_DATA}.crai" ] && cp "${ONT_DATA}.crai" "${ONT_BASENAME}.crai"
    [ -f "${ONT_DATA}.bai" ] && cp "${ONT_DATA}.bai" "${ONT_BASENAME}.bai"
    ONT_DATA="$ONT_BASENAME"
fi

echo ">>> Running pipeline with parameters:" | tee -a "$LOGFILE"
echo "  Dataset: $DATASET" | tee -a "$LOGFILE"
echo "  Bench params: $BENCH_PARAMS" | tee -a "$LOGFILE"
echo "  Test mode: $TEST_MODE" | tee -a "$LOGFILE"
[ -n "$REF_GENOME" ] && echo "  Reference: $REF_GENOME" | tee -a "$LOGFILE"
[ -n "$ONT_DATA" ] && echo "  ONT data: $ONT_DATA" | tee -a "$LOGFILE"
[ -n "$PACBIO_DATA" ] && echo "  PacBio data: $PACBIO_DATA" | tee -a "$LOGFILE"

NXF_CMD="$SCRATCH/nextflow pipeline.nf --bench_params \"$BENCH_PARAMS\" --test $TEST_MODE --dataset $DATASET"
[ -n "$REF_GENOME" ] && NXF_CMD="$NXF_CMD --ref $REF_GENOME"
[ -n "$ONT_DATA" ] && NXF_CMD="$NXF_CMD --ont_data $ONT_DATA"
[ -n "$PACBIO_DATA" ] && NXF_CMD="$NXF_CMD --pacbio_data $PACBIO_DATA"

eval $NXF_CMD | tee -a "$LOGFILE"

mkdir -p "$OUTDIR"
cp -r results/* "$OUTDIR/" | tee -a "$LOGFILE"
echo "Results saved to: $OUTDIR" | tee -a "$LOGFILE"
echo "=== Completed at $(date) ===" | tee -a "$LOGFILE"
# Copy ONT data if provided
if [ -n "$ONT_DATA" ]; then
    echo "  Copying ONT data..." | tee -a "$LOGFILE"
    ONT_BASENAME=$(basename "$ONT_DATA")
    cp "$ONT_DATA" "$ONT_BASENAME"
    # Copy index files (try .crai for CRAM or .bai for BAM)
    [ -f "${ONT_DATA}.crai" ] && cp "${ONT_DATA}.crai" "${ONT_BASENAME}.crai"
    [ -f "${ONT_DATA}.bai" ] && cp "${ONT_DATA}.bai" "${ONT_BASENAME}.bai"
    ONT_DATA="$ONT_BASENAME"
fi

echo ">>> Running pipeline with parameters:" | tee -a "$LOGFILE"
echo "  Dataset: $DATASET" | tee -a "$LOGFILE"
echo "  Bench params: $BENCH_PARAMS" | tee -a "$LOGFILE"
echo "  Test mode: $TEST_MODE" | tee -a "$LOGFILE"
[ -n "$REF_GENOME" ] && echo "  Reference: $REF_GENOME" | tee -a "$LOGFILE"
[ -n "$ONT_DATA" ] && echo "  ONT data: $ONT_DATA" | tee -a "$LOGFILE"
[ -n "$PACBIO_DATA" ] && echo "  PacBio data: $PACBIO_DATA" | tee -a "$LOGFILE"

NXF_CMD="$SCRATCH/nextflow pipeline.nf --bench_params \"$BENCH_PARAMS\" --test $TEST_MODE --dataset $DATASET"
[ -n "$REF_GENOME" ] && NXF_CMD="$NXF_CMD --ref $REF_GENOME"
[ -n "$ONT_DATA" ] && NXF_CMD="$NXF_CMD --ont_data $ONT_DATA"
[ -n "$PACBIO_DATA" ] && NXF_CMD="$NXF_CMD --pacbio_data $PACBIO_DATA"

eval $NXF_CMD | tee -a "$LOGFILE"

mkdir -p "$OUTDIR"
cp -r results/* "$OUTDIR/" | tee -a "$LOGFILE"
echo "Results saved to: $OUTDIR" | tee -a "$LOGFILE"
echo "=== Completed at $(date) ===" | tee -a "$LOGFILE"