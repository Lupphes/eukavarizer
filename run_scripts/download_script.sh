#!/bin/bash
#========================================================================================
# PBS Job: Download and Merge Illumina FASTQ files from GIAB
#========================================================================================
#PBS -N fastq_download_merge
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=800gb
#PBS -l walltime=12:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/illumina_giab_hg002/logs/fastq_download_merge.log

#========================================================================================
# Configuration
#========================================================================================
DATADIR=/storage/brno2/home/luppo
SCRATCH=$SCRATCHDIR

PROJECT_NAME="illumina_giab_hg002_full"
PROJECT_DIR="$DATADIR/$PROJECT_NAME"
LOGDIR="$PROJECT_DIR/logs"
OUTDIR="$PROJECT_DIR/merged_fastq"

mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

LOGFILE="$LOGDIR/fastq_download_merge.log"

FTP_BASE="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads"

# Set to "dry" for first 6 chunks or "full" for all 17 chunks
MODE="${MODE:-dry}"

cleanup() {
    echo "Cleaning up scratch..." | tee -a "$LOGFILE"
    clean_scratch
}
trap cleanup EXIT

echo "=== FASTQ Download and Merge Pipeline ===" | tee -a "$LOGFILE"
echo "Started on $(hostname) at $(date)" | tee -a "$LOGFILE"
echo "Working in scratch: $SCRATCH" | tee -a "$LOGFILE"
echo "Mode: $MODE" | tee -a "$LOGFILE"

cd "$SCRATCH"

if [ "$MODE" = "dry" ]; then
    CHUNKS="001 002 003 004 005 006"
    echo "=== Downloading DRY RUN (first 6 chunks per lane) ===" | tee -a "$LOGFILE"
else
    CHUNKS="001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017"
    echo "=== Downloading FULL RUN (all 17 chunks per lane) ===" | tee -a "$LOGFILE"
fi

#========================================================================================
# Download Files
#========================================================================================
echo "" | tee -a "$LOGFILE"
echo "Step 1: Downloading L001 R1 files..." | tee -a "$LOGFILE"
for chunk in $CHUNKS; do
    FILE="D1_S1_L001_R1_${chunk}.fastq.gz"
    echo "Downloading $FILE..." | tee -a "$LOGFILE"
    wget -q --show-progress "${FTP_BASE}/${FILE}" || curl -O "${FTP_BASE}/${FILE}"
done

echo "" | tee -a "$LOGFILE"
echo "Step 2: Downloading L001 R2 files..." | tee -a "$LOGFILE"
for chunk in $CHUNKS; do
    FILE="D1_S1_L001_R2_${chunk}.fastq.gz"
    echo "Downloading $FILE..." | tee -a "$LOGFILE"
    wget -q --show-progress "${FTP_BASE}/${FILE}" || curl -O "${FTP_BASE}/${FILE}"
done

echo "" | tee -a "$LOGFILE"
echo "Step 3: Downloading L002 R1 files..." | tee -a "$LOGFILE"
for chunk in $CHUNKS; do
    FILE="D1_S1_L002_R1_${chunk}.fastq.gz"
    echo "Downloading $FILE..." | tee -a "$LOGFILE"
    wget -q --show-progress "${FTP_BASE}/${FILE}" || curl -O "${FTP_BASE}/${FILE}"
done

echo "" | tee -a "$LOGFILE"
echo "Step 4: Downloading L002 R2 files..." | tee -a "$LOGFILE"
for chunk in $CHUNKS; do
    FILE="D1_S1_L002_R2_${chunk}.fastq.gz"
    echo "Downloading $FILE..." | tee -a "$LOGFILE"
    wget -q --show-progress "${FTP_BASE}/${FILE}" || curl -O "${FTP_BASE}/${FILE}"
done
echo "" | tee -a "$LOGFILE"
echo "=========================================" | tee -a "$LOGFILE"
echo "Download complete!" | tee -a "$LOGFILE"
echo "=========================================" | tee -a "$LOGFILE"

#========================================================================================
# Merge Files
#========================================================================================
echo "" | tee -a "$LOGFILE"
echo "Step 5: Merging L001 R1 files..." | tee -a "$LOGFILE"
cat D1_S1_L001_R1_0*.fastq.gz > "D1_S1_L001_R1_${MODE}.fastq.gz"
echo "Created D1_S1_L001_R1_${MODE}.fastq.gz" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "Step 6: Merging L001 R2 files..." | tee -a "$LOGFILE"
cat D1_S1_L001_R2_0*.fastq.gz > "D1_S1_L001_R2_${MODE}.fastq.gz"
echo "Created D1_S1_L001_R2_${MODE}.fastq.gz" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "Step 7: Merging L002 R1 files..." | tee -a "$LOGFILE"
cat D1_S1_L002_R1_0*.fastq.gz > "D1_S1_L002_R1_${MODE}.fastq.gz"
echo "Created D1_S1_L002_R1_${MODE}.fastq.gz" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "Step 8: Merging L002 R2 files..." | tee -a "$LOGFILE"
cat D1_S1_L002_R2_0*.fastq.gz > "D1_S1_L002_R2_${MODE}.fastq.gz"
echo "Created D1_S1_L002_R2_${MODE}.fastq.gz" | tee -a "$LOGFILE"

#========================================================================================
# Copy Results Back to DATADIR
#========================================================================================
echo "" | tee -a "$LOGFILE"
echo "Step 9: Copying merged files back to DATADIR..." | tee -a "$LOGFILE"
cp D1_S1_L00*_R*_${MODE}.fastq.gz "$OUTDIR/"
echo "Copied merged files to $OUTDIR/" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "=========================================" | tee -a "$LOGFILE"
echo "All done!" | tee -a "$LOGFILE"
echo "=========================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"


echo "Merged files in OUTDIR:" | tee -a "$LOGFILE"
ls -lh "$OUTDIR/D1_S1_L00*_R*_${MODE}.fastq.gz" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "=== Pipeline completed at $(date) ===" | tee -a "$LOGFILE"