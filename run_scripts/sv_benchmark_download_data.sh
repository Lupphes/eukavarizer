#!/bin/bash
#PBS -N sv_benchmark_download
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=50gb
#PBS -l walltime=12:00:00
#PBS -m abe
#PBS -M ondrej.sloup@protonmail.com
#PBS -j oe
#PBS -o /storage/brno2/home/luppo/logs/sv_benchmark_download.log

DATADIR=/storage/brno2/home/luppo
LOGFILE="$DATADIR/logs/sv_benchmark_download.log"
mkdir -p "$(dirname "$LOGFILE")"

echo "=== SV Benchmark Data Download (~200GB) ===" | tee -a "$LOGFILE"
echo "Started on $(hostname) at $(date)" | tee -a "$LOGFILE"

if [ -d "$DATADIR/SV_Benchmark_CMRG_GIAB" ]; then
    echo "Repository exists, pulling latest changes..." | tee -a "$LOGFILE"
    cd "$DATADIR/SV_Benchmark_CMRG_GIAB"
    git pull | tee -a "$LOGFILE"
else
    cd "$DATADIR"
    git clone https://github.com/Lupphes/SV_Benchmark_CMRG_GIAB.git | tee -a "$LOGFILE"
    cd SV_Benchmark_CMRG_GIAB
fi

echo ">>> Downloading benchmark data..." | tee -a "$LOGFILE"
bash fetch_data.sh | tee -a "$LOGFILE"

echo ">>> Downloading GIAB HG002 PacBio BAM..." | tee -a "$LOGFILE"
mkdir -p "$DATADIR/SV_Benchmark_CMRG_GIAB/giab_data"
cd "$DATADIR/SV_Benchmark_CMRG_GIAB/giab_data"

FTP_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/alignment"
echo "Downloading BAM file (63GB)..." | tee -a "$LOGFILE"
wget -c "${FTP_BASE}/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam" | tee -a "$LOGFILE"

echo "Downloading BAM index..." | tee -a "$LOGFILE"
wget -c "${FTP_BASE}/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam.bai" | tee -a "$LOGFILE"

echo "=== Completed at $(date) ===" | tee -a "$LOGFILE"
echo "Data location: $DATADIR/SV_Benchmark_CMRG_GIAB" | tee -a "$LOGFILE"
echo "BAM file: $DATADIR/SV_Benchmark_CMRG_GIAB/giab_data/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam" | tee -a "$LOGFILE"
