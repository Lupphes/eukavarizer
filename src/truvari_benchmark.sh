#!/bin/bash

# Truvari Benchmarking Script
# This script runs truvari bench to evaluate SV calling performance
# Usage: ./truvari_benchmark.sh <input_vcf> <truth_vcf> <reference_fasta> <bed_file> <output_dir>
# ./src/truvari_benchmark.sh \
#   -i data_collapsed/bcftools_concat.vcf.gz \
#   -t true/HG002_SVs_Tier1_v0.6.vcf.gz \
#   -r data_collapsed/hs37d5.fa \
#   -b true/HG002_SVs_Tier1_v0.6.bed \
#   -o result_3

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" >&2
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" >&2
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" >&2
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Arguments:
    -i, --input VCF          Input VCF file to benchmark (can be .vcf or .vcf.gz)
    -t, --truth VCF          Truth/baseline VCF file (.vcf.gz)
    -r, --reference FASTA    Reference genome FASTA file
    -b, --bed FILE           BED file with regions to include
    -o, --output DIR         Output directory

Optional Arguments:
    -c, --collapse           Run truvari collapse before benchmarking
    --refdist INT            Maximum distance for variants [default: 1000]
    --pctseq FLOAT           Percent sequence similarity [default: 0.0]
    --pctsize FLOAT          Percent size similarity [default: 0.7]
    --dup-to-ins             Convert duplications to insertions
    --passonly               Only consider PASS variants
    --no-collapse-args       Use default collapse parameters (strict)
    -h, --help               Show this help message

Example:
    $0 -i calls.vcf.gz -t truth.vcf.gz -r ref.fa -b regions.bed -o results --collapse

    # With custom parameters
    $0 -i calls.vcf.gz -t truth.vcf.gz -r ref.fa -b regions.bed -o results \\
        --collapse --refdist 500 --pctsize 0.9

EOF
    exit 1
}

sort_and_index_vcf() {
    local input_vcf="$1"
    local output_dir="$2"
    local basename=$(basename "$input_vcf" .vcf.gz)
    local sorted_vcf="${output_dir}/${basename}.sorted.vcf.gz"
    
    log_info "Sorting VCF: $input_vcf"
    bcftools sort "$input_vcf" -O z -o "$sorted_vcf"
    
    log_info "Indexing sorted VCF..."
    tabix -p vcf "$sorted_vcf"
    
    echo "$sorted_vcf"
}

INPUT_VCF=""
TRUTH_VCF=""
REFERENCE=""
BED_FILE=""
OUTPUT_DIR=""
RUN_COLLAPSE=false
REFDIST=1000
PCTSEQ=0.0
PCTSIZE=0.7
DUP_TO_INS="--dup-to-ins"
PASSONLY="--passonly"
COLLAPSE_REFDIST=1000
COLLAPSE_PCTSEQ=0.7
COLLAPSE_PCTSIZE=0.7
USE_STRICT_COLLAPSE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_VCF="$2"
            shift 2
            ;;
        -t|--truth)
            TRUTH_VCF="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -b|--bed)
            BED_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--collapse)
            RUN_COLLAPSE=true
            shift
            ;;
        --refdist)
            REFDIST="$2"
            shift 2
            ;;
        --pctseq)
            PCTSEQ="$2"
            shift 2
            ;;
        --pctsize)
            PCTSIZE="$2"
            shift 2
            ;;
        --dup-to-ins)
            DUP_TO_INS="--dup-to-ins"
            shift
            ;;
        --passonly)
            PASSONLY="--passonly"
            shift
            ;;
        --no-collapse-args)
            USE_STRICT_COLLAPSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

if [[ -z "$INPUT_VCF" ]] || [[ -z "$TRUTH_VCF" ]] || [[ -z "$REFERENCE" ]] || [[ -z "$BED_FILE" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required arguments"
    usage
fi

for file in "$INPUT_VCF" "$TRUTH_VCF" "$REFERENCE" "$BED_FILE"; do
    if [[ ! -f "$file" ]]; then
        log_error "File not found: $file"
        exit 1
    fi
done

for tool in truvari bcftools tabix; do
    if ! command -v $tool &> /dev/null; then
        log_error "$tool is not installed or not in PATH"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"
log_info "Output directory: $OUTPUT_DIR"

if [[ ! -f "${REFERENCE}.fai" ]]; then
    log_info "Indexing reference FASTA..."
    samtools faidx "$REFERENCE"
fi

BENCH_VCF="$INPUT_VCF"

if [[ "$RUN_COLLAPSE" == true ]]; then
    log_info "Running truvari collapse to remove redundant variants..."
    
    if [[ ! -f "${INPUT_VCF}.tbi" ]]; then
        log_warning "Input VCF not indexed, attempting to index..."
        if ! tabix -p vcf "$INPUT_VCF" 2>/dev/null; then
            log_warning "Indexing failed (VCF may be unsorted), sorting first..."
            BENCH_VCF=$(sort_and_index_vcf "$INPUT_VCF" "$OUTPUT_DIR")
        else
            BENCH_VCF="$INPUT_VCF"
        fi
    fi
    
    COLLAPSED_VCF="${OUTPUT_DIR}/collapsed.vcf.gz"
    COLLAPSED_SORTED="${OUTPUT_DIR}/collapsed_sorted.vcf.gz"
    
    if [[ "$USE_STRICT_COLLAPSE" == false ]]; then
        log_info "Using relaxed collapse parameters for best F1 score:"
        log_info "  --refdist $COLLAPSE_REFDIST --pctsize $COLLAPSE_PCTSIZE --pctseq $COLLAPSE_PCTSEQ"
        
        truvari collapse \
            --input "$BENCH_VCF" \
            --output "$COLLAPSED_VCF" \
            --reference "$REFERENCE" \
            --refdist $COLLAPSE_REFDIST \
            --pctsize $COLLAPSE_PCTSIZE \
            --pctseq $COLLAPSE_PCTSEQ \
            --removed-output "${OUTPUT_DIR}/removed.vcf"
    else
        log_info "Using default strict collapse parameters"
        
        truvari collapse \
            --input "$BENCH_VCF" \
            --output "$COLLAPSED_VCF" \
            --reference "$REFERENCE" \
            --removed-output "${OUTPUT_DIR}/removed.vcf"
    fi
    
    log_success "Collapse completed"
    
    log_info "Sorting and indexing collapsed VCF..."
    BENCH_VCF=$(sort_and_index_vcf "$COLLAPSED_VCF" "$OUTPUT_DIR")
    
    log_success "Collapsed VCF prepared: $BENCH_VCF"
else
    if [[ ! -f "${INPUT_VCF}.tbi" ]]; then
        log_warning "Input VCF not indexed, attempting to index..."
        if ! tabix -p vcf "$INPUT_VCF" 2>/dev/null; then
            log_warning "Indexing failed (VCF may be unsorted), sorting first..."
            BENCH_VCF=$(sort_and_index_vcf "$INPUT_VCF" "$OUTPUT_DIR")
            log_success "VCF sorted and indexed: $BENCH_VCF"
        fi
    fi
fi

if [[ ! -f "${TRUTH_VCF}.tbi" ]]; then
    log_info "Indexing truth VCF..."
    tabix -p vcf "$TRUTH_VCF"
fi

log_info "Running truvari bench..."
log_info "Parameters:"
log_info "  Input VCF: $BENCH_VCF"
log_info "  Truth VCF: $TRUTH_VCF"
log_info "  Reference: $REFERENCE"
log_info "  BED file: $BED_FILE"
log_info "  refdist: $REFDIST"
log_info "  pctseq: $PCTSEQ"
log_info "  pctsize: $PCTSIZE"

BENCH_OUTPUT="${OUTPUT_DIR}/truvari_bench"

truvari bench \
    -f "$REFERENCE" \
    -b "$TRUTH_VCF" \
    --includebed "$BED_FILE" \
    -c "$BENCH_VCF" \
    -r $REFDIST \
    -p $PCTSEQ \
    -P $PCTSIZE \
    $DUP_TO_INS \
    $PASSONLY \
    -o "$BENCH_OUTPUT"

log_success "Benchmarking completed!"

if [[ -f "${BENCH_OUTPUT}/summary.json" ]]; then
    log_info "Results Summary:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if command -v jq &> /dev/null; then
        TP=$(jq -r '.["TP-base"]' "${BENCH_OUTPUT}/summary.json")
        FP=$(jq -r '.["FP"]' "${BENCH_OUTPUT}/summary.json")
        FN=$(jq -r '.["FN"]' "${BENCH_OUTPUT}/summary.json")
        PRECISION=$(jq -r '.precision' "${BENCH_OUTPUT}/summary.json")
        RECALL=$(jq -r '.recall' "${BENCH_OUTPUT}/summary.json")
        F1=$(jq -r '.f1' "${BENCH_OUTPUT}/summary.json")
        
        printf "  ${GREEN}True Positives:${NC}  %s\n" "$TP"
        printf "  ${RED}False Positives:${NC} %s\n" "$FP"
        printf "  ${RED}False Negatives:${NC} %s\n" "$FN"
        printf "  ${BLUE}Precision:${NC}       %.4f (%.2f%%)\n" "$PRECISION" "$(echo "$PRECISION * 100" | bc -l)"
        printf "  ${BLUE}Recall:${NC}          %.4f (%.2f%%)\n" "$RECALL" "$(echo "$RECALL * 100" | bc -l)"
        printf "  ${BLUE}F1 Score:${NC}        %.4f\n" "$F1"
    else
        cat "${BENCH_OUTPUT}/summary.json"
    fi
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
fi

log_info "Full results available in: $BENCH_OUTPUT"
log_info "Output files:"
log_info "  - Summary: ${BENCH_OUTPUT}/summary.json"
log_info "  - True Positives: ${BENCH_OUTPUT}/tp-comp.vcf.gz"
log_info "  - False Positives: ${BENCH_OUTPUT}/fp.vcf.gz"
log_info "  - False Negatives: ${BENCH_OUTPUT}/fn.vcf.gz"

if [[ "$RUN_COLLAPSE" == true ]]; then
    log_info "  - Collapsed VCF: $BENCH_VCF"
    log_info "  - Removed variants: ${OUTPUT_DIR}/removed.vcf"
fi

log_success "Done!"