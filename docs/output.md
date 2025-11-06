# nf-core/eukavarizer: Output

## Introduction

This document describes the output produced by **nf-core/eukavarizer**. Most quality control metrics are aggregated in the **MultiQC** report, whereas structural variant (SV) calling results, merges, and final reports are placed in dedicated subfolders. All output directories are placed under the top-level `results/` folder (or path specified with `--outdir`).

---

## Pipeline overview

This pipeline processes data in the following major steps:

1. **Reference Retrieval** (optional, if not using local reference) - Downloads and indexes reference genome
2. **Sequence Retrieval or Local Input** - Downloads from ENA/SRA or processes local FASTQ/BAM/CRAM/FAST5/POD5 files
3. **Basecalling** (optional, for nanopore data) - Dorado basecaller for FAST5/POD5 files
4. **Quality Control** - Pre-filtering FastQC, trimming (fastp/fastplong), contamination removal (BBDuk), post-filtering FastQC, MultiQC reports
5. **Read Alignment** - BWA/BWA-MEM2 for short reads or Minimap2 for long reads
6. **Post-Alignment Processing** (optional) - Duplicate marking, base quality score recalibration (BQSR)
7. **Structural Variant Calling** - Multiple callers: DELLY, Manta, GRIDSS, Dysgu, TIDDIT, SvABA (short reads); Sniffles, CuteSV (long reads)
8. **SV Standardization** - Reheadering and harmonization of VCF files
9. **SV Merging & Filtering** - SURVIVOR and BCFtools merge strategies
10. **Reporting** - Varify HTML report with visualizations, MultiQC aggregate reports

Below are the key result folders, with common file names and descriptions.

> **Note on intermediate files**: The pipeline generates intermediate alignment files (BAM/CRAM) that are kept in Nextflow work directories but not published to the results folder to save disk space. Only final VCF files, QC reports, and reference files are published.

---

## Reference retrieval outputs

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/ref/</code> with subdirectories for indices</summary>

Example structure:

```
results/
└── 4932/
    └── ref/
        ├── GCF_000146045.2_R64_genomic.fna           # Uncompressed FASTA
        ├── GCF_000146045.2_R64_genomic.fna.gz        # BGZipped FASTA
        ├── GCF_000146045.2_R64_genomic.fna.fai       # FASTA index (uncompressed)
        ├── GCF_000146045.2_R64_genomic.fna.gz.fai    # FASTA index (bgzipped)
        ├── GCF_000146045.2_R64_genomic.fna.gz.gzi    # BGZIP index
        ├── bwa/                                       # BWA indices (if using BWA)
        │   ├── *.amb, *.ann, *.bwt, *.pac, *.sa
        ├── bwamem2/                                   # BWA-MEM2 indices (if using BWA-MEM2)
        │   ├── *.amb, *.ann, *.bwt.*, *.pac, *.0123
        └── minimap/                                   # Minimap2 indices (if using Minimap2)
            └── GCF_000146045.2_R64_genomic.mmi
```

**Contents**:

- **Reference FASTA files**:
  - Uncompressed (`.fna`) and BGZipped (`.fna.gz`) versions
  - FASTA indices (`.fai`) for both versions
  - BGZIP index (`.gzi`) for random access to compressed file

- **Alignment indices** (depending on `--bwamem2` and `--minimap2_flag` parameters):
  - `bwa/` - BWA index files (default for short reads)
  - `bwamem2/` - BWA-MEM2 index files (faster variant, use `--bwamem2`)
  - `minimap/` - Minimap2 index files (for long reads, use `--minimap2_flag`)

If a local reference is used (`--reference_fasta`), the pipeline copies or links it into `ref/` for indexing. If using `--taxonomy_id` without a local reference, the pipeline retrieves it from RefSeq via BioDbCore.

</details>

---

## Sequence retrieval outputs

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/ena/</code></summary>

```
results/
└── 4932/
    └── ena/
        ├── SRR123456_1.fastq.gz    # Paired-end read 1
        ├── SRR123456_2.fastq.gz    # Paired-end read 2
        └── SRR789012.fastq.gz      # Single-end or long reads
```

**Contents**:

- **Downloaded reads from ENA/SRA** (if using `--taxonomy_id` with BioDbCore retrieval)
- **Converted reads from samplesheet input** (if using `--input samplesheet.csv`)
  - BAM/CRAM files are converted to FASTQ using `SAMTOOLS_COLLATEFASTQ`
  - SRA files are converted using `SRATOOLS_FASTERQDUMP`

Supported input formats in samplesheet: FASTQ, BAM, CRAM, SRA, FAST5, POD5, BAX_H5

</details>

---

## Basecalling outputs (Oxford Nanopore)

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/dorado/</code></summary>

```
results/
└── 4932/
    └── dorado/
        ├── sample_1.fastq.gz       # Basecalled from FAST5
        ├── sample_2.fastq.gz       # Basecalled from POD5
        └── *.log                   # Dorado log files
```

**Contents** (only if using FAST5 or POD5 input):

- **Basecalled FASTQ files** from Oxford Nanopore raw signal data
- Dorado is used for high-performance basecalling of:
  - `--fast5_dir` - FAST5 format files
  - `--pod5_dir` - POD5 format files

The basecalled FASTQ files are then passed to the QC and alignment steps.

</details>

---

## Quality Control

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/qc/</code> with multiple subfolders</summary>

```
results/
└── 4932/
    └── qc/
        ├── pre_fastqc/          # FastQC before trimming
        ├── pre_multiqc/         # MultiQC aggregating pre-filtering reports
        ├── fastp/               # Fastp trimming (short reads)
        ├── fastplong/           # Fastplong trimming (long reads)
        ├── bbduk/               # BBDuk contamination removal
        ├── seqtk/               # Seqtk read subsampling
        ├── after_fastqc/        # FastQC after trimming
        ├── after_multiqc/       # MultiQC aggregating post-filtering reports
        ├── samtools_stats/      # Alignment statistics
        ├── dedup/               # Deduplication metrics
        └── bqsr/                # Base quality score recalibration
```

**Contents by subdirectory**:

### Pre-filtering QC

- **pre_fastqc/** (if `--fastqc_flag`):
  - `*_fastqc.html` - Interactive HTML QC reports
  - `*_fastqc.zip` - Detailed QC data for each sample

- **pre_multiqc/** (if `--multiqc_flag`):
  - `multiqc_report.html` - Aggregated QC overview across all samples
  - `multiqc_data/` - Raw data from all QC tools in machine-readable format
  - `multiqc_plots/` - Static plot assets

### Trimming and Filtering

- **fastp/** (if `--fastp_flag` for short reads):
  - `*.html` - Fastp HTML reports with QC metrics
  - `*.json` - JSON format statistics (before/after filtering)

- **fastplong/** (if `--fastp_flag` for long reads):
  - `*.html` - Fastplong HTML reports
  - `*.json` - JSON format statistics

- **bbduk/** (if `--bbmap_bbduk_flag`):
  - `*.trim.fastq.gz` - Trimmed/filtered FASTQ files
  - `*.log` - BBDuk contamination removal logs

- **seqtk/** (if `--seqtk_flag`):
  - `*.fastq.gz` - Subsampled FASTQ files
  - Useful for downsampling high-coverage datasets

### Post-filtering QC

- **after_fastqc/** (if `--fastqc_flag`):
  - `*_fastqc.html` - FastQC reports after trimming
  - `*_fastqc.zip` - QC data post-filtering

- **after_multiqc/** (if `--multiqc_flag`):
  - `multiqc_report.html` - Aggregated post-filtering QC
  - `multiqc_data/` - Post-filtering metrics
  - `multiqc_plots/` - Plot assets

### Alignment QC

- **samtools_stats/**:
  - `*.stats` - Comprehensive alignment statistics from Samtools
  - Read mapping rates, insert sizes, error rates, coverage

### Post-alignment Processing

- **dedup/** (if using GATK4 MarkDuplicates):
  - `*.bam` - BAM files with duplicates marked
  - `*.metrics.txt` - Duplication metrics
  - Identifies PCR/optical duplicates

- **bqsr/** (if using GATK4 Base Quality Score Recalibration):
  - `*.table` - BQSR recalibration tables
  - `*.bam` - Recalibrated BAM files
  - Corrects systematic sequencing errors

</details>

---

## Alignment outputs

<details markdown="1">
<summary>Intermediate files (kept in Nextflow work directories)</summary>

> **Important**: Alignment BAM/CRAM files are **NOT published** to the results directory to save disk space. They are kept in Nextflow's work directories and used for SV calling, then can be cleaned up after the pipeline completes.

The pipeline generates the following intermediate alignment files during execution:

**Alignment process**:

1. **BWA/BWA-MEM2** (for short reads) or **Minimap2** (for long reads) aligns reads to reference
2. **SAMtools** sorts and indexes alignments
3. Optionally: **GATK4** marks duplicates and performs BQSR
4. **SV callers** consume these BAM files directly from work directories

**Intermediate directories in work/**:

```
work/
└── [hash]/
    └── intermediate/
        ├── alignments/
        │   ├── bwa/         # BWA alignments
        │   ├── bwamem2/     # BWA-MEM2 alignments
        │   └── minimap2/    # Minimap2 alignments
        ├── samtools_view/   # Format conversions
        ├── samtools_merge/  # Merged BAM files
        └── index/           # BAM indices
```

If you need to retain alignment files, you can:

- Use the `-resume` flag to keep work directories
- Modify `modules.config` to enable publishing for specific alignment processes
- Use `--outdir` to specify a location with sufficient storage

</details>

---

## Structural Variant Calls

<details markdown="1">
<summary>Individual SV caller outputs: <code>{taxonomy_id}/{caller}/</code></summary>

```
results/
└── 4932/
    ├── delly/               # DELLY (short reads)
    ├── manta/               # Manta (short reads)
    ├── gridss/              # GRIDSS (short reads)
    ├── dysgu/               # DYSGU (both read types)
    ├── tiddit/              # TIDDIT (short reads)
    ├── svaba/               # SvABA (short reads, requires BWA)
    ├── sniffles/            # Sniffles (long reads)
    └── cutesv/              # CuteSV (long reads)
```

### Short-read SV callers

**DELLY** (`{taxonomy_id}/delly/`, enabled by `--delly_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF with SV calls
- `*.vcf.gz.tbi` - Tabix index
- `*.vcf.gz.csi` - CSI index
- Detects deletions, duplications, inversions, translocations

**Manta** (`{taxonomy_id}/manta/`, enabled by `--manta_flag`):

- `*.diploidSV.vcf.gz` - Diploid structural variants
- `*.candidateSV.vcf.gz` - Candidate SVs (all evidence)
- `*.candidateSmallIndels.vcf.gz` - Small indels
- `*.tbi` - Tabix indices
- Fast and accurate for germline SVs

**GRIDSS** (`{taxonomy_id}/gridss/`, enabled by `--gridss_flag`):

- `*.vcf.gz` - GRIDSS VCF output
- `*.tbi` - Tabix index
- `annotate/` subdirectory (if `--gridss_annotate`):
  - Annotated SV calls with additional information
- High-resolution breakpoint detection using de Bruijn graphs

**DYSGU** (`{taxonomy_id}/dysgu/`, enabled by `--dysgu_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF
- `*.tbi` - Tabix index
- Works with both short and long reads

**TIDDIT** (`{taxonomy_id}/tiddit/`, enabled by `--tiddit_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF
- `*.tbi` - Tabix index
- Fast SV caller using coverage and discordant pairs

**SvABA** (`{taxonomy_id}/svaba/`, enabled by `--svaba_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF
- `*.tbi` - Tabix index
- `annotate/` subdirectory (optional):
  - Annotated SV calls
- **Note**: Only works with BWA (not BWA-MEM2)
- Uses local assembly for SV detection

### Long-read SV callers

**Sniffles** (`{taxonomy_id}/sniffles/`, enabled by `--sniffles_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF
- `*.tbi` - Tabix index
- Optimized for Pacific Biosciences and Oxford Nanopore data
- Accurate detection of complex SVs

**CuteSV** (`{taxonomy_id}/cutesv/`, enabled by `--cutesv_flag`):

- `*.vcf` - Uncompressed VCF
- `*.vcf.gz` - BGZipped VCF
- `*.tbi` - Tabix index
- Fast and sensitive long-read SV caller

### File formats

All VCF files follow standard VCF 4.2 format with the following typical INFO fields:

- `SVTYPE` - Type of structural variant (DEL, DUP, INV, TRA, INS, BND)
- `SVLEN` - Length of the structural variant
- `END` - End position of the variant
- Caller-specific confidence scores and quality metrics

> **Note**: Directories only appear if the corresponding caller is enabled. Use `--<caller>_flag` parameters to control which callers run. At least one SV caller must be enabled.

</details>

---

## Merged / Filtered SV results

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/results/vcf/</code></summary>

```
results/
└── 4932/
    └── results/
        ├── reheader/                        # Standardized VCF headers
        ├── vcf/                             # Final merged and filtered VCFs
        │   ├── *_survivor_merge.vcf.gz      # SURVIVOR merged
        │   ├── *_survivor_merge.vcf.gz.tbi  # Tabix index
        │   ├── *_survivor_merge_filtered.vcf.gz        # SURVIVOR filtered
        │   ├── *_survivor_merge_filtered.vcf.gz.tbi    # Tabix index
        │   ├── *_survivor_stats.tsv         # SURVIVOR statistics
        │   ├── *_bcftools_concat.vcf.gz     # BCFtools concatenated
        │   ├── *_bcftools_concat.vcf.gz.tbi # Tabix index
        │   ├── *_bcftools_filtered.vcf.gz   # BCFtools filtered
        │   ├── *_bcftools_filtered.vcf.gz.tbi # Tabix index
        │   ├── *_bcftools_stats.txt         # BCFtools statistics
        │   └── sorted/                      # Sorted VCF variants
        ├── merge/                           # Intermediate merge files
        │   └── vcf/
        └── svync/                           # Harmonized SV files (optional)
            └── *.tsv
```

### SV Standardization

**VCF Reheadering** (`results/reheader/`):

- Before merging, all individual caller VCFs are standardized with:
  - `EUK_CALLER` INFO field - Records which SV calling algorithm was used
  - `EUK_PLATFORM` INFO field - Records sequencing platform (Illumina, PacBio, ONT)
  - `EUK_SE` INFO field - Records if data is single-end or paired-end
- Ensures consistent sample names across all VCFs
- Standardized VCFs are passed to merge tools

### SURVIVOR Merge Strategy

**SURVIVOR_MERGE** (`results/vcf/*_survivor_merge.vcf.gz`):

- Merges SVs from multiple callers using breakpoint proximity and SV type
- Parameters controlled by:
  - `--survivor_max_distance` - Maximum distance between breakpoints (default: 1000)
  - `--survivor_min_support` - Minimum number of callers supporting an SV (default: 1)
  - `--survivor_min_size` - Minimum SV size to consider (default: 50)
- Adds INFO fields:
  - `SUPP` - Number of supporting callers
  - `SUPP_VEC` - Binary vector indicating which callers support the variant
  - `SVMETHOD` - Comma-separated list of supporting callers

**SURVIVOR_FILTER** (`results/vcf/*_survivor_merge_filtered.vcf.gz`):

- Applies additional filtering to merged SURVIVOR VCF:
  - Minimum SV size threshold
  - Minimum caller support threshold
  - Quality score filtering
- Produces final high-confidence SV set

**SURVIVOR_STATS** (`results/vcf/*_survivor_stats.tsv`):

- Tab-separated statistics file with:
  - Per-caller SV counts by type
  - Overlap statistics between callers
  - Size distribution summaries

### BCFtools Merge Strategy

**BCFTOOLS_CONCAT** (`results/vcf/*_bcftools_concat.vcf.gz`):

- Concatenates VCFs from all callers
- Maintains individual caller information
- Useful for comparing caller-specific results

**BCFTOOLS_FILTER** (`results/vcf/*_bcftools_filtered.vcf.gz`):

- Applies BCFtools-based filtering:
  - Quality score thresholds
  - Depth filtering
  - Allele frequency filtering
- More flexible than SURVIVOR for custom filters

**BCFTOOLS_STATS** (`results/vcf/*_bcftools_stats.txt`):

- Comprehensive statistics including:
  - SV type counts
  - Length distributions
  - Quality score distributions
  - Per-sample statistics

### SVYNC Harmonization (Optional)

**SVYNC** (`results/svync/*.tsv`):

- Harmonizes SV representations across callers
- Converts VCF to standardized TSV format
- Includes:
  - Chromosome, start, end positions
  - SV type and length
  - Caller information
  - Quality scores
- Useful for downstream analysis and visualization

> **Note**: Both SURVIVOR and BCFtools merge strategies are run in parallel. SURVIVOR is better for combining overlapping calls, while BCFtools preserves individual caller information. Choose the appropriate merged VCF based on your downstream analysis needs.

</details>

---

## Final Reports and Analysis

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/results/</code></summary>

```
results/
└── 4932/
    └── results/
        ├── {taxonomy_id}.html       # Main Varify HTML report
        └── plots/                   # Visualization assets
            ├── *.png                # Static plots
            └── *.html               # Interactive plots
```

### Varify Report

**Main HTML Report** (`results/{taxonomy_id}.html`):

- Comprehensive interactive HTML report with:
  - **Executive Summary**: Key statistics and QC metrics
  - **SV Type Distribution**: Bar charts and pie charts showing SV types (DEL, DUP, INV, etc.)
  - **Size Distribution**: Histograms of SV lengths
  - **Caller Comparison**: Venn diagrams showing overlap between callers
  - **Genomic Distribution**: Circos plots or chromosome ideograms showing SV locations
  - **Quality Metrics**: QC plots from FastQC, MultiQC, and alignment statistics
  - **Filtering Summary**: Before/after filtering statistics
  - **Methods Description**: Automatically generated methods text for publications

**Plots Directory** (`results/plots/`):

- **PNG images**: High-resolution static plots for presentations
- **HTML files**: Interactive Plotly visualizations
- Includes:
  - SV type bar charts
  - Size distribution histograms
  - Caller overlap Venn diagrams
  - Chromosome distribution plots
  - Quality score distributions

### Usage

The Varify report is the primary output for interpreting pipeline results. It aggregates:

- Quality control metrics from all samples
- SV calls from all enabled callers
- Merged and filtered SV results
- Comparison statistics between merge strategies

Open `{taxonomy_id}.html` in a web browser to view the complete analysis.

</details>

---

## Configuration-Dependent Outputs

<details markdown="1">
<summary>Which outputs appear based on pipeline parameters</summary>

The following table shows which output directories are created based on pipeline configuration flags:

| Parameter            | Default  | Output Directory                                                     | Description                    |
| -------------------- | -------- | -------------------------------------------------------------------- | ------------------------------ |
| `--taxonomy_id`      | Required | `{taxonomy_id}/ref/`                                                 | Reference genome retrieval     |
| `--fast5_dir`        | -        | `{taxonomy_id}/dorado/`                                              | Dorado basecalling for FAST5   |
| `--pod5_dir`         | -        | `{taxonomy_id}/dorado/`                                              | Dorado basecalling for POD5    |
| `--fastqc_flag`      | `true`   | `{taxonomy_id}/qc/pre_fastqc/`<br>`{taxonomy_id}/qc/after_fastqc/`   | FastQC reports                 |
| `--multiqc_flag`     | `true`   | `{taxonomy_id}/qc/pre_multiqc/`<br>`{taxonomy_id}/qc/after_multiqc/` | MultiQC aggregated reports     |
| `--fastp_flag`       | `false`  | `{taxonomy_id}/qc/fastp/` or `fastplong/`                            | Read trimming/filtering        |
| `--bbmap_bbduk_flag` | `false`  | `{taxonomy_id}/qc/bbduk/`                                            | Contamination removal          |
| `--seqtk_flag`       | `false`  | `{taxonomy_id}/qc/seqtk/`                                            | Read subsampling               |
| `--bwamem2`          | `false`  | `{taxonomy_id}/ref/bwamem2/`                                         | BWA-MEM2 indices (vs BWA)      |
| `--minimap2_flag`    | `false`  | `{taxonomy_id}/ref/minimap/`                                         | Minimap2 indices (long reads)  |
| `--delly_flag`       | `true`   | `{taxonomy_id}/delly/`                                               | DELLY SV calls                 |
| `--manta_flag`       | `false`  | `{taxonomy_id}/manta/`                                               | Manta SV calls                 |
| `--gridss_flag`      | `true`   | `{taxonomy_id}/gridss/`                                              | GRIDSS SV calls                |
| `--gridss_annotate`  | `false`  | `{taxonomy_id}/gridss/annotate/`                                     | GRIDSS annotations             |
| `--dysgu_flag`       | `false`  | `{taxonomy_id}/dysgu/`                                               | DYSGU SV calls                 |
| `--tiddit_flag`      | `false`  | `{taxonomy_id}/tiddit/`                                              | TIDDIT SV calls                |
| `--svaba_flag`       | `false`  | `{taxonomy_id}/svaba/`                                               | SvABA SV calls (BWA only)      |
| `--sniffles_flag`    | `false`  | `{taxonomy_id}/sniffles/`                                            | Sniffles SV calls (long reads) |
| `--cutesv_flag`      | `false`  | `{taxonomy_id}/cutesv/`                                              | CuteSV SV calls (long reads)   |
| -                    | Always   | `{taxonomy_id}/results/vcf/`                                         | Merged/filtered SV results     |
| -                    | Always   | `{taxonomy_id}/results/{taxonomy_id}.html`                           | Varify HTML report             |
| -                    | Always   | `pipeline_info/`                                                     | Pipeline execution metadata    |

> **Note**: At least one SV caller must be enabled. The default configuration enables DELLY and GRIDSS for a balanced approach to SV detection.

</details>

---

## Pipeline information

<details markdown="1">
<summary>Output directory: <code>pipeline_info/</code></summary>

```
results/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    ├── pipeline_dag.dot
    ├── pipeline_dag.svg
    ├── software_versions.yml
    ├── params.json
    └── ...
```

**These files** help with troubleshooting, provenance, and reproducibility:

- **execution_trace.txt**: table listing processes, resources, and durations
- **execution_report.html** / **execution_timeline.html**: interactive workflow reports from Nextflow
- **pipeline_dag.svg**: directed acyclic graph (DAG) illustrating pipeline steps
- **software_versions.yml**: pinned versions of tools used
- **params.json**: final resolved parameter settings
</details>

[Nextflow trace docs](https://www.nextflow.io/docs/latest/tracing.html) provide more details.

---

## Summary

The **nf-core/eukavarizer** pipeline produces a well-organized output structure:

### Key Output Locations

1. **Reference Files**: `{taxonomy_id}/ref/` - Reference genome and alignment indices
2. **Quality Control**: `{taxonomy_id}/qc/` - Pre/post-filtering QC reports, trimming, deduplication
3. **Individual SV Calls**: `{taxonomy_id}/{caller}/` - Raw VCF files from each enabled caller
4. **Merged Results**: `{taxonomy_id}/results/vcf/` - Unified SV calls using SURVIVOR and BCFtools
5. **Final Report**: `{taxonomy_id}/results/{taxonomy_id}.html` - Interactive Varify HTML report
6. **Pipeline Metadata**: `pipeline_info/` - Execution reports, DAG, software versions

### Important Notes

- **Intermediate alignment files** (BAM/CRAM) are kept in Nextflow work directories, not published to results
- **Output directories** depend on enabled flags - see Configuration-Dependent Outputs section
- **Primary result** is the Varify HTML report, which aggregates all analyses
- **Two merge strategies** (SURVIVOR and BCFtools) are provided - choose based on your needs:
  - SURVIVOR: Better for overlapping SV calls from multiple callers
  - BCFtools: Preserves individual caller information for comparison

### Typical Workflow

1. Start with the **Varify HTML report** (`results/{taxonomy_id}/{taxonomy_id}.html`)
2. Review **MultiQC reports** (`qc/pre_multiqc/` and `qc/after_multiqc/`) for quality assessment
3. Examine **merged VCF files** (`results/vcf/*_survivor_merge_filtered.vcf.gz`) for final SV calls
4. Check **individual caller VCFs** (`{caller}/*.vcf.gz`) if specific validation is needed
5. Use **pipeline_info/** for troubleshooting and reproducibility documentation

For more details on pipeline execution and parameters, see the [usage documentation](usage.md).
