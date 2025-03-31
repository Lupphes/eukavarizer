# nf-core/eukavarizer: Output

## Introduction

This document describes the output produced by **nf-core/eukavarizer**. Most quality control metrics are aggregated in the **MultiQC** report, whereas structural variant (SV) calling results, merges, and final reports are placed in dedicated subfolders. All output directories are placed under the top-level `results/` folder (or path specified with `--outdir`).

---

## Pipeline overview

This pipeline processes data in the following major steps:

1. **Reference Retrieval** (optional, if not using local reference)
2. **Sequence Retrieval or Local Input**
3. **Quality Control** (FastQC, optionally FASTP, BBDuk, etc.)
4. **Read Alignment** (BWA or Minimap2)
5. **Structural Variant Calling** (DELLY, Manta, GRIDSS, Dysgu, TIDDIT, Sniffles, CuteSV, etc.)
6. **SV Merging & Filtering** (SURVIVOR, BCFtools)
7. **Reporting** (MultiQC, final HTML report, pipeline info)

Below are the key result folders, with common file names and descriptions.

---

## Reference retrieval outputs

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/ref/</code> and <code>{taxonomy_id}/ref/minimap</code></summary>

Example structure:
```
results/
└── 4932/
    └── ref/
        ├── GCF_000146045.2_R64_genomic.fna.gz
        ├── GCF_000146045.2_R64_genomic.fna.gz.fai
        ├── GCF_000146045.2_R64_genomic.fna.gz.gzi
        ├── GCF_000146045.2_R64_genomic.mmi       # Minimap2 index
        ├── *.bwt / *.sa / *.ann / *.amb / *.pac  # BWA or BWA-MEM2 indices
        └── minimap/
            ├── GCF_000146045.2_R64_genomic.mmi
```

**Contents**:
- **Reference FASTA** + index files
- **Minimap2**, **BWA**, or **BWA-MEM2** index files
- Possibly compressed (`*.gz`) versions of reference plus `.fai`, `.gzi`

If a local reference is used, the pipeline copies or links it into `ref/` for indexing. If using a remote reference, the pipeline places it here after download.
</details>

---

## Sequence retrieval outputs

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/ena/</code></summary>

```
results/
└── 4932/
    └── ena/
        ├── SRR123456_1.fastq.gz
        ├── SRR123456_2.fastq.gz
        ├── ...
        └── size/
            └── *.tsv  (seqkit file size reports)
```

If the pipeline downloads reads from ENA/SRA, the raw data is placed here. For **local** input (`--sequence_dir`), you may see symbolic links or collated outputs (e.g., from `SAMTOOLS_COLLATEFASTQ`, `SRATOOLS_FASTERQDUMP`).
</details>

---

## Quality Control

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/qc/</code> and subfolders</summary>

```
results/
└── 4932/
    └── qc/
        ├── pre_fastqc/
        ├── after_fastqc/
        ├── fastp/
        ├── fastplong/
        ├── bbduk/
        └── multiqc/
```

**Possible contents**:

- **FastQC**:
  - `*_fastqc.html` and `*_fastqc.zip` for each sample, pre- or post-filtering
- **FASTP**, **Fastplong**, **BBDuk**:
  - Logs/JSON files with trimming stats
- **MultiQC**:
  - `multiqc_report.html`: aggregated QC overview
  - `multiqc_data/`: raw data from each QC tool
  - `multiqc_plots/`: static assets from the HTML report
</details>

---

## Alignment outputs

<details markdown="1">
<summary>Output directories: <code>{taxonomy_id}/ena/bwa</code>, <code>{taxonomy_id}/ena/minimap</code>, <code>{taxonomy_id}/ena/sort</code>, <code>{taxonomy_id}/ena/index</code></summary>

```
results/
└── 4932/
    └── ena/
        ├── bwa/
        ├── minimap/
        ├── sort/
        └── index/
```

**Possible contents**:

- `*.bam`: Aligned reads to reference
- `*.bai` or `*.csi`: Index files for each BAM
- Samtools logs or sorting logs

The pipeline typically sorts and indexes the BAM files for downstream analysis.
</details>

---

## Structural Variant Calls

<details markdown="1">
<summary>Individual SV caller outputs</summary>

```
results/
└── 4932/
    ├── delly/
    ├── manta/
    ├── gridss/
    ├── dysgu/
    ├── tiddit/
    ├── sniffles/
    └── cutesv/
```

Each folder typically contains:

- `*.vcf`, `*.vcf.gz`: raw structural variant calls
- `*.tbi` indexes for compressed VCFs
- Caller logs, intermediate files (e.g., `*.bam` for Manta evidence)

If a caller is disabled via flags (e.g., `--delly_flag false`), that folder might be absent.
</details>

---

## Merged / Filtered SV results

<details markdown="1">
<summary>Output directory: <code>{taxonomy_id}/results/vcf</code></summary>

```
results/
└── 4932/
    └── results/
        └── vcf/
            ├── survivor_merged.vcf.gz
            ├── survivor_merge_filtered.vcf.gz
            ├── bcftools_merged.vcf.gz
            ├── bcftools_filtered.vcf.gz
            ├── *.tbi
            ├── *.stats
            └── ...
```

**Key files**:

- **SURVIVOR_MERGE** / **SURVIVOR_FILTER**:
  - `survivor_merged.vcf.gz` + `.tbi`
  - `survivor_merge_filtered.vcf.gz`
  - `survivor_stats.tsv`
- **BCFTOOLS_MERGE** / **BCFTOOLS_FILTER**:
  - `bcftools_merged.vcf.gz` + `.tbi`
  - `bcftools_filtered.vcf.gz`
  - `bcfmerge_stats.tsv`
- Additional stats, logs, index files

These represent the **unified** SV calls and any filtering performed (minimum size, depth, or caller support).
</details>

---

## Final Reports and Analysis

<details markdown="1">
<summary>Output directories: <code>{taxonomy_id}/results</code> or <code>{taxonomy_id}/cutesv</code>, <code>{taxonomy_id}/varify</code>, etc.</summary>

- **VARIFY** (custom module for final variant report):
  - `varify_report.html`
  - Additional JSON/TSV for reporting
- **SVYNC** (if used):
  - Harmonized structural variant sets, e.g., `*.svync.tsv`

Exact folder names may vary depending on your pipeline configuration and flags (e.g., `sniffles_flag`, `cutesv_flag`).
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

In summary, the **nf-core/eukavarizer** pipeline organizes outputs into subfolders by step or tool. FastQC, MultiQC, and general pipeline logs end up in `qc/` and `pipeline_info/`, while alignment, variant calling, and final SV merges are placed in separate directories. This structure allows you to easily navigate and inspect each stage of your analysis.
