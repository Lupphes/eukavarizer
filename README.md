<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-eukavarizer_logo_dark.png">
    <img alt="nf-core/eukavarizer" src="docs/images/nf-core-eukavarizer_logo_light.png">
  </picture>
</h1>

<!-- [![GitHub Actions CI Status](https://github.com/nf-core/eukavarizer/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/eukavarizer/actions/workflows/ci.yml) -->
<!-- [![GitHub Actions Linting Status](https://github.com/nf-core/eukavarizer/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/eukavarizer/actions/workflows/linting.yml) -->
<!-- [![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/eukavarizer/results) -->
<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23eukavarizer-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/eukavarizer)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/eukavarizer** is a comprehensive, modular, and reproducible bioinformatics pipeline designed for the detection and analysis of **structural variants (SVs)** in **eukaryotic genomes**. The pipeline supports both **short-read and long-read sequencing data**, integrates multiple state-of-the-art SV callers, and provides unified, high-quality variant calls ready for downstream analysis.

### Key Features

- **Multi-caller SV detection**: Integrates 8 SV callers for comprehensive variant discovery
  - **Short-read optimized**: DELLY, Manta, GRIDSS, SVABA, TIDDIT, DYSGU
  - **Long-read optimized**: Sniffles, CuteSV, DYSGU
- **Automated data retrieval**: Downloads reference genomes from RefSeq and sequencing data from ENA using BioDbCore
- **Complete preprocessing pipeline**: Quality control (FastQC, Fastp, fastplong, BBDuk), adapter trimming, and read filtering
- **Flexible alignment options**: BWA-MEM (default), BWA-MEM2 (faster), or Minimap2 (long-read optimized)
- **Intelligent variant merging**: Unifies results using SURVIVOR and BCFtools with SVYNC standardization
- **Base quality recalibration**: Optional GATK4/GATK4Spark BQSR for improved accuracy
- **Comprehensive reporting**: MultiQC quality reports and interactive HTML summaries with Varify
- **Production-ready**: Follows nf-core best practices for reproducibility and scalability

> **Suitable for diverse eukaryotic organisms**: From single-celled yeast to complex mammalian genomes.

## Workflow Overview

1. **Reference genome retrieval** or usage of user-provided genome
   - Automatic download from RefSeq using BioDbCore
   - Generates BWA/BWA-MEM2/Minimap2 indices
   - Creates bgzipped and indexed FASTA files

2. **Read QC and preprocessing**
   - Quality control with FastQC and MultiQC
   - Adapter trimming and filtering (Fastp for short reads, fastplong for long reads)
   - Optional downsampling (Seqtk) and basecalling (Dorado for Nanopore)

3. **Read alignment**
   - BWA-MEM (default for short reads)
   - BWA-MEM2 (faster alternative)
   - Minimap2 (optimized for long reads)

4. **BAM processing**
   - Multi-lane merging (Samtools)
   - Duplicate marking (GATK4 MarkDuplicates)
   - Optional base quality score recalibration (GATK4 BQSR)

5. **SV calling** via one or more tools:
   - **Short-read**: DELLY, Manta, GRIDSS, SVABA, TIDDIT, DYSGU
   - **Long-read**: Sniffles, CuteSV, DYSGU

6. **VCF standardization**
   - SVYNC normalization
   - StructuralVariantAnnotation for breakpoint refinement

7. **SV merging & filtering**
   - SURVIVOR for cross-caller merging
   - BCFtools for concatenation and filtering

8. **Report generation**
   - MultiQC for aggregated QC metrics
   - Varify for interactive SV visualization
   - BCFtools stats for variant statistics

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

## Quick Start

### Automatic data retrieval (minimal example)

```bash
nextflow run nf-core/eukavarizer \
   -profile docker,short_quick \
   --taxonomy_id 4932 \
   --outdir results/yeast_test
```

This automatically downloads the reference genome for *S. cerevisiae* (taxonomy ID 4932) and retrieves sequencing data from ENA.

### Local input with samplesheet

```bash
nextflow run nf-core/eukavarizer \
   -profile docker,short_full \
   --input samplesheet.csv \
   --reference_genome genome.fa.gz \
   --taxonomy_id 4932 \
   --outdir results/
```

### Long-read analysis

```bash
nextflow run nf-core/eukavarizer \
   -profile docker,long_full \
   --input samplesheet.csv \
   --reference_genome genome.fa.gz \
   --taxonomy_id 4932 \
   --minimap2_flag \
   --outdir results/
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Profiles and Parameters

The pipeline provides several pre-defined profiles to optimise analysis based on input read types and analysis depth. Combine these with a compute environment profile such as `docker`, `conda`, or `mamba`.

### Read Type and Analysis Depth Profiles

| Profile | Read Type | Enabled Callers | Use Case |
|---------|-----------|----------------|----------|
| `short_quick` | Short reads | DELLY, Manta | Quick tests, small genomes |
| `short_medium` | Short reads | DELLY, Manta, GRIDSS, TIDDIT | Balanced analysis |
| `short_full` | Short reads | All short-read callers | Comprehensive detection |
| `long_quick` | Long reads | Sniffles | Quick tests |
| `long_medium` | Long reads | Sniffles, CuteSV | Balanced analysis |
| `long_full` | Long reads | Sniffles, CuteSV, DYSGU | Comprehensive detection |
| `mix_quick` | Hybrid | DELLY, Sniffles | Quick hybrid analysis |
| `mix_medium` | Hybrid | DELLY, Manta, Sniffles, CuteSV | Balanced hybrid |
| `mix_full` | Hybrid | All compatible callers | Maximum sensitivity |

Each profile adjusts:
- Enabled SV callers
- Tool-specific arguments and parameters
- Filtering thresholds for SURVIVOR and BCFtools

### Compute Environment Profiles

- `docker`: Uses Docker containers (recommended)
- `conda` / `mamba`: Uses Conda or Mamba for software installation
- `test`: Runs with minimal test dataset

**Example**: Use `-profile docker,short_full` for a comprehensive analysis on short-read data with Docker.

### Important Parameters

#### Required Parameters

| Parameter | Description |
| --------- | ----------- |
| `--taxonomy_id` | NCBI Taxonomy ID for reference retrieval (e.g., 4932 for *S. cerevisiae*) |
| `--outdir` | Output directory for results |

#### Input Options

| Parameter | Description |
| --------- | ----------- |
| `--input` | Path to samplesheet CSV file |
| `--reference_genome` | Path to reference FASTA file (if not auto-downloading) |

#### Alignment Options

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--bwamem2` | Use BWA-MEM2 instead of BWA-MEM | `false` |
| `--minimap2_flag` | Use Minimap2 for long reads | `false` |

#### SV Caller Flags

| Parameter | Description |
| --------- | ----------- |
| `--delly_flag` | Enable DELLY caller |
| `--manta_flag` | Enable Manta caller |
| `--gridss_flag` | Enable GRIDSS caller |
| `--svaba_flag` | Enable SVABA caller (requires BWA, not BWA-MEM2) |
| `--tiddit_flag` | Enable TIDDIT caller |
| `--dysgu_flag` | Enable DYSGU caller |
| `--sniffles_flag` | Enable Sniffles caller |
| `--cutesv_flag` | Enable CuteSV caller |

For a complete list of parameters, run:
```bash
nextflow run nf-core/eukavarizer --help
```

---

## Pipeline output

The pipeline produces the following main outputs in the specified `--outdir`:

### Main Output Files

| File/Directory | Description |
| -------------- | ----------- |
| `multiqc/multiqc_report.html` | Comprehensive quality control summary across all samples and tools |
| `varify/report.html` | Interactive HTML report with SV visualizations and summary statistics |
| `sv_merged/*.vcf.gz` | Final merged and filtered structural variant calls |
| `sv_calling/*/` | Per-caller VCF files and logs |
| `alignment/*.bam` | Processed BAM files (duplicate-marked, optionally BQSR) |
| `pipeline_info/` | Execution reports, timeline, and resource usage |

**Compatibility Notes**:
- SVABA requires BWA indices (not compatible with `--bwamem2`)
- DYSGU supports both short and long reads
- Minimap2 is recommended for long-read data (enable with `--minimap2_flag`)

## Credits

**nf-core/eukavarizer** was originally written by [Ondřej Sloup (@Lupphes)](https://github.com/Lupphes).

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- Add contributors here -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

### Pipeline Citation

If you use nf-core/eukavarizer for your analysis, please cite:

> **Big Data Analysis: Workflow for Analysing Structural Variants in Sequenced Eukaryotic Genomes**
>
> Ondřej Sloup
>
> <!-- DOI: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) -->

### nf-core Framework

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
