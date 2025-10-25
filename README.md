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

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A5v25.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with nf-core template](https://img.shields.io/badge/nf--core%20template-3.4.1-brightgreen.svg)](https://nf-co.re/tools)
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

The pipeline processes data through these main stages:

1. **Input & Reference** → Retrieve or use provided reference genome; process input data (FASTQ/BAM/CRAM/raw formats)
2. **QC & Alignment** → Quality control, trimming, and alignment (BWA-MEM/BWA-MEM2/Minimap2)
3. **SV Calling** → Multi-caller structural variant detection (DELLY, Manta, GRIDSS, Sniffles, CuteSV, etc.)
4. **Merging & Filtering** → Unify and filter variants using SURVIVOR and BCFtools
5. **Reporting** → Generate interactive Varify report and MultiQC quality metrics

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

## Quick Start

### Test the pipeline

```bash
nextflow run nf-core/eukavarizer -profile docker,test --resume
```

### Automatic data retrieval (minimal example)

```bash
nextflow run nf-core/eukavarizer \
   -profile docker,short_quick \
   --taxonomy_id 4932 \
   --outdir results/yeast_test
```

This automatically downloads the reference genome for *S. cerevisiae* (taxonomy ID 4932) and retrieves sequencing data from ENA using BioDbCore.

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

### HPC/Cluster execution

Example PBS/SLURM job scripts are provided in `run_scripts/` for different scenarios:
- Short-read analysis: `run_scripts/eukavarizer_job_short.sh`
- Long-read (PacBio): `run_scripts/eukavarizer_job_pac.sh`
- Long-read (Nanopore): `run_scripts/eukavarizer_job_nano.sh`

See [docs/usage.md](docs/usage.md) for detailed HPC deployment instructions.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Key Parameters

| Parameter | Description | Required |
| --------- | ----------- | -------- |
| `--taxonomy_id` | NCBI Taxonomy ID (e.g., 4932 for yeast, 9606 for human) | Yes |
| `--outdir` | Output directory | Yes |
| `--input` | Samplesheet CSV file | No* |
| `--reference_genome` | Reference FASTA file | No* |

*If `--input` is not provided, BioDbCore automatically retrieves data from ENA/SRA

### Analysis Profiles

Combine read type with compute environment: `-profile docker,short_full`

| Profile | Read Type | SV Callers |
|---------|-----------|------------|
| `short_quick/medium/full` | Short reads | DELLY, Manta, GRIDSS, TIDDIT, SVABA, DYSGU |
| `long_quick/medium/full` | Long reads | Sniffles, CuteSV, DYSGU |
| `mix_quick/medium/full` | Hybrid | Combined callers |
| `test` | Test data | Minimal yeast test |

**Compute profiles**: `docker` (recommended), `singularity`, `conda`, `mamba`

For all parameters: `nextflow run nf-core/eukavarizer --help` or see [docs/usage.md](docs/usage.md)

## Documentation

- **[Usage Guide](docs/usage.md)** - Detailed instructions, input formats, parameters, and HPC deployment
- **[Output Documentation](docs/output.md)** - Complete output file descriptions and directory structure
- **[Parameter Reference](https://nf-co.re/eukavarizer/parameters)** - Full parameter documentation

For samplesheet format and examples, see [docs/usage.md](docs/usage.md#samplesheet-input).

---

## Pipeline output

Primary outputs in `--outdir`:

- **`{taxonomy_id}/results/{taxonomy_id}.html`** - Varify interactive report with SV visualizations
- **`{taxonomy_id}/qc/after_multiqc/multiqc_report.html`** - Comprehensive QC metrics
- **`{taxonomy_id}/results/vcf/`** - Merged and filtered structural variants (SURVIVOR and BCFtools strategies)
- **`{taxonomy_id}/{caller}/`** - Individual SV caller outputs (DELLY, Manta, GRIDSS, Sniffles, CuteSV, etc.)
- **`{taxonomy_id}/ref/`** - Reference genome and indices
- **`pipeline_info/`** - Execution reports (timeline, trace, software versions)

See [docs/output.md](docs/output.md) for complete output documentation and directory structure

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
