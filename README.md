<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-eukavarizer_logo_dark.png">
    <img alt="nf-core/eukavarizer" src="docs/images/nf-core-eukavarizer_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/eukavarizer/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/eukavarizer/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/eukavarizer/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/eukavarizer/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/eukavarizer/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23eukavarizer-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/eukavarizer)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/eukavarizer** is a modular and reproducible bioinformatics pipeline designed for the detection and analysis of **structural variants (SVs)** in **eukaryotic genomes**. The pipeline supports both **short- and long-read data**, integrates a variety of SV callers (e.g., GRIDSS, DELLY, Sniffles, CuteSV), and unifies the results for further analysis. Key features include:

- Reference genome retrieval from RefSeq or use of a local reference
- Automated FASTQ preprocessing and alignment (BWA or Minimap2)
- Modular SV calling with customizable parameters
- Merging and filtering of SVs using SURVIVOR and BCFtools
- Visual and statistical reporting

> Suitable for diverse eukaryotic organisms, from yeast to mammals.

## Workflow Overview

1. **Reference genome retrieval** or usage of user-provided genome
2. **Read QC and preprocessing**
3. **Mapping** using BWA or Minimap2 (long-read aware)
4. **SV calling** via one or more tools:
   - Short-read: DELLY, GRIDSS, TIDDIT, Manta
   - Long-read: Sniffles, CuteSV
5. **SV merging & filtering** with SURVIVOR and BCFtools
6. **Report generation** (including MultiQC and HTML summaries)


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

## Quick Start

### Minimal example:
```bash
nextflow run nf-core/eukavarizer \
   -profile docker,short_quick \
   --taxonomy_id 4932 \
   --outdir results/
```

### Local input example:
```bash
nextflow run nf-core/eukavarizer \
   -profile docker,short_full \
   --taxonomy_id 4932 \
   --reference_genome ./data/4932/ref/genome.fa.gz \
   --sequence_dir ./data/4932/ena/ \
   --outdir results/
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Profiles and Parameters

The pipeline provides several pre-defined profiles to optimise analysis based on input read types and analysis depth. Combine these with a compute environment profile such as `docker`, `conda`, `mamba`, or `kube`.

### Read Type and Analysis Depth Profiles

- `short_quick`, `short_medium`, `short_full`: For short-read data
- `long_quick`, `long_medium`, `long_full`: For long-read data
- `mix_quick`, `mix_medium`, `mix_full`: For hybrid/mixed sequencing

Each level adjusts:
- Enabled SV callers
- Tool-specific arguments
- Filtering thresholds for SURVIVOR and BCFtools

### Compute Environment Profiles

- `docker`: Uses Docker containers
- `conda` / `mamba`: Uses Conda or Mamba for software installation
- `kube`: For Kubernetes environments

Use `-profile docker,short_full` for a full analysis on short-read data with Docker.

### Important Parameters

| Parameter          | Description                                     |
|-------------------|-------------------------------------------------|
| `taxonomy_id`     | NCBI Taxonomy ID for reference retrieval        |
| `sequence_dir`    | Path to directory with raw sequence files       |
| `reference_genome`| Path to a FASTA file for the reference genome   |
| `outdir`          | Output directory for results                    |

---

## Pipeline output

The pipeline produces the following main outputs:

- `multiqc_report.html`: quality control summary
- `report.html`: merged SVs and summary statistics
- `plots/`: additional per-caller plots and images
- VCF files (per caller and merged)

## Credits

nf-core/eukavarizer was originally written by Ondrej Lupphes Sloup.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#eukavarizer` channel](https://nfcore.slack.com/channels/eukavarizer) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

