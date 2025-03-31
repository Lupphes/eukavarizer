# nf-core/eukavarizer: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and is no longer contained in Markdown files._

---

## Introduction

**nf-core/eukavarizer** is a modular and reproducible pipeline built for **structural variant (SV) detection** in **eukaryotic genomes**. It supports both **short and long sequencing reads**, integrates multiple SV calling tools (e.g., **GRIDSS**, **DELLY**, **Sniffles**, **CuteSV**, etc.), and handles reference genome retrieval or uses a local reference.

> **Key features** include:
> - Automated FASTQ detection without a samplesheet
> - Reference-based or reference retrieval from **RefSeq**
> - Quality control using **FastQC**, **MultiQC**, and optional trimming steps
> - SV calling via multiple algorithms
> - Merging and filtering (SURVIVOR, BCFtools)
> - Aggregated reporting (MultiQC, HTML summary)

---

## Input data

Unlike some nf-core pipelines, **eukavarizer** does **not** require a samplesheet. Instead, specify a directory or glob path via `--sequence_dir` (or `--input_dir`) to automatically detect FASTQ files. For example:

```bash
--sequence_dir './fastqs/'
```
or
```bash
--input_dir '/path/to/reads/*.fastq.gz'
```

- **Paired-end** reads should contain `_R1_` / `_R2_` or `.1.fastq.gz` / `.2.fastq.gz`.
- **Single-end** reads can simply match any `.fastq.gz` or `.fq.gz`.
- Multiple FASTQs with the same sample prefix will be merged as replicates.

---

## Running the pipeline

A minimal command:

```bash
nextflow run nf-core/eukavarizer \
  --sequence_dir './fastqs/' \
  --reference_genome './ref/genome.fa.gz' \
  --taxonomy_id 4932 \
  --outdir './results' \
  -profile docker
```

This will:
1. Download or validate the reference genome
2. Preprocess and QC your reads
3. Run the selected SV callers
4. Merge, filter, and report variants

> **Tip**: Provide pipeline parameters via `-params-file` for convenience:

```bash
nextflow run nf-core/eukavarizer \
  -params-file params.yaml \
  -profile docker
```

Where `params.yaml` could be:

```yaml
sequence_dir: "./fastqs/"
reference_genome: "./ref/genome.fa.gz"
taxonomy_id: 4932
outdir: "./results"
gridss_flag: true
delly_flag: true
sniffles_flag: false
cutesv_flag: false
```

---

## Profiles

nf-core pipelines support **profiles** to configure compute environments, container backends, or pipeline behavior. For instance:

- **short_quick**, **short_medium**, **short_full** → short-read data with different depths of analysis
- **long_quick**, **long_medium**, **long_full** → long-read data
- **mix_quick**, **mix_medium**, **mix_full** → hybrid reads
- **docker** → Docker containers
- **conda** or **mamba** → uses Conda/Mamba environments
- **singularity** → Singularity containers
- **wave** → Nextflow Wave containers (≥ 24.03.0-edge)
- **test** → pipeline’s built-in test dataset

Combine them, for example:

```bash
-profile short_full,docker
```

> The **order** matters — later profiles override earlier ones.

---

## Updating the pipeline

Keep the pipeline up to date:

```bash
nextflow pull nf-core/eukavarizer
```

To pin a specific release for reproducibility:

```bash
nextflow run nf-core/eukavarizer -r 1.0.0
```

---

## Reproducibility

For consistent re-runs:

1. Use `-r <version>` to specify pipeline release
2. Use container profiles (e.g., `docker`, `singularity`)
3. Capture your final parameters file (`params.json` / `params.yaml`)
4. Archive `software_versions.yml` and `execution_report.html` from `pipeline_info/`

---

## Running in the background

For lengthy jobs:

```bash
nextflow run ... -bg
```
or wrap with `tmux`, `screen`, or your HPC scheduler system (e.g., SLURM).

---

## Memory settings

The Nextflow JVM can require more memory for large runs. Set an environment variable:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

to allocate between 1–4 GB. Adjust these values to suit your needs.

---

### Questions / Issues

If you encounter problems:

- Check the [nf-core docs](https://nf-co.re/docs)
- Ask on the [Slack `#eukavarizer` channel](https://nfcore.slack.com/channels/eukavarizer)
- File a GitHub Issue in the [nf-core/eukavarizer repository](https://github.com/nf-core/eukavarizer/issues)
