# nf-core/eukavarizer: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and is no longer contained in Markdown files._

---

## Introduction

**nf-core/eukavarizer** is a modular and reproducible pipeline built for **structural variant (SV) detection** in **eukaryotic genomes**. It supports both **short and long sequencing reads**, integrates multiple SV calling tools (e.g., **GRIDSS**, **DELLY**, **Sniffles**, **CuteSV**, etc.), and handles reference genome retrieval or uses a local reference.

> **Key features** include:
>
> - Reference-based or reference retrieval from **RefSeq**
> - Quality control using **FastQC**, **MultiQC**, and optional trimming steps
> - SV calling via multiple algorithms
> - Merging and filtering (SURVIVOR, BCFtools)
> - Aggregated reporting (MultiQC, HTML summary)

---

## Input data

**eukavarizer** requires a **samplesheet** as input using the `--input` parameter. The samplesheet is a CSV file that specifies your input data, sample metadata, and sequencing platform information.

### Samplesheet input

Use a CSV samplesheet with `--input`:

```bash
--input samplesheet.csv
```

**Samplesheet format:**

The pipeline uses a samplesheet format that supports multiple input types:

```csv
patient,sex,status,sample,platform,lane,fastq_1,fastq_2,bam,cram,sra,bax.h5,fast5,pod5
HG002,NA,0,sample1,illumina,L001,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,,,,,,
HG002,NA,0,sample2,pacbio,m54238,/path/to/sample2.fastq,,,,,,,
HG003,F,1,sample3,illumina,L001,,,,/path/to/sample3.bam,,,,,
```

**Required columns:**

- `patient` - Patient/subject identifier (can be "NA" or any identifier)
- `sex` - Sex of patient (M/F/NA)
- `status` - Affected status (0=unaffected, 1=affected, NA)
- `sample` - Unique sample identifier (required)
- `platform` - Sequencing platform: `illumina`, `pacbio`, `nanopore` (required)
- `lane` - Sequencing lane identifier (can be any unique identifier per sample)

**Input data columns (provide ONE per row):**

- `fastq_1` - Path to forward/R1 reads or single-end FASTQ
- `fastq_2` - Path to reverse/R2 reads (optional, for paired-end)
- `bam` - Path to BAM file
- `cram` - Path to CRAM file
- `sra` - SRA accession
- `bax.h5` - Path to PacBio BAX.H5 files
- `fast5` - Path to Oxford Nanopore FAST5 files
- `pod5` - Path to Oxford Nanopore POD5 files

**Example samplesheets:**

**Illumina paired-end short reads:**

```csv
patient,sex,status,sample,platform,lane,fastq_1,fastq_2,bam,cram,sra,bax.h5,fast5,pod5
HG002,NA,0,HG002_sample,illumina,L001,/data/HG002_L001_R1.fastq.gz,/data/HG002_L001_R2.fastq.gz,,,,,,
HG002,NA,0,HG002_sample,illumina,L002,/data/HG002_L002_R1.fastq.gz,/data/HG002_L002_R2.fastq.gz,,,,,,
```

**PacBio long reads:**

```csv
patient,sex,status,sample,platform,lane,fastq_1,fastq_2,bam,cram,sra,bax.h5,fast5,pod5
HG002,NA,0,HG002_sample,pacbio,m54238,/data/pacbio_run1.fastq,,,,,,,
HG002,NA,0,HG002_sample,pacbio,m54239,/data/pacbio_run2.fastq,,,,,,,
```

**Mixed sample types:**

```csv
patient,sex,status,sample,platform,lane,fastq_1,fastq_2,bam,cram,sra,bax.h5,fast5,pod5
patient1,F,1,sample1,illumina,L001,/data/short_R1.fastq.gz,/data/short_R2.fastq.gz,,,,,,
patient1,F,1,sample2,nanopore,minion,/data/long.fastq,,,,,,,
patient2,M,0,sample3,illumina,L001,,,,/data/aligned.bam,,,,,
```

**BAM/CRAM input:**

```csv
patient,sex,status,sample,platform,lane,fastq_1,fastq_2,bam,cram,sra,bax.h5,fast5,pod5
patient1,NA,0,sample1,illumina,L001,,,/data/sample1.bam,,,,,,
patient2,NA,0,sample2,illumina,L001,,,,/data/sample2.cram,,,,,
```

**Important notes:**

- Multiple rows with the same `sample` are automatically merged (useful for multi-lane data)
- Only fill in ONE input column per row (fastq_1/fastq_2, bam, cram, sra, etc.)
- Leave empty columns blank (no quotes needed)
- For paired-end FASTQ, provide both `fastq_1` and `fastq_2`
- For single-end FASTQ, only provide `fastq_1`
- `$DATADIR` can be used as a placeholder and replaced using `sed` (see run_scripts examples)

### Automatic data retrieval (auto-generates samplesheet)

If no samplesheet is provided via `--input`, the pipeline automatically retrieves data from ENA/SRA using BioDbCore and generates a samplesheet internally:

```bash
# Pipeline will automatically call BIODBCORE_ENA when --input is not provided
--taxonomy_id 4932  # Saccharomyces cerevisiae
```

Configure the automatic retrieval with:

```bash
--library_strategy WGS                    # Library strategy (default: WGS)
--instrument_platform Illumina            # Platform (default: Illumina)
--minimum_coverage 30                     # Min coverage (default: 40)
--maximum_coverage 100                    # Max coverage (default: 70)
--max_results 5                           # Max samples to retrieve (default: 1)
--assembly_quality ""                     # Assembly quality filter
```

BioDbCore will:

1. Search ENA/SRA for sequencing data matching the taxonomy ID and filters
2. Download the best matching dataset(s)
3. Generate a samplesheet automatically
4. Pass the samplesheet to the pipeline

### Reference genome input

**Local reference:**

```bash
--reference_genome /path/to/genome.fa.gz
```

- Accepts compressed (`.gz`) or uncompressed FASTA
- Automatically indexed for BWA/BWA-MEM2/Minimap2

**Automatic retrieval from RefSeq:**

```bash
--taxonomy_id 4932
# Pipeline downloads reference from NCBI RefSeq
```

### Example usage

**Samplesheet with local reference:**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --reference_genome ./ref/genome.fa \
    --taxonomy_id 4932 \
    --outdir ./results \
    -profile docker
```

**Samplesheet with automatic reference retrieval:**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 9606 \
    --outdir ./results \
    -profile singularity,short_full
```

**Automatic data and reference retrieval (BioDbCore generates samplesheet):**

```bash
# When --input is not provided, BioDbCore automatically retrieves data
nextflow run nf-core/eukavarizer \
    --taxonomy_id 4932 \
    --outdir ./results \
    -profile docker,short_medium
```

Note: This searches ENA/SRA for data matching the taxonomy ID and automatically generates a samplesheet.

---

## Running the pipeline

A minimal command:

```bash
nextflow run nf-core/eukavarizer \
  --input samplesheet.csv \
  --reference_genome './ref/genome.fa.gz' \
  --taxonomy_id 4932 \
  --outdir './results' \
  -profile docker
```

This will:

1. Process samples from the samplesheet
2. Download or validate the reference genome
3. Preprocess and QC your reads
4. Run the selected SV callers
5. Merge, filter, and report variants

> **Tip**: Provide pipeline parameters via `-params-file` for convenience:

```bash
nextflow run nf-core/eukavarizer \
  -params-file params.yaml \
  -profile docker
```

Where `params.yaml` could be:

```yaml
input: "samplesheet.csv"
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

## Testing the pipeline

### Quick test run

To verify the pipeline is working correctly, use the built-in test profile:

```bash
nextflow run nf-core/eukavarizer -profile docker,test --resume
```

This runs a minimal test with yeast data (taxonomy ID 4932) and should complete in a few minutes.

### Comprehensive testing with nf-test

For developers and contributors, the pipeline includes comprehensive nf-test suites:

**Test all workflows, subworkflows, and local modules:**

```bash
nf-test test workflows/eukavarizer/ subworkflows/local/* modules/local/* tests/default.nf.test --profile docker
```

**Test specific components:**

```bash
# Test main workflow only
nf-test test workflows/eukavarizer/ --profile docker

# Test all subworkflows
nf-test test subworkflows/local/* --profile docker

# Test all local modules
nf-test test modules/local/* --profile docker

# Test a specific module
nf-test test modules/local/varify/ --profile docker
```

**Run tests with verbose output:**

```bash
nf-test test workflows/eukavarizer/ --profile docker --verbose
```

> **Note**: nf-test requires Docker or Singularity containers. Tests validate input/output behavior, parameter handling, and ensure modules work correctly in isolation and integrated workflows.

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

## Example run scripts

The repository includes example batch job scripts in the `run_scripts/` directory for running the pipeline on HPC clusters with PBS/Torque schedulers. These scripts demonstrate best practices for production runs.

### Available scripts

| Script                                   | Description                         | Read Type          | Resources               |
| ---------------------------------------- | ----------------------------------- | ------------------ | ----------------------- |
| `eukavarizer_job_dry_run.sh`             | Minimal test run for validation     | Mixed (short+long) | 64 CPUs, 512GB RAM, 24h |
| `eukavarizer_job_short.sh`               | Full short-read analysis            | Short reads        | 64 CPUs, 1.5TB RAM, 48h |
| `eukavarizer_job_nano.sh`                | Oxford Nanopore long reads          | Nanopore           | 64 CPUs, 768GB RAM, 24h |
| `eukavarizer_job_pac.sh`                 | PacBio long reads                   | PacBio             | 64 CPUs, 768GB RAM, 24h |
| `eukavarizer_job_long_horse.sh`          | Long-read analysis for horse genome | Long reads         | 64 CPUs, 768GB RAM, 24h |
| `eukavarizer_job_short_rice_japonica.sh` | Rice japonica short reads           | Short reads        | Variable                |
| `eukavarizer_job_short_rice_indica.sh`   | Rice indica short reads             | Short reads        | Variable                |
| `eukavarizer_job_benchmark.sh`           | Benchmarking run                    | Mixed              | Variable                |
| `eukavarizer_job_benchmark_hs37d5.sh`    | Benchmark with hs37d5 reference     | Mixed              | Variable                |
| `eukavarizer_job_final.sh`               | Production-ready template           | Configurable       | Variable                |

### Script structure

Each script follows this general pattern:

```bash
#!/bin/bash
#PBS -N job_name
#PBS -l select=1:ncpus=64:mem=512gb:scratch_local=400gb
#PBS -l walltime=24:00:00

# Define paths
DATADIR=/path/to/data
SCRATCH=$SCRATCHDIR

# Load modules
module add openjdk/17
module add mambaforge

# Setup environment
export NXF_WORK=$DATADIR/work
export NXF_HOME="$SCRATCH/.nextflow"

# Run pipeline
nextflow run main.nf \
    -profile mamba,short_full \
    --taxonomy_id 9606 \
    --reference_genome "$DATADIR/ref/genome.fna.gz" \
    --input samplesheet.csv \
    --outdir "$DATADIR/results"
```

### Key features

**Scratch space management:**

- Uses fast local scratch storage (`$SCRATCHDIR`) for temporary files
- Stores final results in permanent storage (`$DATADIR`)
- Automatic cleanup with `clean_scratch`

**Environment setup:**

- Loads required modules (Java, Mamba/Conda)
- Sets Nextflow working directory (`NXF_WORK`)
- Configures conda cache directories
- Enables debug logging (`NXF_LOG_LEVEL=DEBUG`)

**Resource allocation:**

- Requests appropriate CPUs, memory, and scratch space
- Sets realistic walltime limits
- Email notifications on job start/finish/abort

**Samplesheet handling:**

- Uses `sed` to replace `$DATADIR` placeholders in samplesheets
- Formats samplesheets for the current environment

### Adapting for your environment

To use these scripts on your HPC cluster:

1. **Modify PBS directives** for your scheduler (SLURM, SGE, etc.):

   ```bash
   # SLURM example
   #SBATCH --job-name=eukavarizer
   #SBATCH --cpus-per-task=64
   #SBATCH --mem=512G
   #SBATCH --time=24:00:00
   ```

2. **Update paths** to match your filesystem:

   ```bash
   DATADIR=/your/storage/path
   ```

3. **Adjust module names** for your system:

   ```bash
   module load java/17
   module load anaconda3
   ```

4. **Select appropriate profile**:
   - Use `-profile singularity` if Docker isn't available
   - Use `-profile conda` or `-profile mamba` for Conda environments
   - Combine with tool profiles: `short_full`, `long_full`, `mix_medium`

5. **Configure resource requests** based on your data:
   - Small genomes (yeast): 64GB RAM, 8-16 CPUs
   - Medium genomes (rice): 256GB RAM, 32 CPUs
   - Large genomes (human): 512GB-1.5TB RAM, 64+ CPUs

### Example: Running on SLURM

Convert PBS script to SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=eukavarizer_short
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=logs/eukavarizer_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com

# Load modules
module load java/17
module load nextflow

# Run pipeline
nextflow run nf-core/eukavarizer \
    -profile singularity,short_full \
    --input /data/samplesheet.csv \
    --taxonomy_id 4932 \
    --outdir /data/results \
    -resume
```

Submit with: `sbatch run_eukavarizer.sh`

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

## Advanced usage

### Custom SV caller configurations

Enable or disable specific SV callers using flags:

```bash
# Short-read callers
--delly_flag true          # DELLY (default: true)
--manta_flag false         # Manta (default: false)
--gridss_flag true         # GRIDSS (default: true)
--dysgu_flag false         # DYSGU (default: false)
--tiddit_flag false        # TIDDIT (default: false)
--svaba_flag false         # SvABA (default: false, requires BWA)

# Long-read callers
--sniffles_flag false      # Sniffles (default: false)
--cutesv_flag false        # CuteSV (default: false)
```

**Example - Enable all short-read callers:**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    --delly_flag true \
    --manta_flag true \
    --gridss_flag true \
    --dysgu_flag true \
    --tiddit_flag true \
    -profile docker
```

### Quality control options

Control QC and trimming steps:

```bash
--fastqc_flag true         # Run FastQC (default: true)
--multiqc_flag true        # Run MultiQC (default: true)
--fastp_flag false         # Run fastp trimming (default: false)
--bbmap_bbduk_flag false   # Run BBDuk contamination removal (default: false)
--seqtk_flag false         # Run seqtk subsampling (default: false)
--seqtk_size 1.0           # Subsampling fraction (1.0 = no subsampling)
```

**Example - Full QC pipeline:**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    --fastqc_flag true \
    --fastp_flag true \
    --bbmap_bbduk_flag true \
    --multiqc_flag true \
    -profile docker
```

### Alignment options

Choose alignment algorithm:

```bash
--bwamem2 false           # Use BWA-MEM2 instead of BWA (default: false)
--minimap2_flag false     # Use Minimap2 for long reads (default: false)
```

**For short reads (BWA/BWA-MEM2):**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    --bwamem2 true \
    -profile docker,short_full
```

**For long reads (Minimap2):**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet_long.csv \
    --taxonomy_id 4932 \
    --minimap2_flag true \
    --sniffles_flag true \
    --cutesv_flag true \
    -profile docker,long_full
```

### SV merging parameters

Configure SURVIVOR merge behavior:

```bash
--survivor_max_distance 1000   # Max breakpoint distance (bp)
--survivor_min_support 1       # Min callers supporting SV
--survivor_min_size 50         # Min SV size (bp)
```

**Example - Strict merging (high confidence):**

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    --survivor_min_support 2 \
    --survivor_min_size 100 \
    -profile docker
```

### Resource configuration

Override default resources for specific processes:

**Via command line:**

```bash
--max_cpus 32
--max_memory 128.GB
--max_time 48.h
```

**Via custom config file:**

Create `custom.config`:

```groovy
process {
    withName: GRIDSS {
        cpus = 16
        memory = 64.GB
        time = 12.h
    }
    withName: DELLY {
        cpus = 8
        memory = 32.GB
    }
}
```

Run with:

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    -profile docker \
    -c custom.config
```

### Output directory structure

Customize output organization:

```bash
--outdir ./results                  # Main output directory
--publish_dir_mode symlink          # symlink, copy, or rellink
--taxonomy_id 4932                  # Creates {outdir}/{taxonomy_id}/
```

### Resume failed runs

If a pipeline run fails or is interrupted:

```bash
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    -profile docker \
    -resume
```

The `-resume` flag reuses cached results from previous runs, only re-running failed or new processes.

### Debugging options

For troubleshooting:

```bash
# Enable verbose logging
export NXF_LOG_LEVEL=DEBUG

# Trace execution
export NXF_TRACE=true

# Custom log file
export NXF_LOG_FILE=./pipeline.log

# Run with debug output
nextflow run nf-core/eukavarizer \
    --input samplesheet.csv \
    --taxonomy_id 4932 \
    -profile docker \
    -with-report report.html \
    -with-trace trace.txt \
    -with-timeline timeline.html \
    -with-dag flowchart.png
```

### Clean up work directory

After successful completion, clean up intermediate files:

```bash
# Remove all work files for a specific run
nextflow clean -f

# Remove work files older than a certain date
nextflow clean -before 2024-01-01

# List what would be removed (dry run)
nextflow clean -n
```

### Pipeline-specific profiles

Use pre-configured tool combinations:

**Short-read profiles:**

- `short_quick` - DELLY + GRIDSS only
- `short_medium` - DELLY + GRIDSS + Manta
- `short_full` - All short-read callers

**Long-read profiles:**

- `long_quick` - Sniffles only
- `long_medium` - Sniffles + CuteSV
- `long_full` - All long-read callers

**Mixed profiles:**

- `mix_quick` - Minimal callers for both
- `mix_medium` - Moderate callers for both
- `mix_full` - All callers

**QC profiles:**

- `qc_off` - Disable QC steps (for testing)

**Example combinations:**

```bash
# Quick short-read analysis
-profile docker,short_quick

# Full long-read analysis with QC disabled
-profile docker,long_full,qc_off

# Mixed data with medium caller set
-profile singularity,mix_medium
```

---

## Questions / Issues

If you encounter problems:

- Check the [nf-core docs](https://nf-co.re/docs)
- Ask on the [Slack `#eukavarizer` channel](https://nfcore.slack.com/channels/eukavarizer)
- File a GitHub Issue in the [nf-core/eukavarizer repository](https://github.com/nf-core/eukavarizer/issues)
