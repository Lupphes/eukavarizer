# nf-core/eukavarizer: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/eukavarizer, created with the [nf-core](https://nf-co.re/) template.

### `Added`

#### Pipeline Features

- Multi-caller structural variant detection supporting 8 SV callers (DELLY, Manta, GRIDSS, SVABA, TIDDIT, DYSGU, Sniffles, CuteSV)
- Support for both short-read and long-read sequencing data
- Automated reference genome retrieval from RefSeq via BioDbCore
- Automated sequencing data retrieval from ENA via BioDbCore
- Optional GATK4/GATK4Spark base quality score recalibration (BQSR)
- Interactive HTML reporting with Varify
- Comprehensive QC reporting with MultiQC

#### Modules (Local)

- `VARIFY`: Unified structural variant reporting and visualization tool
- `BIODBCORE_ENA`: Automated ENA sequencing data retrieval
- `BIODBCORE_REFSEQ`: Automated RefSeq reference genome retrieval
- `SEQKIT_SIZE`: Fast read length calculation
- `DORADO`: Nanopore basecalling from raw signals
- `SVABA_ANNOTATE`: SVABA VCF annotation with StructuralVariantAnnotation
- `SAMPLE_REHEADER`: VCF sample name reheadering
- `GRIDSS_ANNOTATE`: GRIDSS breakend annotation with StructuralVariantAnnotation

#### Modules (nf-core)

- Quality Control: FastQC, MultiQC, Fastp, fastplong, BBDuk, Seqtk
- Alignment: BWA, BWA-MEM2, Minimap2
- BAM/VCF Processing: Samtools (view, merge, index, collatefastq, faidx, stats), BCFtools (concat, sort, merge, stats, filter), GATK4, GATK4Spark, Tabix
- SV Callers: DELLY, Manta, GRIDSS, SVABA, TIDDIT, DYSGU, Sniffles, CuteSV
- SV Processing: SURVIVOR, SVYNC, StructuralVariantAnnotation
- Data Retrieval: SRATools

#### Subworkflows

- `REFERENCE_RETRIEVAL`: Automated reference genome download and indexing
- `SEQUENCE_PROCESSOR`: Complete read processing pipeline (QC, trimming, alignment, BAM processing)
- `QUALITY_CONTROL`: Read quality control and filtering
- `SEQUENCE_ALIGNER`: Multi-aligner support (BWA/BWA-MEM2/Minimap2)
- `SEQUENCE_MERGER`: BAM merging, duplicate marking, and BQSR
- `SV_CALLING_*`: Individual subworkflows for each SV caller (8 total)
- `SV_UNIFICATION`: Cross-caller variant merging and filtering
- `REPORT_GENERATION`: MultiQC and Varify report generation

#### Configuration

- Pre-configured analysis profiles: `short_quick`, `short_medium`, `short_full`, `long_quick`, `long_medium`, `long_full`, `mix_quick`, `mix_medium`, `mix_full`
- Optimized resource allocation for each process
- ENA data retrieval parameters (library strategy, platform, coverage filters)
- Comprehensive parameter validation via nf-schema
- Updated nf-schema plugin to v2.5.1 (latest stable release)
- Added maxErrValSize parameter for better error display control
- Added includeConfig for biodbcore.config to properly expose ENA retrieval parameters
- Removed deprecated minimap2_profile parameter from schema
- Fixed nf-test configuration: changed testsDir from "tests" to "." for proper test discovery
- Added includeConfig for conf/modules.config to load process-specific configurations
- Added includeConfig for conf/biodbcore.config for ENA data retrieval profiles
- Added requires statement to nf-test.config (minimum version 0.9.0)
- Moved ENA retrieval parameters to main nextflow.config with null defaults
- Updated tests/config/nextflow.config with ENA parameter overrides for test isolation

#### Documentation

- Comprehensive README with workflow overview, usage examples, and parameter documentation
- SV caller compatibility matrix
- Enhanced tower.yml for Nextflow Tower/Seqera Platform integration
- Complete tool citations in CITATIONS.md

### `Fixed`

#### Configuration and Testing

- Fixed .nf-core.yml validation error: removed invalid update.exclude list format
- Removed global includeConfig for biodbcore.config to prevent test interference
- ENA retrieval parameters now have null defaults in main config, with specific values in biodbcore.config for relevant profiles

#### Module Modernization

- Added `when` clause for conditional execution to all local modules
- Added `versions.yml` output to all local modules for reproducibility
- Added `task.ext.prefix` and `task.ext.args` support for configurability
- Standardized output patterns using globs instead of hardcoded filenames
- Updated stub sections for all local modules
- Renamed `main.yml` to `meta.yml` for consistency with nf-core standards

#### Meta.yaml Updates

- Updated all local module meta.yaml files with proper input/output documentation
- Added comprehensive descriptions for all modules
- Standardized meta.yaml structure across all modules

#### Test Updates

- Added versions.yml assertions to all local module tests
- Ensured all tests follow nf-core testing standards
- Added nf-test configuration

#### Workflow Headers

- Modernized headers for main EUKAVARIZER workflow
- Updated headers for all 8 SV calling subworkflows
- Updated headers for all 9 sequence processing subworkflows
- Standardized header format with comprehensive documentation

#### Citations

- **Added missing tools**: Tabix/bgzip (Li 2011), SeqKit (Shen et al. 2016)
- **Fixed incorrect citations**:
  - seqtk: Changed to GitHub software citation (no formal publication)
  - fastplong: Changed to GitHub software citation (no separate publication)
  - DYSGU: Corrected authors (Cleal and Baird), journal (Nucleic Acids Research), DOI (10.1093/nar/gkac039)
  - TIDDIT: Corrected DOI (10.12688/f1000research.11168.2)
  - BWA-MEM2: Changed to arXiv citation (arXiv:1907.12931)
  - StructuralVariantAnnotation: Corrected page numbers (2046-2048) and DOI (btac042)
  - BioDbCore: Corrected year (2023)
  - Varify: Corrected year (2024)
- **Organized citations** into logical categories (Basecalling, QC, Alignment, BAM/VCF Processing, SV Calling, etc.)
- Added complete bibliography entries with DOIs and full article titles

### `Dependencies`

#### Core Dependencies

- Nextflow >= 24.04.2
- nf-schema plugin for parameter validation

#### Container/Environment Management

- Docker support
- Singularity/Apptainer support
- Conda/Mamba support
- Podman support

#### Key Software Versions

- Dorado (Oxford Nanopore Technologies)
- FastQC, MultiQC
- Fastp, fastplong, BBDuk, Seqtk, SeqKit
- BWA, BWA-MEM2, Minimap2
- Samtools, BCFtools, GATK4, GATK4Spark, Tabix
- DELLY, Manta, GRIDSS, SVABA, TIDDIT, DYSGU, Sniffles, CuteSV
- StructuralVariantAnnotation, SURVIVOR, SVYNC
- Varify, BioDbCore

(Specific version numbers are defined in module environment.yml files)

### `Deprecated`

- None
