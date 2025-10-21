# nf-core/eukavarizer: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

### Basecalling

- [Dorado](https://github.com/nanoporetech/dorado)

> Oxford Nanopore Technologies (2024). Dorado: High-performance basecaller for Oxford Nanopore sequencing data. https://github.com/nanoporetech/dorado

### Quality Control

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

> Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

> Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19):3047-8. doi: 10.1093/bioinformatics/btw354.

- [fastp](https://doi.org/10.1093/bioinformatics/bty560)

> Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17):i884-i890. doi: 10.1093/bioinformatics/bty560.

> Chen S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta, 2:e107. doi: 10.1002/imt2.107.

- [fastplong](https://github.com/OpenGene/fastplong)

> Chen S. fastplong: Ultra-fast preprocessing and quality control for long-read sequencing data. https://github.com/OpenGene/fastplong

- [BBDuk (BBMap)](https://sourceforge.net/projects/bbmap/)

> Bushnell B. (2014). BBMap: A fast, accurate, splice-aware aligner. https://sourceforge.net/projects/bbmap/

- [seqtk](https://github.com/lh3/seqtk)

> Li H. seqtk: Toolkit for processing sequences in FASTA/Q formats. https://github.com/lh3/seqtk

- [SeqKit](https://doi.org/10.1371/journal.pone.0163962)

> Shen W, Le S, Li Y, Hu F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10):e0163962. doi: 10.1371/journal.pone.0163962.

### Read Alignment

- [BWA](https://doi.org/10.1093/bioinformatics/btp324)

> Li H, Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14):1754–1760. doi: 10.1093/bioinformatics/btp324.

- [BWA-MEM2](https://arxiv.org/abs/1907.12931)

> Vasimuddin M, Misra S, Li H, Aluru S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), pp. 314-324. arXiv:1907.12931.

- [Minimap2](https://doi.org/10.1093/bioinformatics/bty191)

> Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.

### BAM/VCF Processing

- [SAMtools](https://doi.org/10.1093/bioinformatics/btp352)

> Li H, Handsaker B, Wysoker A, et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16):2078–2079. doi: 10.1093/bioinformatics/btp352.

- [BCFtools](https://doi.org/10.1093/gigascience/giab008)

> Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2):giab008. doi: 10.1093/gigascience/giab008.

- [Tabix](https://doi.org/10.1093/bioinformatics/btq671)

> Li H. (2011). Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5):718-719. doi: 10.1093/bioinformatics/btq671.

- [GATK4](https://doi.org/10.1101/gr.107524.110)

> McKenna A, Hanna M, Banks E, et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9):1297-1303. doi: 10.1101/gr.107524.110.

> DePristo MA, Banks E, Poplin R, et al. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics, 43(5):491-498. doi: 10.1038/ng.806.

> Van der Auwera GA, O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media. ISBN: 978-1491975183.

### SV Calling (Short Reads)

- [DELLY](https://doi.org/10.1093/bioinformatics/bts378)

> Rausch T, Zichner T, Schlattl A, Stüzer A, Benes V, Korbel JO. (2012). DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics, 28(18):i333-i339. doi: 10.1093/bioinformatics/bts378.

- [Manta](https://doi.org/10.1093/bioinformatics/btv710)

> Chen X, Schulz-Trieglaff O, Shaw R, et al. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32(8):1220–1222. doi: 10.1093/bioinformatics/btv710.

- [GRIDSS](https://doi.org/10.1101/gr.222109.117)

> Cameron DL, Schröder J, Penington JS, et al. (2017). GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. Genome Research, 27(12):2050-2060. doi: 10.1101/gr.222109.117.

- [SvABA](https://doi.org/10.1101/gr.221028.117)

> Wala JA, Bandopadhayay P, Greenwald NF, et al. (2018). SvABA: genome-wide detection of structural variants and indels by local assembly. Genome Research, 28(4):581–591. doi: 10.1101/gr.221028.117.

- [TIDDIT](https://doi.org/10.12688/f1000research.11168.2)

> Eisfeldt J, Vezzi F, Olason P, Nilsson D, Lindstrand A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6:664. doi: 10.12688/f1000research.11168.2.

### SV Calling (Long Reads)

- [Sniffles](https://doi.org/10.1038/s41592-018-0001-7)

> Sedlazeck FJ, Rescheneder P, Smolka M, et al. (2018). Accurate detection of complex structural variations using single-molecule sequencing. Nature Methods, 15(6):461-468. doi: 10.1038/s41592-018-0001-7.

- [CuteSV](https://doi.org/10.1186/s13059-020-02107-y)

> Jiang T, Liu Y, Jiang Y, et al. (2020). Long-read-based human genomic structural variation detection with cuteSV. Genome Biology, 21(1):189. doi: 10.1186/s13059-020-02107-y.

### SV Calling (Both Read Types)

- [DYSGU](https://doi.org/10.1093/nar/gkac039)

> Cleal K, Baird DM. (2022). Dysgu: efficient structural variant calling using short or long reads. Nucleic Acids Research, 50(9):e53. doi: 10.1093/nar/gkac039.

### SV Processing and Annotation

- [StructuralVariantAnnotation](https://doi.org/10.1093/bioinformatics/btac042)

> Cameron DL, Dong R, Papenfuss AT. (2022). StructuralVariantAnnotation: a R/Bioconductor foundation for a caller-agnostic structural variant software ecosystem. Bioinformatics, 38(7):2046-2048. doi: 10.1093/bioinformatics/btac042.

- [SURVIVOR](https://doi.org/10.1038/ncomms14061)

> Jeffares DC, Jolly C, Hoti M, et al. (2017). Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast. Nature Communications, 8:14061. doi: 10.1038/ncomms14061.

- [SVYNC](https://github.com/nvnieuwk/svync)

> Vannieuwkerke N. (2023). SVYNC: A simple structural variant standardization tool. https://github.com/nvnieuwk/svync

### Data Retrieval and Reporting

- [SRA-Tools](https://doi.org/10.1093/nar/gkq1019)

> Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. Nucleic Acids Research, 39(Database issue):D19–D21. doi: 10.1093/nar/gkq1019.

- [BioDbCore](https://github.com/luppo/biodbcore)

> Sloup O. (2023). BioDbCore: Python package for retrieving sequencing data from ENA and reference genomes from RefSeq. https://github.com/luppo/biodbcore

- [Varify](https://github.com/Lupphes/Varify)

> Sloup O. (2024). Varify: Unified structural variant reporting and visualization tool. https://github.com/Lupphes/Varify

### Utilities

- [GNU Coreutils](https://www.gnu.org/software/coreutils/)

> GNU Project. GNU Coreutils: Core utilities for GNU operating systems. https://www.gnu.org/software/coreutils/

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

> Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

> Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

> da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx179. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

> Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

> Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
