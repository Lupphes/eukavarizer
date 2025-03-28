qsub -I -l select=1:ncpus=32:mem=128gb:scratch_local=200gb -l walltime=24:00:00
cd $SCRATCHDIR

git clone https://github.com/Lupphes/eukavarizer.git

module add openjdk/17
curl -s https://get.nextflow.io | bash

module add mambaforge
mkdir -p ./.conda_pkgs
mkdir -p ./.conda_envs

export NXF_CONDA_CACHEDIR="$(pwd)/.conda_next"
export NXF_LOG_LEVEL=DEBUG
export NXF_TRACE=true
export NXF_WORK=/storage/brno2/home/luppo/work
export NXF_LOG_FILE=/storage/brno2/home/luppo/logs/.nextflow.log

cd eukavarizer
./../nextflow run main.nf -profile mamba,mix_medium,qc_off -resume
./../nextflow run main.nf -profile mamba,mix_medium,qc_off --taxonomy_id 9606 --reference_genome "/storage/brno2/home/luppo/data/9606/ref/hg38.fa.gz" --sequence_dir "/storage/brno2/home/luppo/data/9606/gib" --outdir /storage/brno2/home/luppo/out -resume



./../nextflow run main.nf -profile mamba,mix_medium,qc_off --taxonomy_id 4932 --reference_genome "/storage/brno2/home/luppo/data/4932/ref/GCF_000146045.2_R64_genomic.fna.gz" --sequence_dir "/storage/brno2/home/luppo/data/4932/ena/" --outdir /storage/brno2/home/luppo/out
