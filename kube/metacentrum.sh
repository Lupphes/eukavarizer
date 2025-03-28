qsub -I -l select=1:ncpus=32:mem=128gb:scratch_local=200gb -l walltime=24:00:00
cd $SCRATCHDIR

git clone https://github.com/Lupphes/eukavarizer.git

module add openjdk/17
curl -s https://get.nextflow.io | bash

module add mambaforge
mkdir -p ./.conda_pkgs
mkdir -p ./.conda_envs

export CONDA_PKGS_DIRS=./.conda_pkgs
export CONDA_ENVS_PATH=./.conda_envs

cd eukavarizer
./../nextflow run main.nf -profile mamba,mix_medium,qc_off -resume
./../nextflow run main.nf -profile mamba,mix_medium,qc_off --taxonomy_id 9606 --reference_genome "/storage/brno2/home/luppo/data/9606/ref/hg38.fa.gz" --sequence_dir "/storage/brno2/home/luppo/data/9606/gib" --outdir /storage/brno2/home/luppo/out
