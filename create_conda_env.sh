#!/bin/bash
#SBATCH --job-name=pyRenv
#SBATCH --output=pyEnv.txt
#SBATCH --time=30-99:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB

cd /group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/conda_envs
conda create -p pyRenv python=3.11.8 r-base=4.4 r-essentials=4.4 -y
conda activate /group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/conda_envs/pyRenv

conda install -c conda-forge igraph r-png r-reticulate r-xml r-leidenbase r-rcurl r-rcurl zlib liblzma-devel gcc libgcc cairo r-cairo r-ggrastr r-usethis r-lme4 r-nloptr r-gert r-sf hdf5 meson ninja hatchling hatch-vcs cmake arrow-cpp anndata leidenalg gxx boost -y

conda install -c bioconda bioconductor-rtracklayer bioconductor-cner gtfparse bedtools htslib harmonypy -y

pip install pyarrow==14.0.2 --only-binary :all:
pip install umap-learn diopy numpy jupyter ipykernel ipython decoupler pyarrow pybedtools

R CMD installPackages.R
