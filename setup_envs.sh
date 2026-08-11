#!/bin/bash -i
source config.yaml

#Install MiniForge, saying yes to all recommended options
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -y
#conda init bash
source ~/.bashrc

#Download SeqTK 
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make

#Create directory for databases
mkdir ~/databases

#Create blast environment
conda env create --file envs/blast-env.yaml

#Download and unzip GTDB database (latest)
mkdir -p ~/databases/blast-db; cd ~/databases/blast-db
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz
gunzip ~/databases/blast-db/bac120_ssu_reps.fna.gz

#Make blast database
conda activate blast-env
makeblastdb -dbtype nucl -in ~/databases/blast-db/bac120_ssu_reps.fna

#Download GTDB Metadata TSV file
mkdir ~/databases/GTDB-db/; cd ~/databases/GTDB-db
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz

#Download NCBI databases (RefSeq and GenBank)
mkdir -p ~/databases/NCBI; cd ~/databases/NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

#Download BAKTA and BAKTA_db
conda env create --file envs/bakta-2.yaml
conda activate bakta-2
mkdir -p ~/databases/
bakta_db download --output ~/databases/  --type full

#Create gRodon environment
conda env create --file envs/gRodon.yaml
