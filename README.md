# ASV2GTDB2gRodon

Table of contents:

1. [Quickstart](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#1-quickstart)
2. [Overview of pipeline steps](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#2-overview-of-pipeline-steps)
3. [Expected output files](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#3-expected-output-files)
4. [Setting up conda](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#4-setting-up-conda)
5. [Cloning repo and adding raw data](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#5-cloning-repo-and-adding-raw-data)
6. [Downloading tools and databases](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#6-downloading-tools-and-databases)
7. [Running pipeline](https://github.com/stfx-microeco-lab/ASV2GTDB2gRodon#7-running-pipeline)


## 1. Quickstart

The following instructions assume you have git installed.

1. First, clone the repo and `cd` into the directory:
```
git clone git@github.com:stfx-microeco-lab/ASV2GTDB2gRodon.git
cd ASV2GTDB2gRodon
```

2. Here is a script to quickly setup databases and the pipeline's dependecies. Assuming you do not have miniforge installed, it will install miniforge and download the necessary databases for you. These databases are very large, please check if you already have them available in a shared folder.
```
./setup_env_db.sh
```

3. Modify the configuration `config/tutorial/config.yaml` file to include the correct paths to your databases. Please specify the name of your organism of interest and the identity percentage you want GTDB to use.

4. Run the pipeline using your raw data. Make sure the configuration file contains the path to the original TSV file. You can run it with:

```
./ASV-2-Growth.sh
```

## 2. Overview of pipeline steps

A pipeline to take an ASV table, find related genomes from GTDB, and then predict their growth rates using gRodon

1. Grep for the organism of interest and grab ids
2. Use ids to get ASV sequences with seqtk
3. BLAST ASV sequences against GTDB (bacteria) 16S reference sequences
4. Use matched IDs to get genome sequences from NCBI (raw fastas)
5. Use bakta pipeline to predict ORFs and proteins
6. Take output of bakta pipeline and input into gRodon2 to infer predicted maximal growth rate

## 3. Expected output files

The expected output should look like the following:

```
#> There were 8290 genes either with lengths not multiples of 3 or not above
#> length threshold (default 240bp), these genes have been ignored
#> Warning in filterSeq(genes = genes, highly_expressed = highly_expressed, :
#> There were 8290 genes either with lengths not multiples of 3 or not above
#> length threshold (default 240bp), these genes have been ignored
#> $CUBHE
#> [1] 0.8409072
#>
#> $GC
#> [1] 0.5800937
#>
#> $GCdiv
#> [1] 0.08009369
#>
#> $ConsistencyHE
#> [1] 0.636041
#>
#> $CUB
#> [1] 0.5100059
#>
#> $CPB
#> [1] NA
#>
#> $FilteredSequences
#> [1] 8290
#>
#> $nHE
#> [1] 106
#>
#> $dCUB
#> [1] -0.6488187
#>
#> $d
#> [1] 0.4878424
#>
#> $LowerCI
#> [1] 0.387939
#>
#> $UpperCI
#> [1] 0.6181996
```

However, the actual output is looking like this:

```
Loading required package: IRanges
Loading required package: XVector
Loading required package: Seqinfo

Attaching package: ‘Biostrings’

The following object is masked from ‘package:base’:

    strsplit

> ffn_files <- Sys.glob("bakta-results/*/*.ffn")
> genes <- readDNAStringSet(ffn_files)
> cds_files <- Sys.glob("currdir/bakta-results/*_CDS_names.txt")
> CDS_IDs <- lapply(cds_files, readLines)
> gene_IDs <- gsub(" .*","",names(genes))
> genes <- genes[gene_IDs %in% CDS_IDs]
> highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
> results <- predictGrowth(genes, highly_expressed)
Error in getCodonStatistics(genes = genes, highly_expressed = highly_expressed,  :
  No highly expressed genes?
Calls: predictGrowth -> getCodonStatistics
In addition: Warning message:
In predictGrowth(genes, highly_expressed) :
  Less than 10 highly expressed genes provided, performance may suffer
Execution halted
```

## 4. Setting up Conda

```
#Install MiniForge, saying yes to all recommended options
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -y
#conda init bash
source ~/.bashrc
```
## 5. Cloning repo and adding your raw data

```
git clone git@github.com:stfx-microeco-lab/ASV2GTDB2gRodon.git
cd ASV2GTDB2gRodon
```

Add to the `rawdata` folder your TSV (to obtain ids) and FASTA file (to obtain ASV sequences).

## 6. Downloading tools and databases

### 6.1 Downloading SeqTK to obtain ASV sequences

```
#Download SeqTK
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
```

### 6.2 Downloading BLAST env and database to obtain GTDB identifiers

```
#Create directory for databases
mkdir ~/databases

#Create blast environment
conda env create --file envs/blast-env.yaml

#Download and unzip GTDB database (latest)
mkdir -p ~/databases/blast-db; cd ~/databases/blast-db
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz
gunzip ~/databases/blast-db/bac120_ssu_reps.fna.gz

#Make BLAST database
conda activate blast-env
makeblastdb -dbtype nucl -in ~/databases/blast-db/bac120_ssu_reps.fna
```

### 6.3 Downloading GTDB Metdata TSV file to obtain genome assembly accession numbers

```
#Download GTDB Metadata TSV file
mkdir ~/databases/GTDB-db/; cd ~/databases/GTDB-db
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
```

### 6.4 Downloading NCBI databases to obtain links to genomes

```
#Download NCBI databases (RefSeq and GenBank)
mkdir -p ~/databases/NCBI; cd ~/databases/NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

### 6.5 Downloading BAKTA env and database to predict ORFs and proteins

```
#Download BAKTA and BAKTA_db
conda env create --file envs/bakta-2.yaml
conda activate bakta-2
mkdir -p ~/databases/
bakta_db download --output ~/databases/  --type full
```

### 6.6 Create gRodon environment to infer predicted maximal growth rate

```
#Create gRodon environment
conda env create --file envs/gRodon.yaml
```

## 7. Running pipeline

### Before starting, do not forget to modify the configuration file!!

The code is organized into six steps, following the structure outlined in the Overview of pipeline steps section.

### 7.1 Grep for the organism of interest and grab ids

In my case, the IDs were in the first column of the TSV file. To use a different column, change the 1 in `-f1` to the corresponding column number. An intermediate file will be created and placed in the `intermediate` file.

```
grep "$organism_of_interest" rawdata/*.tsv  | cut -f1 > intermediate/seqs_$organism_of_interest.ids
```

### 7.2 Use ids to get ASV sequences with seqtk

```
cd $currdir
seqtk subseq rawdata/*.fasta intermediate/seqs_$organism_of_interest.ids > intermediate/seqs_$organism_of_interest.fasta
```

### 7.3 BLAST ASV sequences against GTDB (bacteria) 16S reference sequences

```
#Activate BLAST environment
conda activate blast-env

#BLAST ASV sequences against GTDB (bacteria) 16S reference sequences (query coverage set to 100% and identity to 95%)
echo "Identity is set to $perc_identity%. Results will be saved to intermediate folder, with the percentage identity appended to the filename."
blastn -query intermediate/seqs_$organism_of_interest.fasta  -db ${blastDBpath}bac120_ssu_reps.fna  -outfmt 6 -perc_identity $perc_identity -qcov_hsp_perc 100 > intermediate/blastresults_${organism_of_interest}_${perc_identity}

#Obtain GTDB identifiers
cut -f2 intermediate/blastresults_${organism_of_interest}_${perc_identity} > intermediate/GTDB_identifier_${perc_identity}
```

### 7.4 Use matched IDs to get genome sequences from NCBI (raw fastas)

```
#Run GTDB_identifiers against GTDB Metadata to obtain genome assembly accession numbers
grep -f intermediate/GTDB_identifier_${perc_identity} ${gtdbDB}bac120_metadata.tsv | grep "Complete Genome" > intermediate/subsetted_metadata_${perc_identity}.tsv
```
Accession numbers can be found in either GenBank or RefSeq, which is why we run the pipeline against both databases. The code grabs field 66, which contains either the GCA (GenBank) or GCF (RefSeq) assembly accession, and combines it with field 58, which contains the second part of the accession number. These values are joined with an underscore to create the complete accession number.
The completed accession numbers are stored in a temporary file. This file is used to query both databases and obtain the genome links. Once its finished, the temporary file is deleted, and the genome links are concatenated into an intermediate file.

```
#Run GenBank accession numbers (GCA*_ASM*) against database
grep -f <(awk '{print $66 "_" $58}' intermediate/subsetted_metadata_${perc_identity}.tsv) ${ncbiDBpath}assembly_summary_genbank.txt | cut -f 20 >> intermediate/links_genomes_${organism_of_interest}_${perc_identity}

#Run RefSeq accesion numbers (GCF*_ASM*) against database
grep -f <(awk '{print $66 "_" $58}' intermediate/subsetted_metadata_${perc_identity}.tsv) ${ncbiDBpath}assembly_summary_refseq.txt | cut -f 20 |  >> intermediate/links_genomes_${organism_of_interest}_${perc_identity}
```
The genome links are modified to construct the correct urls for the FASTA files (.fna.gz).

```
#Edit links to get the genomic data files
awk -F/ -v out="intermediate/links_genomes_${organism_of_interest}_${perc_identity}" '{print $0 $(NF-1) "_genomic.fna.gz" > out}' intermediate/links_genomes_${organism_of_interest}_${perc_identity}


#Use links to download genome sequences from NCBI (raw fastas)
mkdir -p genomes ; cd genomes
wget -i ../intermediate/links_genomes_${organism_of_interest}_${perc_identity}
gunzip -r .
```

### 7.5 Use bakta pipeline to predict ORFs and proteins

Once the BAKTA pipeline is run, the output for each genome is stored in a separate directory. The directory is named based on the corresponding genome accession. This step might take a while!

```
#Activate BAKTA environment
conda activate bakta-2

#Use bakta pipeline to predict ORFs and proteins
for genome in *_genomic.fna
        do base=${genome%%_genomic.fna}
        bakta -d ${BAKTAdb} -o ../bakta-results/${base}/ ${genome} --force
done
```

### 7.6 Take output of bakta pipeline and input into gRodon2 to infer predicted maximal growth rate

```
cd ../bakta-results/
for gff in */*.gff3
        do base=${gff%.gff3}
        sed -n '/##FASTA/q;p' "$gff" | awk '$3=="CDS"' | awk {'print $9'} | awk 'gsub(";.*","")' | awk 'gsub("ID=","")' > ../intermediate/${base}_CDS_names.txt
done
cd ${currdir}
conda activate gRodon
R --vanilla << 'EOF'
library(gRodon)
library(Biostrings)
ffn_files <- Sys.glob("bakta-results/*/*.ffn")
genes <- readDNAStringSet(ffn_files)
cds_files <- Sys.glob("currdir/bakta-results/*_CDS_names.txt")
CDS_IDs <- lapply(cds_files, readLines)
gene_IDs <- gsub(" .*","",names(genes))
genes <- genes[gene_IDs %in% CDS_IDs]
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
results <- predictGrowth(genes, highly_expressed)
print(results)
EOF
```







