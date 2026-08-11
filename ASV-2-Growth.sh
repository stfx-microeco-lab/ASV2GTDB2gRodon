#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
source config.yaml

#Grep out your ASVs of interest, grab ids
grep "$organism_of_interest" rawdata/*.tsv  | cut -f1 > intermediate/seqs_$organism_of_interest.ids

#Use ids to get ASV sequences with seqtk
cd $currdir
seqtk subseq rawdata/*.fasta intermediate/seqs_$organism_of_interest.ids > intermediate/seqs_$organism_of_interest.fasta

#Activate BLAST environment
conda activate blast-env

#BLAST ASV sequences against GTDB (bacteria) 16S reference sequences (query coverage set to 100% and identity to 95%)
echo "Identity is set to $perc_identity%. Results will be saved to intermediate folder, with the percentage identity appended to the filename."
blastn -query intermediate/seqs_$organism_of_interest.fasta  -db ${blastDBpath}bac120_ssu_reps.fna  -outfmt 6 -perc_identity $perc_identity -qcov_hsp_perc 100 > intermediate/blastresults_${organism_of_interest}_${perc_identity}

#Obtain GTDB identifiers
cut -f2 intermediate/blastresults_${organism_of_interest}_${perc_identity} > intermediate/GTDB_identifier_${perc_identity}

#Run GTDB_identifiers against GTDB Metadata to obtain genome assembly accession numbers
grep -f intermediate/GTDB_identifier_${perc_identity} ${gtdbDB}bac120_metadata.tsv | grep "Complete Genome" > intermediate/subsetted_metadata_${perc_identity}.tsv

#Run GenBank accession numbers (GCA*_ASM*) against database
grep -f <(awk '{print $66 "_" $58}' intermediate/subsetted_metadata_${perc_identity}.tsv) ${ncbiDBpath}assembly_summary_genbank.txt | cut -f 20 >> intermediate/links_genomes_${organism_of_interest}_${perc_identity}

#Run RefSeq accesion numbers (GCF*_ASM*) against database
grep -f <(awk '{print $66 "_" $58}' intermediate/subsetted_metadata_${perc_identity}.tsv) ${ncbiDBpath}assembly_summary_refseq.txt | cut -f 20 |  >> intermediate/links_genomes_${organism_of_interest}_${perc_identity}

#Edit links to get the genomic data files
awk -F/ -v out="intermediate/links_genomes_${organism_of_interest}_${perc_identity}" '{print $0 $(NF-1) "_genomic.fna.gz" > out}' intermediate/links_genomes_${organism_of_interest}_${perc_identity}


#Use links to download genome sequences from NCBI (raw fastas)
mkdir -p genomes ; cd genomes
wget -i ../intermediate/links_genomes_${organism_of_interest}_${perc_identity}
gunzip -r .

#Activate BAKTA environment
conda activate bakta-2

#Use bakta pipeline to predict ORFs and proteins
for genome in *_genomic.fna
	do base=${genome%%_genomic.fna}
	bakta -d ${BAKTAdb} -o ../bakta-results/${base}/ ${genome} --force
done

#Take output of bakta pipeline and input into gRodon2 to infer predicted maximal growth rate
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

