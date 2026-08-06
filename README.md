# ASV2GTDB2gRodon

Your mission, should you choose to accept it is:

1. grep out your ASVs of interest, grab ids
2. Use ids to get ASV sequences with seqtk
3. BLAST ASV sequences against GTDB (bacteria) 16S reference sequences, I recommend setting query coverage to a high level (maybe up to 100%) and identity to 95%
4. Use matched IDs to get genome sequences from GTDB (raw fastas)
5. Use bakta pipeline (search on google or github) to predict ORFs and proteins
6. Take output of bakta pipeline and input into gRodon2 (also check google or github) to infer predicted maximal growth rate and optimal temperature
