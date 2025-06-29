#!/bin/bash

ID=$1

cellsnp-lite -s 10X_dir/${ID}/outs/possorted_genome_bam.bam \
      -b 10X_dir/${ID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
      -O working_dir/${ID} \
      -R working_dir/${ID}/Heteroplasmy.vcf \
      -p 10 --minMAF 0 --minCOUNT 0 --UMItag Auto --genotype
