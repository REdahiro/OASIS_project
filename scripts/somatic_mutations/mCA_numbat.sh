#!/bin/bash

ID=$1

Rscript numbat_dir/inst/bin/pileup_and_phase.R \
    --label ${ID} \
    --samples ${ID} \
    --bams 10X_dir/${ID}/outs/possorted_genome_bam.bam \
    --barcodes 10X_dir/${ID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outdir phase \
    --eagle Eagle_dir/eagle \
    --gmap Eagle_dir/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    --paneldir 1000G_hg38/ \
    --ncores 10


