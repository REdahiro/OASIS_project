# single-cell deconvolution of mCAs

library(Seurat)
library(numbat)

data<-readRDS("seurat_obj.rds")
count_mat<-data@assays$RNA@counts
df_allele=read.table("allele_counts.tsv.gz",sep="\t",header=T)
existCNV=read.table("existCNV.txt",sep="\t",header=T)

out = run_numbat(
    count_mat,
    ref_hca,
    df_allele,
    segs_consensus_fix=existCNV,
    genome = "hg38",
    t = 1e-5,
    min_LLR=5,
    ncores = 10,
    init_k=10,
    plot = TRUE,
    max_iter = 2,
    max_entropy=1,
    out_dir = './output'
)

