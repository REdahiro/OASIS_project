# single-cell deconvolution of LOY

import numpy as np
import pandas as pd
import scanpy as sc

adata=sc.read_h5ad("anndata.h5ad")

df=pd.read_table("genes.table.csv",sep=",")
Y=df[df["chromosome_name"]=="Y"]["external_gene_name"].tolist()
gene=adata.var.index.tolist()
Ygene=(set(Y) & set(gene))
adata.var["Ygene"]=adata.var_names.isin(Ygene)

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["Ygene"], inplace=True, log1p=False
)
adata.obs["LOY"] = adata.obs["pct_counts_Ygene"] == 0

