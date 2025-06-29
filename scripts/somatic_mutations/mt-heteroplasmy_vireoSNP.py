# single-cell deconvolution of mt-heteroplasmy

import vireoSNP
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from vireoSNP import BinomMixtureVB

AD = mmread("cellSNP.tag.AD.mtx").tocsc()
DP = mmread("cellSNP.tag.DP.mtx").tocsc()

_model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=2)
_model.fit(AD, DP, min_iter=30, n_init=50)

with open('cellSNP.samples.tsv', 'r') as f:
    sample_id = f.read().splitlines()

conf = []
for i in _model.ID_prob:
    if max(i) >= 0.8:
        conf.append('True')
    else:
        conf.append('False')

clone_id = np.argmax(_model.ID_prob, axis=1)
clones_df = pd.DataFrame(data={'sample_id':sample_id, 'clone_id':clone_id, 'confident':conf})
clones_df[(clones_df["clone_id"]==0)&(clones_df["confident"]=="True")]
clones_df.to_csv("clone.tsv", index=None, sep="\t")


