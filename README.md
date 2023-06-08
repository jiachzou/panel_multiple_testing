# Panel Multiple Testing
R and python software for new inference tool in [panel multiple testing](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4315891).  

We conduct inference on panels where the number of units and the number of features for each unit are large. Our novel method method:
- provides Family-Wise Error Rate (FWER) control
- allows _arbitrary_ cross-unit covariance
- allows each unit-level model with arbitrary pattern of missing values, feature selections, etc. so their support sets are varying

Codes are in the directories:
* R/funs.R: Inference for panel hypotheses (to be available soon)
* python/funs.py: Inference for panel hypotheses

## Installation
There is no need for installation, other than dependent libraries, since our code is fully contained in the two scripts respectively for R and python.

- R: no dependencies
- python: pandas, numpy

## R demo
Available soon

## python demo
When there are $N$ units and $J$ features, the evidence of unit-level regressions can be stored in a matrix: 
- a $P$ matrix $J \times N$ of log p-values;
- whenever $P_{jn}$ is missing, the $j$th feature is not in the support set of $n$th unit-level model.

To run the code, we can select features subject to FWER target of $\alpha$:
```
import numpy as np
import pandas as pd

J, N = log_pval_matrix.shape
alpha_vec = [0.00001,0.01,0.05] # the FWER thresholds you want to try
pmt_rejection_table =panel_unordered(log_pval_matrix)
rho=pmt_rejection_table['rho'].unique()[0] # the panel cohesiveness coefficient

for alpha in alpha_vec:

	selected_panel_multiple_testing =np.sort(pmt_rejection_table.index[pmt_rejection_table['rho_inv.N.p_1']<=alpha]).tolist()

	selected_Bonferroni_multiple_testing =np.sort(pmt_rejection_table.index[pmt_rejection_table['p_1']<=alpha/(J*N)]).tolist()

```

## Usage

To cite this code, please use

```
@article{pelger2022inference,
  title={Inference for Large Panel Data with Many Covariates},
  author={Pelger, Markus and Zou, Jiacheng},
  journal={arXiv preprint arXiv:2301.00292},
  year={2022}
}
```