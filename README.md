# Panel Multiple Testing
R and python software for [panel multiple testing](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4315891).  

New tools for panel inference when different units have different support sets. This is designed for panels where the number of units and the number of features for each unit are large. The inference result is a common feature set with selection false discovery control by Family-Wise Error Rate (FWER) guarantee.

Code is in the directory.
* R/funs.R: Inference for panel hypotheses, implemented 
* python/funs.py: Inference for LASSO at a fixed, deterministic value of lambda.

## Installation
There is no need for installation since our code is fully contained in the two scripts respectively for R and python.

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