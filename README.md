# Panel Multiple Testing
R and python software for [panel multiple testing](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4315891).  

New tools for panel inference when different units have different support sets. This is designed for panels where the number of units and the number of features for each unit are large. The inference result is a common feature set with selection false discovery control by Family-Wise Error Rate (FWER) guarantee.

Code is in the directory.
* R/funs.R: Inference for panel hypotheses, implemented 
* python/funs.py: Inference for LASSO at a fixed, deterministic value of lambda.

## Installation
There is no need for installation since our code is fully contained in the two scripts respectively for R and python.

## R demo

```

panel_unordered

```
## python demo

```

panel_unordered

```

The simulation code is available upon request.

To cite this code, please use

```
@article{pelger2022inference,
  title={Inference for Large Panel Data with Many Covariates},
  author={Pelger, Markus and Zou, Jiacheng},
  journal={arXiv preprint arXiv:2301.00292},
  year={2022}
}
```