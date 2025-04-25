# Hypothesis-free discovery from epidemiological data by automatic detection and local inference for tree-based nonlinearities and interactions - Supplementary code
_Giorgio Spadaccini (Vrije Universiteit Amsterdam, Leiden University),_

_Marjolein Fokkema (Leiden University),_

_Mark van de Wiel (Vrije Universiteit Amsterdam)_

## Content

This GitHub repository contains the codes run to perform the analyses described in the homonymous paper. It is organized as follows:

- **HELIUS**: folder containing the experiments on the dataset from the HELIUS study
  - _experiment_chol.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for the outcome of cholesterol. Includes codes to plot results
  - _experiment_sbp.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for the outcome of systolic blood pressure. Includes codes to plot results
  - _final_fit_chol.R_: Code to perform one final fit on a larger sample size (n=10000), for the outcome of cholesterol. Includes codes to plot results and analyze interactions
  - _final_fit_sbp.R_: Code to perform one final fit on a larger sample size (n=10000), for the outcome of systolic blood pressure. Includes codes to plot results and analyze interactions
- **Simulations**: folder containing the experiments on simulated data to compare different methods
  - _Performance_p10_gaussian.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for continuous outcome and for p=10 features. Includes codes to plot results
  - _Performance_p30_gaussian.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for continuous outcome and for p=30 features. Includes codes to plot results
  - _Performance_p10_logistic.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for binarized outcome and for p=10 features. Includes codes to plot results
  - _Performance_p30_logistic.R_: Code to perform the (comparative) experiments on smaller sample sizes (n=100,300,500,1000,5000), for binarized outcome and for p=30 features. Includes codes to plot results
  - _TypeI_error.R_: Code to estimate the type I error rate of the different methods discussed in the main text of the paper.
  - _Supplementary_simulations.R_: Code to perform the simulations discussed in the supplementary material but not in the main text of the paper.
  - **Fragmented_simulation**: folder containing an example of how the code may be divided into small batches when the experiment involves many replicates. Covers the example of splitting the runs in the file _TypeI_error.R_
    - _run_experiments_block1.R_: Example of one of the files that need to be created to split the replicates into smaller batches. One file corresponds to one batch.
    - _recombine_and_plot.R_: File that needs to be run last, to recombine the output of each batch and plot the results.
- **RuleSHAP**: folder containing the R package to fit a ruleSHAP model. The package may be installed by using the following snippet of code:
  ```
  install.packages('devtools')
  devtools::install_github('GiorgioSpadaccini/ruleSHAP/ruleSHAP')
  ```



## Example

Here is a quick example of how one may use the ruleSHAP package:

```
#install and load the package
install.packages('ruleSHAP')
library('ruleSHAP')

#simulate the data
data=gendata.friedman1(n=100,p=10)

#define the formula
formula=y~.

#Fit the model
fitted_model=ruleSHAP(formula,data)

#Compute Shapley values
shapley_dfs=compute_SHAP(RS_fit,data[1:3e3,1:p],interactions=T)
marginal_shapley=shapley_dfs$marginal
interaction_shapley=shapley_dfs$interactions
```
