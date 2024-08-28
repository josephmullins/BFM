# Replication code for "Family Law Effects on Divorce, Fertility and Child Investment" by Brown, Flinn and Mullins (2024)

In this repository you will find code to replicate the model estimates, figures and tables from "Family Law Effects on Divorce, Fertility and Child Investment" by Brown, Flinn and Mullins (2024).

## Initial Data Cleaning

The script [`make_BFM_data`](R/make_BFM_data.R) compiles the (almost) raw PSID data in `data/data-psid` and `data/data-cds` into a panel of marriage and labor supply outcomes for parents (`data/MarriagePanel.csv`) and a panel of child skill and investment outcomes (`data/KidPanelv2.csv`) which are used to construct the moments for estimating the model as well as the illustrative facts in the paper.

This script additionally uses the BLS' CPI (`data/CPI-U.csv`) to deflate dollar values to year 2000 USD and merges households with prevailing marital law using (`data/StateCodesDivorce.csv`).

## Figures and Moments

The R notebook [Facts and Moments](`R/FactsMoments.Rmd`) uses the compiled panel datasets to produce the descriptive tables, figures, and regressions used in the text. In addition, it uses IPUMS CPS data (Flood et al. [2023](https://doi.org/10.18128/D030.V11.0)) on custody arrangements (`data/data-cps/cps_00038.csv`) to calculate the auxiliary custody allocation moments used when estimating the model. See the paper for more details.


## Estimation of the Model

The model is estimated in `julia`. All source code can be found in `src`, with code to solve the model in `src/model`, code to run estimation in `src/estimation` and code to run counterfactuals in `src/counterfactuals.jl`.

All results found in the paper are generateed by the following scripts:

- The script `scripts/run_estimation.jl` runs all five stages of the estimation routine. Among these, stage 4 is the most computationally burdensome. The code to solve the model (`src/model/solve.jl`) has been written to make use of multiple threads if they are available.
- The script `scripts/run_bootstrap.jl` samples the data with replacement and re-estimates the model to produce a bootstrap sample of estimates
- The script `scripts/estimate_tables.jl` produces tables of estimates found in the paper
- The script `scripts/run_counterfactuals.jl` runs all the counterfactuals found in the paper, including over the bootstrap sample to provide confidence intervals for all calculations.
