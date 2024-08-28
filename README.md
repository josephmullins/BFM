# Replication code for "Family Law Effects on Divorce, Fertility and Child Investment" by Brown, Flinn and Mullins (2024)

## Introduction

In this repository you will find code to replicate the estimates, figures and tables

## Initial Data Cleaning

The script [`make_BFM_data`](R/make_BFM_data.R) compiles the (almost) raw PSID data in `data/data-psid` and `data/data-psid` into a panel of marriage and labor supply outcomes for parents (`data/MarriagePanel.csv`) and a panel of child skill and investment outcomes (`data/KidPanelv2.csv`) which are used to construct the moments for estimating the model as well as the illustrative facts in the paper.

This script additionally uses the BLS' CPI (`data/CPI-U.csv`) to deflate dollar values to year 2000 USD and merges households with prevailing marital law using (`data/StateCodesDivorce.csv`).

## Figures and Moments

The R notebook [Facts and Moments](`R/FactsMoments.Rmd`) uses the compiled panel datasets to produce the descriptive tables, figures, and regressions used in the text. In addition, it uses IPUMS CPS data (Flood et al. [2023](https://doi.org/10.18128/D030.V11.0)) on custody arrangements (`data/data-cps/cps_00038.csv`) to calculate the auxiliary custody allocation moments used when estimating the model. See the paper for more details.

### Other Data

- CPS data on custody allocations are stored in `data/data-cps/cps_00038.csv` and compiled into moments used in estimation 

## Estimation of the Model

The model is estimated in `julia`.

##
