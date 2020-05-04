# An early analysis of the COVID-19 pandemic

## Dataset

This dataset is collected from public agencies or news media,
containing detailed information about some 1400 COVID-19 cases
confirmed in and outside China. This dataset is free to use and share
given that appropriate credits are given under the [CC-BY-4.0
license](https://github.com/qingyuanzhao/bets.covid19/blob/master/LICENSE.md). It
can be loaded in R as a package:
```r
devtools::install_github("qingyuanzhao/bets.covid19")
library(bets.covid19)
head(covid19_data)
```
More details about the dataset can be found in
```r
help(covid19_data)
```
and in [this arXiv preprint](https://arxiv.org/abs/2004.07743).

## Statistical inference: the BETS model

We have developed a generative model for four key epidemiological events: Beginning of exposure, End of exposure, time of Transmission, and time of Symptom onset (BETS). This package implements a likelihood inference for the BETS model. Try:
```r
help(bets.inference)
example(bets.inference)
```
Details of the model and methodology can be found in [this
preprint](https://arxiv.org/abs/2004.07743) on arXiv. In short, we
find that several published early analyses were severely biased by
sample selection. All our analyses, regardless of which subsample and
model were being used, point to **an epidemic doubling time of 2 to
2.5 days** during the early outbreak in Wuhan.

A Bayesian
nonparametric analysis further suggests that **5% of the symptomatic
cases may not develop symptoms within 14 days since infection**. Code
for the Bayesian model and MCMC sampler can be found under the
[bayesian
folder](https://github.com/qingyuanzhao/bets.covid19/tree/master/bayesian).

## Reference

- First report: Qingyuan Zhao, Yang Chen, Dylan S Small. Analysis of the epidemic growth of the early 2019-nCoV outbreak using internationally confirmed cases. medRxiv 2020.02.06.20020941; doi: https://doi.org/10.1101/2020.02.06.20020941
- Full model: Qingyuan Zhao, Niaoqiao Ju, Sergio Bacallado, Rajen Shah. BETS: The dangers of selection bias in early analyses of the coronavirus disease (COVID-19) pandemic. [arXiv:2004.07743](https://arxiv.org/abs/2004.07743).

## Acknowledgement

Many people have contributed to the data collection and given helpful suggestions. We thank Rajen Shah, Yachong Yang, Cindy Chen, Yang Chen, Dylan Small, Michael Levy, Hera He, Zilu Zhou, Yunjin Choi, James Robins, Marc Lipsitch, Andrew Rosenfeld.

## Earlier work

This project first started from a preliminary analysis of some
international COVID-19 cases exported from Wuhan. The report of the
first analysis can be found on
[medRxiv](https://www.medrxiv.org/content/10.1101/2020.02.06.20020941v1). Code
for that analysis can be found in the [report1](https://github.com/qingyuanzhao/2019-nCov-Data/tree/report1) branch.
