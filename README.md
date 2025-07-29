# Unjustified Poisson assumptions lead to overconfident estimates of the effective reproductive number

This repository contains scripts and figures for a paper, addressing the problem of overdispersion in the [EpiEstim](https://mrc-ide.github.io/EpiEstim/index.html) method Cori et al. (2013) for estimating the effective reproductive number R. We demonstrate that relaxing the Poisson distributional assumption by using rather quasi-Poisson or negative binomial models increase the width of confidence intervals.

## Data
For our case study we analyze 3 datasets that have been used for estimating the effective reproduction number R before

- [Influenza among the active military personnel 2009-2010](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/data/flu/daily_flu.csv), used in [Nash et al. 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439), retrieved from the corresponding [reproducibility repository](https://github.com/rebeccanash/EM_EpiEstim_Nash2023/tree/main)
- [COVID-19 in Austria 2021-2022](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/data/covid/covid_ecdc.csv), used in [Richter et al. 2022](https://www.ages.at/fileadmin/Corona/Epidemiologische-Parameter/updates/Update_Epidemiologische_Parameter_des_COVID19_Ausbruchs_2022-04-01.pdf), retrieved from the [ECDC website](https://www.ecdc.europa.eu/en/publications-data/data-daily-new-cases-covid-19-eueea-country)
- [Ebola in Guinea 2014-2015](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/data/ebola/green2022_data_frame.csv), used in [Green et al. 2022](https://royalsocietypublishing.org/doi/10.1098/rsif.2021.0429), retrieved from the [`stemr` package repository](https://github.com/fintzij/stemr/tree/master)

## Code
File [`analyse_Rt.R`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/scripts/analyse_Rt.R) contains a function performing the analysis on a specified dataset. [`Rt_estimation_all.R`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/scripts/Rt_estimation_all.R) script generates the figures

- [`composite_plot.pdf`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/figure/composite_plot.pdf) 
- [`overdispersion_parameters.pdf`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/figure/overdispersion_parameters.pdf)
- [`nbin_L_vs_qpis_plot.pdf`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/figure/nbin_L_vs_qpis_plot.pdf)
- [`nbin_Q_approximation_plot.pdf`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/figure/nbin_Q_approximation_plot.pdf)

File [`llik_illustration`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/scripts/llik_illustration.R) uses one selected estimation window of the Influenza data to illustrate log-likelihood contributions under different count data models. It produces the [`Incidence_and_llik_bw_identity.pdf`](https://github.com/barbora-sobolova/epiestim_overdispesion/blob/main/figure/Incidence_and_llik_bw_identity.pdf) figure.

Folder [`archived_scripts`](https://github.com/barbora-sobolova/epiestim_overdispesion/tree/main/archived_scripts) and [`archived_data`](https://github.com/barbora-sobolova/epiestim_overdispesion/tree/main/data/archived_data) contains old files from the initial phase of the project.

