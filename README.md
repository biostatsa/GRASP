
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <u>G</u>uide for <u>R</u>OC-<u>A</u>UC <u>S</u>ample-size <u>P</u>lanning (GRASP)

`GRASP` is an R package designed to estimate the required sample size
for prediction model development based on the Area Under the ROC Curve
(AUC) value. It can be employed for two study designs. In Study 1, the
goal is to determine the required sample size to ensure the length of
the confidence interval (CI) for the empirical AUC estimator remains
within a specific threshold. Study 2 aims to estimate the required
sample size to test a hypothesis regarding AUC, with controlled type I
error and predetermined statistical power levels. Additionally, the
package contains a function for simulation experiments to assess the
performance of various SSD methods for prediction model development.

## Installation

You can install the development version of `GRASP` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("biostatsa/GRASP")
```

## Usage

``` r


library(GRASP)

# Estimate the required sample size based on AUC for Study 1
GRASP_L(theta = 0.8,h = 3,alpha = 0.05, L = 0.1)
#>            group lowerbound upperbound
#> 1     disease:n1         69        246
#> 2 non-disease:n0        207        738
#> 3          all:n        276        984

# Estimate the required sample size based on AUC for Study 2  
GRASP_test(theta0 = 0.75, theta1 = 0.8, h = 3, alpha = 0.05, beta = 0.2)
#>            group lowerbound upperbound
#> 1     disease:n1        159        563
#> 2 non-disease:n0        477       1689
#> 3          all:n        636       2252
```

## Shiny App

We have developed a complementary Shiny app to provide an interactive
interface for the using of `GRASP` to determine the sample size under 2
study designs.

Access the `GRASP` Shiny app here: [GRASP Shiny
App](https://sjbiostat.shinyapps.io/GRASP/)

## Contact

- Author: Yiwang Zhou [ORCID](https://orcid.org/0000-0002-8023-205X)
- Maintainer: Yiwang Zhou <yiwang.zhou@stjude.org>
