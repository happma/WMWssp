# WMWssp

[![CRANstatus](https://www.r-pkg.org/badges/version/WMWssp)](https://cran.r-project.org/package=WMWssp)
<a href="https://www.rpackages.io/package/WMWssp"><img src="https://www.rpackages.io/badge/WMWssp.svg" /></a>
[![](https://cranlogs.r-pkg.org/badges/WMWssp)](https://cran.r-project.org/package=WMWssp)

Calculates the minimal sample size for the Wilcoxon-Mann-Whitney test that is needed for a given power and two sided type I error rate. The method works for metric data with and without ties, count data, ordered categorical data, and even dichotomous data. But data is needed for the reference group to generate synthetic data for the treatment group based on a relevant effect.
For details, see

Brunner, E., Bathke A. C. and Konietschke, F: Rank- and Pseudo-Rank Procedures in Factorial Designs - Using R and SAS, Springer Verlag, to appear,

Happ, M., Bathke, A. C., & Brunner, E. (2018). Optimal Sample Size Planning for the Wilcoxon-Mann-Whitney-Test. <a href="https://arxiv.org/abs/1805.12249">arXiv preprint arXiv:1805.12249.</a>


To install the current development version:

``` r
## install devtools package if it's not already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("happma/WMWssp")
library(WMWssp)
```

To calculate the sample size we need prior information from one group. Let us call this group reference group. Based on this reference group, we can create artifical data according to an interpretable effect. Note that we have to specify, how many subjects should be assigned to the first and how many to the second group.
``` r
# Prior information for the reference group
x <- c(315,375,356,374,412,418,445,403,431,410,391,475,379)
# generate data for treatment group based on a shift effect
y <- x - 20

# calculate sample size
ssp <- WMWssp(x, y, alpha = 0.05, power = 0.8, t = 1/2)
summary(ssp)
```
It is also possible to vary the allocation rate to even further reduce the sample size.
``` r
# calculate optimal allocation rate t
ssp <- WMWssp_minimize(x, y, alpha = 0.05, power = 0.8)
summary(ssp)
```
