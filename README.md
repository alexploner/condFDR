
# cFDR

This package allows the calculation of the conditional and conjunctional conditional 
false discovery rate (cFDR, ccFDR) from two parallel sets of p-values. The typical
use case is for these two be p-values across SNPs for two different (but presumably
somewhat genetically correlated) phenotypes.

Note: the code is didactic rather than efficient, and has some weird defaults
for p-value pre-filtering for ccFDR.

## Installation

You can install the current version of cFDR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexploner/cFDR")
```

## Example

This is how you can estimate the conditional FDR for p1, conditioned on p2:

``` r
library(cFDR)
data(psynth)
res1 = cFDR(psynth, "p1", "p2", p2_threshold = 1E-5)
head(res1)
head(subset(res1, cFDR < 0.01))
```

This is how you can estimate the conjuctional conditional FDR for p1 and p2, 
corresponding to a false discovery rate for SNPs that are associated with both 
phenotypes:

``` r
res2 = ccFDR(psynth, "p1", "p2", p_threshold = 1E-5)
head(res2)
head(subset(res2, ccFDR < 0.01))
```


