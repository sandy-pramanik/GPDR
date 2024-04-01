# GPDR
Codes to implement Gaussian Process distribution regression (GPDR) as proposed in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) and reproduce figures/results therein.

## Required Package Installation

Before sourcing the GPDR functions, we recommend installing latest versions of required packages `LaplacesDemon`, `ar.matrix`, `tidyverse`, and `doParallel` from `GitHub`

``` r
remotes::install_github("LaplacesDemonR/LaplacesDemon")
remotes::install_github("cran/ar.matrix")
remotes::install_github("tidyverse/tidyverse")
remotes::install_github("cran/doParallel")
```
or `CRAN`

``` r
install.packages(c("LaplacesDemon", "ar.matrix", "tidyverse", "doParallel"))
```

## Sourcing GPDR Functions

Specify path to the `GPDR_sourcecode` folder available from this `GitHub` repository. The functions can then be sourced as

``` r
# sourcecode.path = ...    # path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"
```

## What Are The Source Codes For?

In the `GPDR_sourcecode` folder,
* `GPDR_Matern.R` contains functions to fit GPDR and output estimates from it
* `GPDR_functions.R` sources required libraries and `GPDR_Matern.R`. It also contains functions for setting the true regression function `beta1()` in simulations, providing Kernel Density estimates (KDE), drawing random samples from Dirichlet process (DP), generating data for simulation where independent covariates are simulated from DP and independent and dependent covariates are drawn from a continuous model.

