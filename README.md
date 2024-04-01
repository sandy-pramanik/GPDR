# GPDR
Codes to implement Gaussian Process distribution regression (GPDR) as proposed in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) and reproduce figures/results therein.

## Required Package Installation

Before sourcing the GPDR functions, we recommend installing latest versions of required packages `LaplacesDemon`, `ar.matrix`, `tidyverse`, and `doParallel`, `rstan`, `kernlab`, `glmnet`.

## Sourcing GPDR Functions

Specify path to the `GPDR_sourcecode` folder available from this `GitHub` repository. The functions can then be sourced as

``` r
sourcecode.path = ...    # specifies path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"
```

## What Are The Source Codes For?

In the `GPDR_sourcecode` folder,
* `GPDR_Matern.R` contains functions to fit GPDR and output its estimates
* `GPDR_functions.R` sources required libraries and `GPDR_Matern.R`. It also contains functions for setting the true regression function `beta1()` in simulations, providing Kernel Density estimates (KDE), drawing random samples from Dirichlet process (DP), generating data for simulation where independent covariates are simulated from DP and independent and dependent covariates are drawn from a continuous model.

## Reproducing Simulation Results

* Figures 1, 4 and 6 uses empirical posterior risks of GPDR, KDE, and BDR from simulations.

  * Run `SimulationDP_0.1.R` and `SimulationDP_25.R` for simulation results of GPDR and KDE for independent DP-distributed covariates with concentration parameters 0.1 and 25. Change `alpha_DP` in any of the files to get performances for other concentration parameter values.
  * Run `BDR_Simulation.R` for simulation results of BDR for independent DP-distributed covariates with concentration parameters 25. Change `alpha` for `generate_data_dp` function in the file to get performances for other concentration parameter values.

* The empirical CDFs in Figure 2 is obtained by drawing random samples from DP using `rdp()`

`plot_simulation.R` has codes for creating all the figures mentioned above.

* `MWE_DP.R` and `MWE_CP_dep.R` have minimal working examples (MWEs) under independent DP-distributed covariates and dependent continuous covariates. Figures 3 and 5 are generated within these files. `MWE_CP.R` has MWE under independent continuous covariates. It generates a similar figure and is omitted in the manuscript.

* Figures 7 and 8 compare GPDR fits when true data generative model is linear and non-linear.
  * `gp_dist_reg_nonlin_linear_data.R` fits linear and non-linear GPDRs when true data generative model is linear
  * `gp_dist_reg_nonlin.R` fits linear and non-linear GPDRs when true data generative model is non-linear
The files also generates figures within.
