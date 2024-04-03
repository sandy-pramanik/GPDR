# GPDR
Codes to implement Gaussian Process distribution regression (GPDR) as proposed in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) and reproduce figures/results therein.

## Required Package Installation

Before sourcing the GPDR functions, we recommend installing latest versions of required packages `LaplacesDemon`, `ar.matrix`, `tidyverse`, `doParallel`, `rstan`, `kernlab`, `glmnet`, `ggplot2`, `dplyr`, and `tidyr`.

## Sourcing GPDR Functions

Specify the path to the `GPDR_sourcecode` folder available from this `GitHub` repository. The functions can then be sourced as

``` r
sourcecode.path = ...    # specifies path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"
```

## What Are The Source Codes For?

In the `GPDR_sourcecode` folder,
 
 * `GPDR_Matern.R` contains functions to fit GPDR and output its estimates
 * `GPDR_functions.R` sources required libraries and `GPDR_Matern.R`. It also contains
   * `GPfit()` for fitting GPDR,
   * `rdp()` for drawing random samples from the Dirichlet process (DP),
   * `generate_data_dp()` for simulating data as in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) where independent covariates are distributed according to DP,
   * `generate_data_cp_indep()` to simulate data as in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) for independently distributed continuous covariates,
   * `generate_data_cp_dep()` for simulating data as in [Tang et al. 2023+](https://arxiv.org/abs/2303.06434) where continuous covariates within each subject are dependent,
   * `beta1()` to set the true regression function in simulations, and
   * `kde_fit()` to provide Kernel Density estimates (KDE).

## Reproducing Simulation Results

* Figures 1, 4, and 6 use empirical posterior risks of GPDR, KDE, and BDR from simulations.
  * Run `SimulationDP.R` for simulation results of GPDR and KDE for independent DP-distributed covariates. Change `alpha_DP` on line 12 there to set the concentration parameter value.
  * Run `BDR_Simulation.R` for simulation results of BDR ([Law et al. 2018](https://proceedings.mlr.press/v84/law18a/law18a.pdf)) for independent DP-distributed covariates with concentration parameters 25. Change `alpha` in `generate_data_dp()` on line 142 to set the concentration parameter value.
  * Run `SimulationCP_indep.R` and `SimulationCP_dep.R` for simulation results under independent and dependent continuous covariates. Change `rhoX` on line 13 in `SimulationCP_dep.R` to control the degree of dependence.
* The empirical CDFs in Figure 2 are obtained by drawing random samples from DP using `rdp()`
* `plot_simulation.R` loads the simulation outputs and creates all the Figures 1, 2, 4, and 6.
* `MWE_DP.R`, `MWE_CP_indep.R`, and `MWE_CP_dep.R` have minimal working examples (MWEs) under independent DP-distributed, and independent and dependent continuous covariates. Figures like 3 and 5 are generated within these files.
* Figures 7 and 8 compare linear and non-linear GPDR when the true data generative model is linear and non-linear.
  * `gp_dist_reg_nonlin_linear_data.R` considers the case where the true data generative model is linear
  * the true data generative model is non-linear in `gp_dist_reg_nonlin.R`
  The files also generate the figures within.
