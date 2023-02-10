# File descriptions

* `gauss_ar1_0miss_datasets.rds`: 1000 simulated AR(1) time series with 2 covariates to mimic the basic structure of the GPP dataset. Each time series is stored in a list and named `y`. The parameters used to simulate the data are stored in a sublist called `sim_params` and include `phi` (the AR parameter), `beta` the linear coefficients for the covariates, and `X`, the model matrix for the covariates. 

* `ricker_0miss_datasets.rds`: 1000 simulated integer-valued time series generated using a Ricker population model with Poisson-distributed demographic stochasticity. Each time series is stored in a list and named `y`. The parameters used to simulate the data are stored in a sublist called `sim_params` and include `r`, the intrinsic growth rate, and `alpha`, the intraspecific competition coefficient.