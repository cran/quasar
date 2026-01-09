<img src="sticker.svg" align="right" alt="" width="200" />

# quasar

**quasar** is an R package that provides valid inference procedures when multiple quantile regressions are fitted simultaneously.  
It implements the **rank-score–based closed testing** approach proposed in:

> **De Santis, F., Vesely, A., and Andreella, A. (2025)**  
> *Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach.*  

## Installation

To install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("angeella/quasar")
```

## Toy Example

```r
library(quasar)
library(quantreg)

# Set the dimension of the covariates
p <- 3
Sigma <- diag(p)

# Simulate data
dat_n <- simulateData(
  n            = 200,
  beta         = 0,
  gamma        = c(0.2, -0.1),
  mu           = 4,
  Sigma        = Sigma,
  sigma.y      = 0.5,
  distribution = "t"
)

# Fit quantile regressions
mod <- rq(
  y ~ X + Z1 + Z2,
  tau  = c(0.1, 0.25, 0.5, 0.75, 0.9),
  data = dat_n
)

# Perform closed testing based on rank-score statistics
closedTesting(mod, X = "X")
```

## Simulation Code

The `R` scripts used to reproduce the simulation studies (Figures 2, 3, 4 and 5 in the paper)  
are available in the [`simulations/`](https://github.com/angeella/quasar/tree/main/simulations) folder of this repository.

## Bugs and Feedback

Did you find some bugs?  
Please write to **angela.andreella[at]unive[dot]it**  
or open an issue on the [GitHub Issues page](https://github.com/angeella/quasar/issues),  
preferably including a reproducible example created with the [`reprex`](https://reprex.tidyverse.org/) package.

## License

GPL (≥ 3)

## Citation

If you use **quasar** in your research, please cite:

> De Santis, F., Vesely, A., and Andreella, A. (2026).  
> *Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach.*  
