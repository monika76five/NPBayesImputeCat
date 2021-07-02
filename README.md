# NPBayesImputeCat
Non-parametric Bayesian Multiple Imputation for Categorical Data

These routines create multiple imputations of missing at random categorical data, with or without structural zeros. Imputations are based on Dirichlet process mixtures of multinomial distributions, which is a non-parametric Bayesian modeling approach that allows for flexible joint modeling.

### Installation

The ```NPBayesImputeCat``` package can be installed from CRAN as follows:

```{r, eval = FALSE}
install.packages("NPBayesImputeCat")
```

The latest version can be installed from GitHub as follows:

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github(repo = "monika76five/NPBayesImputeCat")
```

Please report any bugs you have found under the Issues tab.
