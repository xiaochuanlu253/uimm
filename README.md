# uimm

**uimm** (Universe of Individual-level Mortality Models) is an R package for **mortality modelling** at the individual level.
It provides a unified S4-based framework for fitting, evaluating, and comparing classical and extended mortality models, including parametric, tree-based, and proportional-hazards approaches.

---

## 🌐 Overview

The package integrates actuarial, demographic, and machine learning approaches to model human mortality using individual-level data.  
It supports baseline models, tree partitioning, model-based recursive partitioning, and proportional-hazard extensions, all accessible via a single formula-based interface.

---

## ✨ Features

- **Unified S4 framework** for all parametric age-only models  
- **Built-in mortality laws:**
  - Gompertz, Makeham, Perks, Beard, Kannisto, Thatcher
- **Optimisation-based fitting** via global search (`GenSA`)  
- **Evaluation metrics:** AIC, log-likelihood, integrated Brier score (IBS)  
- **Extended model types:**
  - `BaselineModel`" Age-only PAC
  - `STModel`: Survival tree (LTRCART)
  - `STPModel`: Tree partitioning with PAC at each node
  - `MTPModel`: Model-based recursive partitioning (mob + Gompertz)
  - `GPHModel`: Gompertz proportional hazards (via `flexsurvreg`)
- **Unified formula interface:** one call to `mort_fit()` fits any model

---

## 🧩 Installation

You can install directly from GitHub (recommended during development):

```r
# install.packages("devtools")
devtools::install_github("xiaochuanlu253/uimm")
# load package
library(uimm)
```

---

## ⚙️ Dependencies

The package automatically checks and loads required libraries on initialisation, including:

`tidyverse`, `survival`, `flexsurv`, `GenSA`, `LTRCtrees`, `partykit`, and others.

---

## 🚀 Quick Start

Here is a minimal reproducible workflow:

```r
library(uimm)
set.seed(1)

# Example data
df <- data.frame(
  birthage     = runif(100, 80, 90),
  censoredage  = runif(100, 90, 105),
  DthIndicator = rbinom(100, 1, 0.4)
)

# Fit a baseline Gompertz model
fit <- mort_fit("baseline",
  formula = Surv(birthage, censoredage, DthIndicator) ~ 1,
  data = df, AC = "Gompertz", patience = 30
)

# Show summary
fit
AIC(fit)

# Evaluate integrated Brier score
eval <- evaluate(fit, newdata = df)
eval$IBS
```

---

## 🌳 Combined Models

```r
# Survival Tree Partitioning (STP)
fit_stp <- mort_fit("stp",
  formula = Surv(birthage, censoredage, DthIndicator) ~ Gender + Education,
  data = df, AC = "Beard"
)
fit_stp

# Model-based Recursive Partitioning (MTP)
fit_mtp <- mort_fit("mtp",
  formula = Surv(birthage, censoredage, DthIndicator) ~ 1 |
    Gender + Education + HealthStatus,
  data = df
)
fit_mtp

# Gompertz Proportional Hazards (GPH)
fit_gph <- mort_fit("gph",
  formula = Surv(birthage, censoredage, DthIndicator) ~ Gender + Education,
  data = df
)
fit_gph
```

---

## 📊 Model Evaluation

```r
evaluate(fit, newdata = df,
         metrics = c("AIC", "loglik_sum", "IBS"),
         ibs_method = "predict")

# Plot age-specific IBS contributions
plot(as.numeric(names(eval$IBS_by_age)), eval$IBS_by_age, type = "l")
```

---

## 🧱 S4 Class Overview

| Class           | Description                                                                   |
| --------------- | ----------------------------------------------------------------------------- |
| `PAC`           | Parametric age-only model definition (hazard, survival, cumulative hazard).   |
| `PACFit`        | Fitted PAC object containing parameters and optimisation output.              |
| `BaselineModel` | Single PAC model fitted to all observations.                                  |
| `STModel`       | Survival tree (LTRCART).                                                      |
| `STPModel`      | Survival tree with PAC model per node.                                        |
| `MTPModel`      | Model-based recursive partitioning (mob + Gompertz).                          |
| `GPHModel`      | Gompertz proportional hazards model (flexsurvreg).                            |

---

## 📚 Documentation

All exported functions include detailed help pages with runnable examples.

```r
?mort_fit
?evaluate
?GompertzPAC
```

You can also explore tutorials and extended examples using:

```r
browseVignettes("uimm")
```

---

## 🧪 Development

To regenerate docs, run tests, and build locally:

```r
devtools::document()
devtools::check()
devtools::build_vignettes()
```

---

## 🤝 Citation

If you use **uimm** in research or publications, please cite:

```
Lu, X., Hanewald, K., & Villegas, A. (2025). A Universe of Individual-Level Mortality Models. Available at SSRN 5256615.
```

---

## 🧭 Roadmap

* [ ] Add covariate-extended PAC models
* [ ] Implement parallel fitting for large datasets
* [ ] Add visualisation functions for model diagnostics
* [ ] Extend tree-based models to nonparametric hazards

---

**uimm** brings together actuarial theory, demography, and survival analysis —
offering a unified interface for individual-level mortality modelling.