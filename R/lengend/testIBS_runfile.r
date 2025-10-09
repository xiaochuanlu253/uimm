# load dataset
load("../data/cvdata.RData")
load("../data/data_sb_nrc.RData")
nhis = data_sb_nrc$df_complete
# load code from age only model implement.r
source("./testIBS.r")



# Baseline model
fit_base <- mort_fit(
  model   = "baseline",
  formula = Surv(birthage, censoredage, DthIndicator) ~ 1,
  data    = cvdata_9800$train,
  data_test = cvdata_9800$test,
  AC = "Beard"
)

# Survival tree
fit_st <- mort_fit(
  model   = "st",
  formula = Surv(birthage, censoredage, DthIndicator) ~
    ADLCount + InterviewerRateHealth + Gender +
    LanguageAbility + CalculationAbility + waves,
  data    = cvdata_9800
)

# Survival tree partitioning
fit_stp <- mort_fit(
  model   = "stp",
  formula = Surv(birthage, censoredage, DthIndicator) ~
    ADLCount + InterviewerRateHealth + Gender +
    LanguageAbility + CalculationAbility + waves,
  data    = cvdata_9800,
  AC = "Gompertz"
)

# Mob tree partitioning
fit_mtp <- mort_fit(
  model   = "mtp",
  formula = Surv(birthage, censoredage, DthIndicator) ~ 1 |
    ADLCount + InterviewerRateHealth + Gender +
    LanguageAbility + CalculationAbility + waves,
  data    = cvdata_9800
)

# Gompertz proportional hazards (flexsurvreg)
fit_gph <- mort_fit(
  model   = "gph",
  formula = Surv(birthage, censoredage, DthIndicator) ~
    ADLCount + InterviewerRateHealth + Gender + LanguageAbility,
  data    = cvdata_9800
)
