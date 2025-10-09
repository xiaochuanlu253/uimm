# =========================
# MAIN EXECUTION: PAC.r
# =========================

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# load dataset
load("../data/cvdata.RData")
# load code from age only model implement.r
source("./PAC.r")

train14 = cvdata_1418$train
pf_Gompertz= pac_fit(Surv(birthage, censoredage, DthIndicator) ~ 1,
                     data = train14,
                     AC = "Gompertz")

predict(pf_Gompertz, newdata = cvdata_1418$test)
evaluate(pf_Gompertz, newdata = cvdata_1418$test)

load("../data/data_sb_nrc.RData")
nhis = data_sb_nrc$df_complete

km_fit = survfit(Surv(precise_age, Dth_age, DthIndicator) ~ 1, data = nhis)
pf_Gompertz_nhis = pac_fit(Surv(precise_age, Dth_age, DthIndicator) ~ 1,
                             data = nhis,
                             AC = "Gompertz")

s_vec = survival_gompertz(pf_Gompertz_nhis@fit_results$par, x = seq(0, 100, by = 0.1))
df_s_vec = data.frame(time = seq(0, 100, by = 0.1), surv = s_vec)

ggplot() +
  geom_step(aes(x = km_fit$time, y = km_fit$surv), color = "blue") +
  geom_line(data = df_s_vec, aes(x = time, y = surv), color = "red") +
  labs(title = "Kaplan-Meier vs Gompertz Survival Curve",
       x = "Age",
       y = "Survival Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 100))