# =========================
# 1. Load Libraries & Setup
# =========================

# List of required libraries
libs <- c(
  "tidyverse", "grid", "libcoin", "mvtnorm", "partykit",
  "LTRCtrees", "labelled", "future.apply", "tictoc", "parallel", "doParallel",
  "foreach", "survival", "flexsurv", "DEoptim", "GenSA"
)
# Install any missing packages
installed_packages <- libs %in% rownames(installed.packages())
if (any(!installed_packages)) install.packages(libs[!installed_packages])
# Load all libraries
invisible(lapply(libs, require, character.only = TRUE))

# =========================
# 3. S4 Class Definitions
# =========================

setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))

# ---- 3.1 MortalityModel (Parent) ----
# S4 class for basic mortality model
setClass(
  "MortalityModel",
  slots = list(
    data = "list",
    data_train = "data.frameOrNULL",
    data_test = "data.frameOrNULL"
  )
)

# Constructor for MortalityModel
MortalityModel <- function(data = NULL, data_train = NULL, data_test = NULL) {
  if (is.data.frame(data)) {
    data_train <- data
    data_test <- NULL
    data <- list(train = data_train, test = data_test)
  }

  if (is.null(data)){
    data <- list(train = data_train, test = data_test)
  }

  new("MortalityModel",
      data = data,
      data_train = data$train,
      data_test = data$test
  )
}
# Show method for MortalityModel
setMethod("show", "MortalityModel", function(object) {
  if (is.null(object@data_train)) {
    cat("Empty Training data\n")
  } else cat("Training data size:", nrow(object@data_train), "\n")
  if (is.null(object@data_test)) {
    cat("Empty Test data\n")
  } else cat("Test data size:", nrow(object@data_test), "\n")
})

# Generics for IBS and logScore
setGeneric("getIBS", function(object) standardGeneric("getIBS"))
setGeneric("getlogScore", function(object) standardGeneric("getlogScore"))

# ---- 3.2 STModel (Survival Tree) ----
# S4 class for survival tree model
setOldClass("rpart")
setClass(
  "STModel",
  slots = list(
    IMM = "MortalityModel",
    model = "rpart",
    tnodes = "character"
  )
)
# Helper to fit survival tree
get_tree <- function(formula, data) {
  # data <- data %>% select(
  #   birthage, censoredage, DthIndicator,
  #   ADLCount, InterviewerRateHealth, Gender,
  #   LanguageAbility, CalculationAbility,
  #   waves
  # )
  data <- as.data.frame(data)
  # formula <- Surv(birthage, censoredage, DthIndicator) ~
  #   ADLCount + InterviewerRateHealth + Gender +
  #   LanguageAbility + CalculationAbility + waves
  LTRCART(
    formula,
    data = data
  )
}
# Constructor for STModel
STModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL, IMM=NULL) {
  if (is.null(IMM)) {
    # If IMM is not provided, create a new MortalityModel
    IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  }
  train_copy <- IMM@data_train
  test_copy <- IMM@data_test
  tree <- get_tree(formula, train_copy)
  train_copy$node <- as.vector(tree$where)
  tnodes <- train_copy$node %>% factor() %>% levels()
  test_copy$group <- predict(tree %>% as.party(), newdata = test_copy, type = "node")
  IMM@data <- list(train = train_copy, test = test_copy)
  IMM@data_train <- train_copy
  IMM@data_test <- test_copy
  new("STModel", IMM = IMM, model = tree, tnodes = tnodes)
}
# Show method for STModel
setMethod("show", "STModel", function(object) {
  cat("Training data size:", nrow(object@IMM@data_train), "\n")
  cat("Test data size:", nrow(object@IMM@data_test), "\n")
  cat("Terminal nodes:", paste(object@tnodes, collapse = ", "), "\n")
})

# ---- 3.3 BaselineModel ----
# S4 class for baseline parametric model
setClass(
  "BaselineModel",
  slots = list(
    IMM = "MortalityModel",
    model_fit = "PACFit"
  )
)
# Constructor for BaselineModel
BaselineModel <- function(data = NULL, data_train = NULL, data_test = NULL, AC=BeardPAC()) {
  IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  set.seed(5346615)
  pac_fit_obj <- new("PACFit", model = AC, data = IMM@data)
  pac_fit_obj <- fit(pac_fit_obj, IMM@data_train, patience = 500)
  new("BaselineModel", IMM = IMM, model_fit = pac_fit_obj)
}
# Show method for BaselineModel
setMethod("show", "BaselineModel", function(object) {
  cat("Training data size:", nrow(object@IMM@data_train), "\n")
  cat("Test data size:", nrow(object@IMM@data_test), "\n")
  cat("Model fit parameters:", object@model_fit@model@params, "\n")
})
# IBS and logScore methods for BaselineModel
setMethod("getIBS", "BaselineModel", function(object) {
  BS(object@model_fit, object@IMM@data_test, object@IMM@Surv_info)
})
setMethod("getlogScore", "BaselineModel", function(object) {
  sum(get_logllh(object@model_fit@model, object@IMM@data_test))
})

setMethod("predict", "BaselineModel", function(object, newdata=NULL, type = "survival") {
  if (type != "survival") {
    stop("Only 'survival' type is supported.")
  }
  if (is.null(newdata)) {
    newdata <- object@IMM@data_test
  }
  # Predict survival probabilities
  pac <- object@model_fit@model
  S <- pac %>% survival(newdata$censoredage)
  S_0 <- pac %>% survival(newdata$birthage)
  pred_probs <- S / S_0
  return(pred_probs)
})

# ---- 3.4 STPModel (Survival Tree Partitioning) ----
# S4 class for partitioned tree model
setClass(
  "STPModel",
  slots = list(
    STModel = "STModel",
    AC = "PAC",
    par_train = "matrix"
  )
)
# Constructor for STPModel
STPModel <- function(data = NULL, data_train = NULL, data_test = NULL, AC = BeardPAC(), ST = NULL) {
  if (is.null(ST)) {
    # If ST is not provided, create a new STModel
    ST <- STModel(data = data, data_train = data_train, data_test = data_test)
  }
  par_train <- c()
  for (nodei in ST@tnodes) {
    set.seed(5346615)
    pac_fit_obj <- new("PACFit", model = AC, data = ST@IMM@data)
    pac_fit_obj = fit(pac_fit_obj, ST@IMM@data_train %>% filter(node == nodei), patience = 500)
    par_train <- par_train %>% rbind(pac_fit_obj@model@params)
  }
  new("STPModel", STModel = ST, AC = AC, par_train = par_train)
}
# Show method for STPModel
setMethod("show", "STPModel", function(object) {
  cat("Training data size:", nrow(object@STModel@IMM@data_train), "\n")
  cat("Test data size:", nrow(object@STModel@IMM@data_test), "\n")
  cat("Terminal nodes:", paste(object@STModel@tnodes, collapse = ", "), "\n")
})
# IBS and logScore methods for STPModel
setMethod("getIBS", "STPModel", function(object) {
  tnodes = object@STModel@tnodes
  test_copy = object@STModel@IMM@data_test
  sum_count = 0
  for (node_num in seq_along(tnodes)) {
    sub_test_copy <- test_copy %>% filter(test_copy$group == tnodes[node_num])
    PACfit <- new("PACFit", model = object@AC, data = object@STModel@IMM@data)
    PACfit@model@params <- object@par_train[node_num, ]
    PACfit@fitted <- TRUE
    BS_node <- BS(PACfit, sub_test_copy, object@STModel@IMM@Surv_info)
    sum_count = sum_count + nrow(sub_test_copy) * BS_node
  }
  return(sum_count / nrow(test_copy))
})
setMethod("getlogScore", "STPModel", function(object) {
  logllh <- 0
  tnodes = object@STModel@tnodes
  test_copy = object@STModel@IMM@data_test
  for (node_num in seq_along(tnodes)) {
    sub_test_copy <- test_copy %>% filter(test_copy$group == tnodes[node_num])
    PACfit <- new("PACFit", model = object@AC, data = object@STModel@IMM@data)
    PACfit@model@params <- object@par_train[node_num, ]
    PACfit@fitted <- TRUE
    logllh_node <- sum(get_logllh(PACfit@model, sub_test_copy))
    logllh <- logllh + logllh_node
  }
  return(logllh)
})
# setMethod("predict", "STPModel", function(object, newdata = NULL, type = "survival") {
#   if (type != "survival") {
#     stop("Only 'survival' type is supported.")
#   }
#   if (is.null(newdata)) {
#     newdata <- object@STModel@IMM@data_test
#   }
#   else{
#     newdata$group <- predict(tree %>% as.party(), newdata = newdata, type = "node")
#   }
#   # Predict survival probabilities for each node
#   df_S <- data.frame()
#   for (nodei in object@STModel@tnodes) {
#     pac <- object@AC
#     pac@params <- object@par_train[match(nodei, object@STModel@tnodes), ]
#     S <- pac %>% survival(newdata$censoredage)
#     S_0 <- pac %>% survival(newdata$birthage)
#     df_S <- df_S %>% rbind(S/S_0)
#   }
#   df_S <- t(df_S)
#   colnames(df_S) <- object@STModel@tnodes
#   # get individual survival given predicted group using vectorized approach
#   group_idx <- match(newdata$group, object@STModel@tnodes)
#   pred_probs <- df_S[cbind(seq_len(nrow(newdata)), group_idx)]
#
#   return(pred_probs)
# })

# ---- 3.5 MortalityModel_FillNA ----
# S4 class for model with missing values filled
setClass(
  "MortalityModel_FillNA",
  slots = list(
    IMM = "MortalityModel",
    train_fillna = "data.frame",
    test_fillna = "data.frame"
  )
)
# Constructor for MortalityModel_FillNA
MortalityModel_FillNA <- function(data = NULL, data_train = NULL, data_test = NULL) {
  IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  train_fillna <- IMM@data_train
  # Fill NAs with "Missing" for categorical variables
  for (i in 5:ncol(train_fillna)) {
    train_fillna[[i]] <- train_fillna[[i]] %>%
      as.character() %>% replace_na("Missing") %>% as.factor()
  }
  test_fillna <- IMM@data_test
  for (i in 5:ncol(test_fillna)) {
    test_fillna[[i]] <- test_fillna[[i]] %>%
      as.character() %>% replace_na("Missing") %>% as.factor()
  }
  new("MortalityModel_FillNA",
      IMM = IMM,
      train_fillna = train_fillna,
      test_fillna = test_fillna
  )
}
# Show method for MortalityModel_FillNA
setMethod("show", "MortalityModel_FillNA", function(object) {
  cat("Training data size:", nrow(object@train_fillna), "\n")
  cat("Test data size:", nrow(object@test_fillna), "\n")
})

# ---- 3.6 MTPModel ----
# S4 class for model-based recursive partitioning (mob)
setOldClass("modelparty")
setClass(
  "MTPModel",
  slots = list(
    IMM_filled = "MortalityModel_FillNA",
    formula = "formula",
    mobfit = "modelparty",
    mob_coef = "matrix",
    mob_tnodes_test = "numeric"
  )
)
# Constructor for MTPModel
MTPModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL) {
  IMM_filled <- MortalityModel_FillNA(data = data, data_train = data_train, data_test = data_test)
  # Define partitioning formula
  # mobformula <- Surv(birthage, censoredage, DthIndicator) ~ 1 |
  #   ADLCount + InterviewerRateHealth + Gender +
  #   LanguageAbility + Illness + GeneralAbility + Education +
  #   SmokeEver + CalculationAbility +
  #   ReactionAbility + MemorisingAbility + DrinkEver +
  #   incomesource + SelfratedHealth
  # if (data_name == "doverall") {
  #   mobformula <- Surv(birthage, censoredage, DthIndicator) ~ 1 |
  #     ADLCount + InterviewerRateHealth + Gender +
  #     LanguageAbility + Illness + GeneralAbility + Education +
  #     SmokeEver + CalculationAbility +
  #     ReactionAbility + MemorisingAbility + DrinkEver +
  #     incomesource + SelfratedHealth + waves
  # }
  # Custom fit function for mob
  gpreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
    flexsurvreg(y ~ 0 + x, weights = weights, dist = "gompertz", ...)
  }
  # Fit mob model
  mobfit <- mob(
    formula,
    data = IMM_filled@train_fillna,
    fit = gpreg,
    control = mob_control(minsize = nrow(IMM_filled@train_fillna) / 50, maxdepth = 10)
  )
  mob_coef <- coef(mobfit)
  mob_tnodes_test <- predict(mobfit, IMM_filled@test_fillna, type = "node")
  # If only one node, fallback to baseline
  if(length(unique(mob_tnodes_test))==1){
    IMM <- IMM_filled@IMM
    GompertzFit <- new("PACFit", model = GompertzPAC(), data = IMM@data)
    GompertzFit@model@params <- c(mob_coef[2], mob_coef[1])
    GompertzFit@fitted <- TRUE
    new("BaselineModel", IMM = IMM, model_fit = GompertzFit)
  } else {
    new("MTPModel",
        IMM_filled = IMM_filled,
        formula = formula,
        mobfit = mobfit,
        mob_coef = mob_coef,
        mob_tnodes_test = mob_tnodes_test
    )
  }
}
# Show method for MTPModel
setMethod("show", "MTPModel", function(object) {
  cat("Training data size:", nrow(object@IMM_filled@train_fillna), "\n")
  cat("Test data size:", nrow(object@IMM_filled@test_fillna), "\n")
  cat("Terminal nodes in MTP model:", nrow(object@mob_coef), "\n")
})
# IBS and logScore methods for MTPModel
setMethod("getIBS", "MTPModel", function(object) {
  tnodes = object@mob_coef %>% rownames()
  test_copy <- object@IMM_filled@test_fillna
  sum_count <- 0
  for (i in seq_along(tnodes)) {
    sub_test_copy <- test_copy %>% filter(object@mob_tnodes_test == tnodes[i])
    pari <- object@mob_coef[match(tnodes[i], rownames(object@mob_coef)), ]
    pari <- pari[2:1]
    PACfit <- new("PACFit", model = GompertzPAC(), data = object@IMM_filled@IMM@data)
    PACfit@model@params <- pari
    PACfit@fitted <- TRUE
    BS_node <- BS(PACfit, sub_test_copy, object@IMM_filled@IMM@Surv_info)
    sum_count <- sum_count + nrow(sub_test_copy) * BS_node
  }
  ibs_mob <- sum_count / nrow(test_copy)
  return(ibs_mob)
})
setMethod("getlogScore","MTPModel", function(object) {
  logllh <- 0
  tnodes = object@mob_coef %>% rownames()
  test_copy <- object@IMM_filled@test_fillna
  for (i in seq_along(tnodes)) {
    sub_test_copy <- test_copy %>% filter(object@mob_tnodes_test == tnodes[i])
    pari <- object@mob_coef[match(tnodes[i], rownames(object@mob_coef)), ]
    pari <- pari[2:1]
    PACfit <- new("PACFit", model = GompertzPAC(), data = object@IMM_filled@IMM@data)
    PACfit@model@params <- pari
    PACfit@fitted <- TRUE
    logllh_node <- sum(get_logllh(PACfit@model, sub_test_copy))
    logllh <- logllh + logllh_node
  }
  return(logllh)
})
# setMethod("predict", "MTPModel", function(object, newdata = NULL, type = "survival") {
#   if (type != "survival") {
#     stop("Only 'survival' type is supported.")
#   }
#   if (is.null(newdata)) {
#     newdata <- object@IMM_filled@test_fillna
#     newdata$group <- object@mob_tnodes_test
#   }
#   else{
#     newdata$group <- predict(mobfit, newdata, type = "node")
#   }
#   tnodes = object@mob_coef %>% rownames()
#   mob_coef <- object@mob_coef[,2:1] # switch order to match Gompertz parameters
#   # Predict survival probabilities for each node
#   df_S <- data.frame()
#   for (nodei in tnodes) {
#     pac <- GompertzPAC()
#     pac@params <- mob_coef[match(nodei, tnodes), ]
#     S <- pac %>% survival(newdata$censoredage)
#     S_0 <- pac %>% survival(newdata$birthage)
#     df_S <- df_S %>% rbind(S/S_0)
#   }
#   df_S <- t(df_S)
#   colnames(df_S) <- tnodes
#   # get individual survival given predicted group using vectorized approach
#   group_idx <- match(newdata$group, tnodes)
#   pred_probs <- df_S[cbind(seq_len(nrow(newdata)), group_idx)]
#
#   return(pred_probs)
# })

# ---- 3.7 GPHModel ----
# S4 class for Gompertz proportional hazards model
setOldClass("flexsurvreg")
setClass(
  "GPHModel",
  slots = list(
    IMM_filled = "MortalityModel_FillNA",
    formula = "formula",
    gpfit = "flexsurvreg",
    gpfit_pred = "numeric",
    gpfit_shape = "numeric"
  )
)
# Constructor for GPHModel
GPHModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL) {
  IMM_filled <- MortalityModel_FillNA(data = data, data_train = data_train, data_test = data_test)
  set.seed(5346615)
  # gpformula <- Surv(birthage, censoredage, DthIndicator) ~
  #   ADLCount + InterviewerRateHealth + Gender +
  #   LanguageAbility + Illness + GeneralAbility + Education +
  #   SmokeEver + CalculationAbility +
  #   ReactionAbility + MemorisingAbility + DrinkEver +
  #   incomesource + SelfratedHealth
  # if (data_name == "doverall") {
  #   gpformula <- update(gpformula, . ~ . + waves)
  # }
  gpfit <- flexsurvreg(formula, data = IMM_filled@train_fillna, dist = "gompertz")
  gpfit_pred <- gpfit %>% predict(IMM_filled@test_fillna, type = "lp")
  gpfit_pred <- gpfit_pred$.pred_link
  gpfit_pred <- log(gpfit_pred)
  gpfit_shape <- gpfit$coefficients[1]
  new("GPHModel",
      IMM_filled = IMM_filled,
      formula = formula,
      gpfit = gpfit,
      gpfit_pred = gpfit_pred,
      gpfit_shape = gpfit_shape
  )
}
# Show method for GPHModel
setMethod("show", "GPHModel", function(object) {
  cat("Training data size:", nrow(object@IMM_filled@train_fillna), "\n")
  cat("Test data size:", nrow(object@IMM_filled@test_fillna), "\n")
  cat("Shape parameter:", object@gpfit_shape, "\n")
})
# IBS and logScore methods for GPHModel
setMethod("getIBS", "GPHModel", function(object) {
  test_fillna <- object@IMM_filled@test_fillna
  Surv_info <- object@IMM_filled@IMM@Surv_info
  ibs_recorder <- rep(0, nrow(test_fillna))
  PACfit <- new("PACFit", model = GompertzPAC(), data = object@IMM_filled@IMM@data)
  PACfit@fitted <- TRUE
  for (i in seq_len(nrow(test_fillna))) {
    PACfit@model@params <- c(object@gpfit_pred[i], object@gpfit_shape)
    ibs_recorder[i] <- BS(PACfit, test_fillna[i, ], Surv_info = Surv_info)
  }
  ibs_gph <- mean(ibs_recorder)
  return(ibs_gph)
})
setMethod("getlogScore", "GPHModel", function(object) {
  logllh = 0
  test_fillna <- object@IMM_filled@test_fillna
  Surv_info <- object@IMM_filled@IMM@Surv_info
  PACfit <- new("PACFit", model = GompertzPAC(), data = object@IMM_filled@IMM@data)
  PACfit@fitted <- TRUE
  for (i in seq_len(nrow(test_fillna))) {
    PACfit@model@params <- c(object@gpfit_pred[i], object@gpfit_shape)
    logllh = logllh + get_logllh(PACfit@model, test_fillna[i,])
  }
  return(logllh)
})
setMethod("predict", "GPHModel", function(object, newdata = NULL, type = "survival") {
  if (type != "survival") {
    stop("Only 'survival' type is supported.")
  }
  if (is.null(newdata)) {
    newdata <- object@IMM_filled@test_fillna
  }
  # Predict survival probabilities
  getS <- function(alpha, ca, ba){
    pac <- GompertzPAC()
    pac@params <- c(alpha, object@gpfit_shape)
    S <- pac %>% survival(ca)
    S_0 <- pac %>% survival(ba)
    return(S / S_0)
  }

  pred_probs <- mapply(getS, object@gpfit_pred, newdata$censoredage, newdata$birthage)

  return(pred_probs)
})

# =========================
# 4. Utility Functions
# =========================

# Log-likelihood extractor for flexsurvreg
logLik.flexsurvreg <- function(object, ...) {
  structure(object$loglik, df = object$npars, class = "logLik")
}
# Empirical estimating function for flexsurvreg (Gompertz only)
estfun.flexsurvreg <- function(object, ...) {
  if (object$dlist$name == "gompertz") {
    Y <- object$data$Y[, 1:3]
    Li <- Y[, 1]
    Ci <- Y[, 2]
    di <- Y[, 3]
    a <- object$coefficients[1]
    log_b <- object$coefficients[2]
    b <- exp(log_b)
    est_shape <- di * Ci - b / (a^2) * (exp(a * Ci) * (a * Ci - 1) - exp(a * Li) * (a * Li - 1))
    est_rate <- di / b - (exp(a * Ci) - exp(a * Li)) / a
    return(cbind(est_shape, est_rate))
  }
}

# =========================
# Unified formula-first fitters
# =========================

# Helper: map AC family input -> PAC object
.AC_registry <- list(
  Gompertz  = GompertzPAC,
  Makeham   = MakehamPAC,
  Perks     = PerksPAC,
  Beard     = BeardPAC,
  Kannisto  = KannistoPAC,
  Thatcher  = ThatcherPAC
)
AC_from <- function(AC = "Beard") {
  if (inherits(AC, "PAC")) return(AC)
  if (is.character(AC)) {
    ctor <- .AC_registry[[AC]]
    if (is.null(ctor)) stop("Unknown AC family: ", AC)
    return(ctor())
  }
  if (is.function(AC)) return(AC())
  stop("AC must be a PAC object, a family name, or a family constructor.")
}

# Small guard: baseline requires ~ 1
.ensure_age_only <- function(formula) {
  tl <- attr(stats::terms(formula), "term.labels")
  if (length(tl) > 0) {
    warning("Baseline model are age-only; RHS covariates will be ignored.")
  }
}

# 1) Baseline age-only parametric model (PAC)
baseline_fit <- function(formula,
                         data = NULL,
                         data_train = NULL,
                         data_test = NULL,
                         AC = "Beard",
                         patience = 500, tol = 1e-6, verbose = FALSE) {
  .ensure_age_only(formula)
  # Reuse your formula-based PAC fitter
  # (expects Surv(birthage, censoredage, DthIndicator) ~ 1)
  pf <- pac_fit(formula, data = data_train,
                AC = if (inherits(AC, "PAC")) AC@name else AC,
                patience = patience, tol = tol, verbose = verbose)
  IMM <- MortalityModel(data = data)
  methods::new("BaselineModel", IMM = IMM, model_fit = pf)
}

# 2) Survival Tree (LTRCART)
st_fit <- function(formula,
                   data = NULL,
                   data_train = NULL,
                   data_test = NULL) {
  STModel(formula = formula, data = data)
}

# 3) Survival Tree Partitioning (tree+PAC per node)
stp_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL,
                    AC = "Beard") {
  st <- STModel(formula = formula, data = data)
  STPModel(ST = st, AC = AC_from(AC))
}

# 4) Model-based recursive partitioning (mob + Gompertz)
mtp_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL) {
  MTPModel(formula = formula, data = data)
}

# 5) Gompertz proportional hazards (flexsurvreg)
gph_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL) {
  GPHModel(formula = formula, data = data)
}

# 6) One front door for everything
mort_fit <- function(model = c("baseline","st","stp","mtp","gph"),
                     formula,
                     data = NULL,
                     data_train = NULL,
                     data_test = NULL,
                     ...) {
  if (is.data.frame(data)) {
    data_train <- data
    data_test <- NULL
    data <- list(train = data_train, test = data_test)
  }

  if (is.null(data)){
    data <- list(train = data_train, test = data_test)
  }

  model <- match.arg(model)
  switch(model,
         baseline = baseline_fit(formula, data, data_train, data_test, ...),
         st       = st_fit(formula, data, data_train, data_test),
         stp      = stp_fit(formula, data, data_train, data_test, ...),
         mtp      = mtp_fit(formula, data, data_train, data_test),
         gph      = gph_fit(formula, data, data_train, data_test),
         stop("Unknown model type.")
  )
}


# ============================================================
# Shared helpers
# ============================================================

# Predict-based IBS from a precomputed S_pred matrix (n x |ages|)
.ibs_from_Spred <- function(newdata, S_pred, ages, normalize = TRUE) {
  n <- nrow(newdata); m <- length(ages)
  CA  <- matrix(as.numeric(newdata$censoredage),  n, m)
  DI  <- matrix(as.integer(newdata$DthIndicator), n, m)
  AGE <- matrix(rep(ages, each = n), n, m)

  # S_emp: 0 if x <= CA & death; 1 if x < CA; NA otherwise
  S_emp <- ifelse(AGE <= CA & DI == 1L, 0,
                  ifelse(AGE < CA, 1, NA_real_))
  # weights: 0 if x > CA & censored; 1 otherwise
  w <- ifelse(AGE > CA & DI == 0L, 0, 1)

  diff2 <- (S_pred - S_emp)^2
  diff2[!is.finite(diff2)] <- NA_real_

  contrib <- diff2 * w
  denom_total <- sum(w, na.rm = TRUE)
  IBS_sum     <- sum(contrib, na.rm = TRUE)
  IBS_mean    <- if (normalize) IBS_sum / max(1, denom_total) else IBS_sum

  # Age-weighted IBS_by_age (sums to IBS when normalized)
  denom_age   <- colSums(w, na.rm = TRUE)
  by_age_mean <- colSums(contrib, na.rm = TRUE) / pmax(1, denom_age)
  if (normalize) {
    weights_age <- denom_age / max(1, denom_total)
    IBS_by_age  <- by_age_mean * weights_age
  } else {
    IBS_by_age  <- by_age_mean * denom_age
  }
  names(IBS_by_age) <- as.character(ages)

  list(IBS = IBS_mean, IBS_by_age = IBS_by_age)
}

# ============================================================
# Common evaluate() generic
# ============================================================
setGeneric("evaluate", function(object, newdata, ...) standardGeneric("evaluate"))

# ============================================================
# BaselineModel: predict/evaluate (delegate to contained PACFit)
# ============================================================

setMethod("predict", "BaselineModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik"),
                   age_range = 80:115) {
            type <- match.arg(type)
            if (is.null(newdata)) newdata <- object@IMM@data_test
            predict.PACFit(object@model_fit, newdata = newdata, type = type, age_range = age_range)
          })

setMethod("evaluate", "BaselineModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   ibs_method = c("predict","ipcw"),
                   age_range = 80:115,
                   normalize_IBS = TRUE,
                   ...) {
            if (is.null(newdata)) newdata <- object@IMM@data_test
            ibs_method <- match.arg(ibs_method)
            # reuse evaluate.PACFit
            evaluate.PACFit(object@model_fit, newdata = newdata, metrics = metrics,
                            ibs_method = ibs_method, age_range = age_range,
                            normalize_IBS = normalize_IBS, ...)
          })

# ============================================================
# STPModel: predict/evaluate (per-node PAC parameters)
# ============================================================

setMethod("predict", "STPModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "group", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            if (is.null(newdata)) {
              newdata <- object@STModel@IMM@data_test
              if (is.null(newdata$group)) {
                newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                         newdata = newdata, type = "node")
              }
            } else {
              newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                       newdata = newdata, type = "node")
            }

            tnodes <- object@STModel@tnodes
            # Precompute per-node functions at 'ages'
            if (type %in% c("survival","hazard")) {
              S_or_mu <- matrix(NA_real_, nrow(newdata), length(ages))
              for (nodei in tnodes) {
                idx <- which(newdata$group == nodei)
                if (!length(idx)) next
                pac <- object@AC
                pac@params <- object@par_train[match(nodei, tnodes), ]
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)             # |ages|
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[idx])
                  M     <- matrix(rep(S_age, each = length(idx)), nrow = length(idx))
                  M     <- sweep(M, 1, S_ba, "/")
                  # ages < birthage -> 1
                  for (k in seq_along(idx)) {
                    younger <- which(ages < newdata$birthage[idx[k]])
                    if (length(younger)) M[k, younger] <- 1
                  }
                  S_or_mu[idx, ] <- M
                } else {
                  mu_age <- pac@mu_(pac@params, ages)
                  S_or_mu[idx, ] <- matrix(rep(mu_age, each = length(idx)), nrow = length(idx))
                }
              }
              colnames(S_or_mu) <- as.character(ages)
              return(S_or_mu)
            }

            # type == "loglik"
            # Sum per-row log-likelihood using node-specific PAC
            ll <- rep(NA_real_, nrow(newdata))
            for (nodei in tnodes) {
              idx <- which(newdata$group == nodei)
              if (!length(idx)) next
              pac <- object@AC
              pac@params <- object@par_train[match(nodei, tnodes), ]
              ll[idx] <- logllh(pac,
                                par = pac@params,
                                x = newdata$birthage[idx],
                                r = newdata$censoredage[idx],
                                i = newdata$DthIndicator[idx])
            }

            if (type == "group") return( newdata$group)
            if (type == "par")   return( object@par_train[match(newdata$group, tnodes), , drop = FALSE])

            names(ll) <- rownames(newdata)
            ll
          })



setMethod("evaluate", "STPModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) {
              newdata <- object@STModel@IMM@data_test
              if (is.null(newdata$group)) {
                newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                         newdata = newdata, type = "node")
              }
            }
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC not defined for STP (no single parametric likelihood)
            if ("AIC" %in% metrics) out$AIC <- NA_real_

            # loglik
            if ("loglik_vec" %in% metrics || "loglik_sum" %in% metrics) {
              ll_vec <- predict(object, newdata, type = "loglik")
              if ("loglik_vec" %in% metrics) out$loglik_vec <- ll_vec
              if ("loglik_sum" %in% metrics) out$loglik_sum <- sum(ll_vec, na.rm = TRUE)
            }

            # IBS via predict-based method

            if ("IBS" %in% metrics) {
              S_pred <- predict(object, newdata, type = "survival", age_range = ages)  # n x |ages|
              ibs <- .ibs_from_Spred(newdata, S_pred, ages, normalize = normalize_IBS)
              out$IBS <- ibs$IBS
              out$IBS_by_age <- ibs$IBS_by_age
            }

            out
          })

# ============================================================
# MTPModel: predict/evaluate (mob nodes -> Gompertz PAC per node)
# ============================================================

setMethod("predict", "MTPModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "group", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            tnodes <- rownames(object@mob_coef)
            if (is.null(newdata)) {
              newdata <- object@IMM_filled@test_fillna
              newdata$group <- object@mob_tnodes_test
            } else {
              newdata$group <- predict(object@mobfit, newdata, type = "node")
            }
            mob_coef <- object@mob_coef[, 2:1]  # reorder to (a, b)

            if (type %in% c("survival","hazard")) {
              M <- matrix(NA_real_, nrow(newdata), length(ages))
              for (nodei in tnodes) {
                idx <- which(newdata$group == nodei)
                if (!length(idx)) next
                pac <- GompertzPAC()
                pac@params <- mob_coef[match(nodei, tnodes), ]
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[idx])
                  B     <- matrix(rep(S_age, each = length(idx)), nrow = length(idx))
                  B     <- sweep(B, 1, S_ba, "/")
                  for (k in seq_along(idx)) {
                    younger <- which(ages < newdata$birthage[idx[k]])
                    if (length(younger)) B[k, younger] <- 1
                  }
                  M[idx, ] <- B
                } else {
                  mu_age <- pac@mu_(pac@params, ages)
                  M[idx, ] <- matrix(rep(mu_age, each = length(idx)), nrow = length(idx))
                }
              }
              colnames(M) <- as.character(ages)
              return(M)
            }

            # type == "loglik"
            ll <- rep(NA_real_, nrow(newdata))
            for (nodei in tnodes) {
              idx <- which(newdata$group == nodei)
              if (!length(idx)) next
              pac <- GompertzPAC()
              pac@params <- mob_coef[match(nodei, tnodes), ]
              ll[idx] <- logllh(pac,
                                par = pac@params,
                                x = newdata$birthage[idx],
                                r = newdata$censoredage[idx],
                                i = newdata$DthIndicator[idx])
            }

            if (type == "group") return( newdata$group)
            if (type == "par")   return( mob_coef[match(newdata$group, tnodes), , drop = FALSE])

            names(ll) <- rownames(newdata)
            ll
          })

setMethod("evaluate", "MTPModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) {
              newdata <- object@IMM_filled@test_fillna
              newdata$group <- object@mob_tnodes_test
            }
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC not defined for the piecewise-node model; return NA
            if ("AIC" %in% metrics) out$AIC <- NA_real_

            if ("loglik_vec" %in% metrics || "loglik_sum" %in% metrics) {
              ll_vec <- predict(object, newdata, type = "loglik")
              if ("loglik_vec" %in% metrics) out$loglik_vec <- ll_vec
              if ("loglik_sum" %in% metrics) out$loglik_sum <- sum(ll_vec, na.rm = TRUE)
            }

            if ("IBS" %in% metrics) {
              S_pred <- predict(object, newdata, type = "survival", age_range = ages)
              ibs <- .ibs_from_Spred(newdata, S_pred, ages, normalize = normalize_IBS)
              out$IBS <- ibs$IBS
              out$IBS_by_age <- ibs$IBS_by_age
            }
            out
          })

# ============================================================
# GPHModel: predict/evaluate (per-individual Gompertz PH)
# ============================================================

setMethod("predict", "GPHModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            if (is.null(newdata)) newdata <- object@IMM_filled@test_fillna

            # shape (b) common; alpha_i varies via gpfit_pred (stored on object)
            b <- as.numeric(object@gpfit_shape)
            # If predicting on newdata different from test_fillna, recompute lp quickly:
            if (!identical(newdata, object@IMM_filled@test_fillna)) {
              # best-effort: use flexsurv predict with type = "lp"
              lp <- predict(object@gpfit, newdata = newdata, type = "lp")$.pred_link
              alpha_vec <- log(lp)
            } else {
              alpha_vec <- object@gpfit_pred
            }

            if (type %in% c("survival","hazard")) {
              M <- matrix(NA_real_, nrow(newdata), length(ages))
              for (i in seq_len(nrow(newdata))) {
                pac <- GompertzPAC()
                pac@params <- c(alpha_vec[i], b)
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[i])
                  v     <- S_age / S_ba
                  v[ages < newdata$birthage[i]] <- 1
                  M[i, ] <- v
                } else {
                  M[i, ] <- pac@mu_(pac@params, ages)
                }
              }
              colnames(M) <- as.character(ages)
              return(M)
            }

            # type == "loglik" (per-individual)
            ll <- rep(NA_real_, nrow(newdata))
            for (i in seq_len(nrow(newdata))) {
              pac <- GompertzPAC()
              pac@params <- c(alpha_vec[i], b)
              ll[i] <- logllh(pac,
                              par = pac@params,
                              x = newdata$birthage[i],
                              r = newdata$censoredage[i],
                              i = newdata$DthIndicator[i])
            }

            if (type == "par") return( cbind(alpha = alpha_vec, shape = rep(b, nrow(newdata))) )

            names(ll) <- rownames(newdata)
            ll
          })

setMethod("evaluate", "GPHModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) newdata <- object@IMM_filled@test_fillna
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC available via flexsurvreg
            if ("AIC" %in% metrics) out$AIC <- tryCatch(stats::AIC(object@gpfit), error = function(e) NA_real_)

            if ("loglik_vec" %in% metrics || "loglik_sum" %in% metrics) {
              ll_vec <- predict(object, newdata, type = "loglik")
              if ("loglik_vec" %in% metrics) out$loglik_vec <- ll_vec
              if ("loglik_sum" %in% metrics) out$loglik_sum <- sum(ll_vec, na.rm = TRUE)
            }


            if ("IBS" %in% metrics) {
              S_pred <- predict(object, newdata, type = "survival", age_range = ages)
              ibs <- .ibs_from_Spred(newdata, S_pred, ages, normalize = normalize_IBS)
              out$IBS <- ibs$IBS
              out$IBS_by_age <- ibs$IBS_by_age
            }

            out
          })

# ============================================================
# (Optional) STModel: node predictions only
# ============================================================
setMethod("predict", "STModel",
          function(object, newdata = NULL,
                   type = c("survival", "node"),
                   age_range = 80:115) {
            type <- match.arg(type)
            if (is.null(newdata)) newdata <- object@IMM@data_test
            if (type == "survival") {

            }
            if (type == "node")
              return(predict(object@model %>% partykit::as.party(), newdata = newdata, type = "node"))
          })

setMethod("evaluate", "STModel",
          function(object, newdata = NULL, ...) {
            if (is.null(newdata)) newdata <- object@IMM@data_test
            preds <- predict(object, newdata, type = "node")
            list(nodes = object@tnodes,
                 n_by_node = as.integer(table(factor(preds, levels = object@tnodes))))
          })
