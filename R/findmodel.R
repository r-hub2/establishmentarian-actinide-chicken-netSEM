
# Helper function for netSEMp1
## FF-1: Linear

netSEM_linear2 <- function(Resp, Var, x) {
  # Define significant digits
  sig_dig <- 3
  # Define formula for executing model
  formula <- glue::glue('{Resp} ~ {Var}')
  # Define a simple linear model
  model <- purrr::exec("lm", formula, data = x)

  # Reduce the size
  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-2: Quadratic with Linear Term

netSEM_quad2 <- function(Resp, Var, x) {
  # Define a general quadratic model
  sig_dig <- 3
  # Define a simple linear model
  formula <-
    glue::glue('{Resp} ~ {Var} + I({Var}^2)')
  model <- purrr::exec("lm", formula, data = x)

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-3: Quadratic with No Linear Term

netSEM_quadNoLinear2 <- function(Resp, Var, x) {
  # Define a general quadratic model
  sig_dig <- 3
  # Define a
  formula <-
    glue::glue('{Resp} ~ {Var}^2')
  model <- purrr::exec("lm", formula, data = x)

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-4: Exponential

netSEM_exponential2 <- function(Resp, Var, x) {
  sig_dig <- 3
  # Define a
  formula <-
    glue::glue('{Resp} ~ exp({Var})')
  #test for condition where Var data has negative or INF values which will cause lm(exp) model to fail
  Var_flag <- x[[enexpr(Var)]]
  model_flag <- all(exp(Var_flag[!is.na(Var_flag)]) != Inf)
  #if flag is TRUE, will return model statistics
  if(model_flag){
  # Define the exponential model
    model <- purrr::exec("lm", formula, data = x)
    model <- model %>%
      butcher::axe_data(verbose = FALSE) %>%
      butcher::axe_env(verbose = FALSE)
  } else {
    model <- NA
  }


  # Output the model
  return(model)
}

## FF-5: Logarithmic

netSEM_log2 <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ log({Var})')
  # Subset variable column for model validity check
  Var_flag <- x[[enexpr(Var)]]
  # Test for condition where Var data has negative or INF values which will cause lm(log) model to fail
  model_flag <- all(Var_flag > 0 & Var_flag < Inf, na.rm = TRUE)
  # If flag is TRUE, will return model statistics
  if(model_flag){
  # Define the log model
    model <- purrr::exec("lm", formula, data = x)
  } else {
    model <- NA
  }

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-6: Square Root (non-linear)

netSEM_SqRt2 <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ sqrt({Var})')
  # Define data to check for fail condition
  Var_flag <- x[[enexpr(Var)]]
  # Test for condition where Var data has negative or INF values which will cause lm(sqrt) model to fail
  model_flag <- all(as.logical(Var_flag > 0 &
                                 Var_flag < Inf, na.rm = TRUE))
  # If flag is TRUE, will return model statistics
  if(model_flag){
    # Define the SqRoot model
    model <- purrr::exec("lm", formula, data = x)
  } else {
    model <- NA
  }

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-7: Inverse Square Root (non-linear)

netSEM_InvSqRt2 <- function(Resp, Var, x) {
  # Define Significant Figures
  sig_dig <- 3
  # Define formula for generating model
  formula <-
    glue::glue('{Resp} ~ I(1/sqrt({Var}))')
  # Define data to check for fail condition
  Var_flag <- x[[enexpr(Var)]]
  # Test for condition where Var data has negative or INF values which will cause lm(sqrt) model to fail
  model_flag <- all(as.logical(Var_flag > 0 & Var_flag < Inf,
                               na.rm = TRUE))
  # If flag is TRUE, will return model statistics
  if(model_flag){
    # Define the Inverse SqRoot model
    model <- purrr::exec("lm", formula, data = x)
  } else {
    model <- NA
  }

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

## FF-8: NLS

netSEM_nls2 <- function(Resp, Var, x) {
  tryCatch({
  nls_formula <-
    glue::glue('{(Resp)} ~ a1 + a2 * exp(a3 * {(Var)})')

  resp <- x[[enexpr(Resp)]]
  var <- x[[enexpr(Var)]]

  #test for condition where Var data has negative or INF values which will cause nls model to fail
  Resp_flag <- x[[enexpr(Resp)]]
  model_flag <- all(as.logical(Resp_flag > 0 & Resp_flag < Inf,
                               na.rm = TRUE))
  #if flag is TRUE, will return model statistics
  if(model_flag){

  a1_0 <- min(resp) /2

  # Take the log of nls_formula
  nls_log_formula <-
    glue::glue('log(({(Resp)} - {a1_0})) ~ {(Var)}')
  model_0 <- purrr::exec("lm", nls_log_formula, data = x)

  # NLS model
  model <- stats::nls(
    formula = nls_formula,
    data = x,
    start = (tibble::lst(a1 = a1_0, a2 = exp(coef(model_0)[1]), a3 = coef(model_0)[2])), #arbitrary starting points
    control = stats::nls.control(maxiter = 50, minFactor = 1/2^20, warnOnly = TRUE) #either converges or stops at 1000 iterations
  )

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  } else {
    model <- NA
  }

  # Output the model
  return(model)

  },
  warning = function(w){
    message('Warning. Model',' Resp = ',Resp,', Var = ',Var,' not good fit for nls.')
    model_summary <- tibble(NA, model_type = "NLS")
  }
  )
}

## FF-9: Change Point

netSEM_CP2 <- function(Resp, Var, x) {

  sig_dig <- 3

  resp <- x[[{{Resp}}]]
  var <- x[[{{Var}}]]
  # formula <- glue::glue('{Resp} ~ {Var}')
  # seg.formula <- glue::glue('~ {var}')
  # model.obj <- purrr::exec("lm", formula, data = x)
  # Define the change point model
  # model <- purrr::exec(segmented.lm,
  #                      obj = model.obj,
  #                      seg.Z = ~ var,
  #                      psi = NA,
  #                      control = segmented::seg.control(
  #                        K = 1,
  #                        fix.npsi = FALSE,
  #                        n.boot = 0,
  #                        it.max = 20
  #                      ))

  # Define the Change Point model
  model <- segmented::segmented(
    obj = lm(resp ~ var, data = x),
    seg.Z = ~ var,
    psi = NA,
    control = segmented::seg.control(
      K = 1,
      fix.npsi = F,
      n.boot = 0,
      it.max = 20
    )
  )

  model <- model %>%
    butcher::axe_data(verbose = FALSE) %>%
    butcher::axe_env(verbose = FALSE)

  # Output the model
  return(model)
}

# Function to save each model
# list_netSEM_models <- function(Resp, Var, mod.name, x) {
#
#   list1 <- netSEM_linear2(Resp ,Var, x)
#   list2 <- netSEM_quad2(Resp, Var, x)
#   list3 <- netSEM_quadNoLinear2(Resp, Var, x)
#   list4 <- netSEM_exponential2(Resp, Var, x)
#   list5 <- netSEM_log2(Resp, Var, x)
#   list6 <- netSEM_SqRt2(Resp, Var, x)
#   list7 <- netSEM_InvSqRt2(Resp, Var, x)
#   list8 <- netSEM_nls2(Resp, Var, x)
#   list9 <- netSEM_CP2(Resp, Var, x)
#
#   list_final <-
#     list(Var, Resp, list1, list2, list3, list4, list5, list6, list7, list8, list9)
#
#   # Rename list names
#   names(list_final) <- c("exogenous", "endogenous", "linear", "quad", "quadNoLinear", "exponential",
#                          "log", "SqRt", "InvSqRt", "nls", "CP")
#
#
#   return(list_final)
# }

list_netSEM_models <- function(Resp, Var, mod.name, x) {
  if (ensym(mod.name) == "Linear") {
    list.mod <- netSEM_linear2(Resp , Var, x)
  }
  else if (ensym(mod.name)  == "Quad") {
    list.mod <- netSEM_quad2(Resp , Var, x)
  }
  else if (ensym(mod.name)  == "SQuad") {
    list.mod <- netSEM_quadNoLinear2(Resp, Var, x)
  }
  else if (ensym(mod.name)  == "Exponential") {
    list.mod <- netSEM_exponential2(Resp, Var, x)
  }
  else if (ensym(mod.name)  == "Log") {
    list.mod <- netSEM_log2(Resp, Var, x)
  }
  else if (ensym(mod.name)  == "SqRoot") {
    list.mod <- netSEM_SqRt2(Resp, Var, x)
  }
  else if (ensym(mod.name) == "InvSqRoot") {
    list.mod <- netSEM_InvSqRt2(Resp, Var, x)
  }
  else if (ensym(mod.name)  == "NLS") {
    list.mod <- netSEM_nls2(Resp, Var, x)
  }
  else if (ensym(mod.name)  == "CP") {
    list.mod <- netSEM_CP2(Resp, Var, x)
  }
  list_final <- list(Var, Resp, list.mod)
  names(list_final) <- c("exogenous", "endogenous",ensym(mod.name))
  return(list_final)
}