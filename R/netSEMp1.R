##' This function carries out netSEM using markovian principle.
##'
##' netSEMp1 determines the univariate relationships in the spirit of the Markovian process and builds a network model. In this case, the relationship between each pair of system variables, including predictors and the system level response is determined with the Markovian property that assumes the value of the current predictor is sufficient in relating to the next level variable, i.e., the relationship is independent of the specific value of the preceding-level variable to the current predictor, given the current value.
##' Each pair of variables is tested for sensible (in the domain knowledge sense) paring relation chosen from pre-selected common functional forms in linear regression settings.
##' Adjusted R-squared is used for model selection for every pair.
##'
##' P-values reported in the "res.print" field of the return list contains the P-values of estimators of linear regression coefficients.
##' The P-values are ordered in the common order of coefficients, i.e. in the order of increasing exponents.
##' For example, in the quadratic functional form \eqn{y \sim b_0 + b_1 x + b_2 x^2}, the three P-values correspond to those of \eqn{\hat{b}_0}, \eqn{\hat{b}_1} and \eqn{\hat{b}_2}, respectively.
##' If there are less than 3 coefficients to estimate, the extra P-value field is filled with NA's.
##'
##' @title network Structural Equation Modeling (netSEM)
##'
##' @param x A dataframe. By default it considers all columns as exogenous variables, except the first column which stores the system's endogenous variable.
##' @param exogenous, by default it considers all columns as exogenous variables except column number 1, which is the main endogenous response.
##' @param endogenous A character string of the column name of the main endogenous OR a numeric number indexing the column of the main endogenous.
##' @return An object of class netSEM, which is a list of the following items:
##'
##' \itemize{
##' \item "data": data frame used in netSEMp1 analysis.
##' \item "exogenous": Main stressor variable.
##' \item "endogenous": Main response variable.
##' \item "bestModels": A dataframe.
##' \item "allModels": A list of model objects.
##' }
##'
##'
##' @export
##'
##' @examples
##' \dontrun{
##' ## Load the sample acrylic data set
##' data(acrylic)
##'
##' ## Run netSEM_markovian
##' ans <- netSEMp1(acrylic, "IrradTot", "YI")
##' }

####### v0.7.0 ##################
netSEMp1 <- function(x, exogenous, endogenous) {
  # place holder for mfExt for 'no visible binding' note
  mfExt <- NULL
  # Check for exogenous & endogenous input. Error if missing.
  # Rearrange dataframe based on exogenous & endogenous specification.
  # Transformed dataframe structure is: Stressor, Response,...
  if (missing(exogenous) | missing(endogenous)) {
    stop("exogenous and endogenous variables need to be specified")
  } else {
    x <- x %>%
      dplyr::relocate(all_of(endogenous), all_of(exogenous))
  }
  x <<- x # save it in the global environment for pseudoR2

  # Create a vector for column names
  name.vec <- colnames(x)
  mfExt <<- x # save it in the global environment for segmented.lm
  # Create a response, variable data frame to map function
  df.names <- data.frame(expand.grid(var = name.vec, resp = name.vec)) %>%
    dplyr::filter(var != resp, # restrict self relationship
                  resp != {{exogenous}}, # restrict setting stressor as response
                  var != {{endogenous}}) %>% # restrict setting response as stressor
    dplyr::mutate_all(as.character)

  res_all_mod <-
    suppressWarnings(purrr::map2_dfr(
      .x = df.names$var,
      .y = df.names$resp,
      ~ collect_netSEM_models(Resp = .y, Var = .x, x = x)
    ))
  # remove unnecessary columns
  #"NA" and "term" columns don't always show up
  if ("NA" %in% names(res_all_mod)) {
    res_all_mod <- res_all_mod %>%
      dplyr::select(-c("NA"))
    if ("term" %in% names(res_all_mod)) {
      res_all_mod <- res_all_mod %>%
        dplyr::select(-c("term"))
    }
  } else {
    if ("term" %in% names(res_all_mod)) {
      res_all_mod <- res_all_mod %>%
        dplyr::select(-c("term"))
    }
  }

  # Find a best model for each response and variable combination
  best_model <- data.frame()
  for (i in 1:length(df.names$resp)) {
    temp <- (
      res_all_mod %>%
        dplyr::filter(Resp == df.names$resp[i], Var == df.names$var[i]) %>%
        dplyr::slice_max(adj.r.squared)
    )[1,]
    best_model <- dplyr::bind_rows(best_model, temp)
  }

  best_model <- best_model %>%
    dplyr::rename(best_model = model_type)

  # Save models in a list
  # model_list <-
  #   suppressWarnings(map2(
  #     .x = df.namesl$var,
  #     .y = df.nmaes$resp
  #     ~ list_netSEM_models(Resp = .y, Var = .x, x = x)
  #   ))

  # Save models in a list
  model_list <-
    suppressWarnings(tibble::tibble(Var = best_model$Var,
                            Resp = best_model$Resp,
                            mod = best_model$best_model) %>%
                       purrr::pmap(
      ~ list_netSEM_models(Resp = ..2, Var = ..1, mod = ..3, x = x)
    ))

  # Final output
  res <- list(
    data = x,
    exogenous = exogenous,
    endogenous = endogenous,
    bestModels = best_model,
    allModels = model_list
  )
  class(res) <- c("netSEMp1", "list")
  # remove "x" variable
  rm(x, envir = .GlobalEnv )
  # remove "mfExt" variable
  rm(mfExt, envir = .GlobalEnv )
  return(res)
}
