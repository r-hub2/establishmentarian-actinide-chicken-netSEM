##' This function selects the variable that best represents a specific degradation mechanism with the greatest adjusted R^2 of <S|M|.
##'
##' @title network Structural Equation Modeling (netSEM)
##'
##' @param a A dataframe. By default it considers all columns as exogenous variables, except the first column which stores the system's endogenous variable.
##' @param exogenous, by default it considers all columns as exogenous variables except column number 1, which is the main endogenous response.
##' @param endogenous A character string of the column name of the main endogenous OR a numeric number indexing the column of the main endogenous.
##' @return A character string of the variable that best represents a specific degradation mechanism.
##'
##' @export
##'
##' @examples
##' \dontrun{
##' ## Load the sample acrylic data set
##' data(acrylic)
##'
##' ## Run netSEMp1_predict
##' ans <- var_selection(a = acrylic, exogenous = "IrradTot", endogenous = "YI")
##' }

####### v0.8.0 ##################
var_selection <- function(a, exogenous, endogenous) {

  # Check for exogenous & endogenous input. Error if missing.
  # Rearrange dataframe based on exogenous & endogenous specification.
  # Transformed dataframe structure is: Stressor, Response,...
  if (missing(exogenous) | missing(endogenous)) {
    stop("exogenous and endogenous variables need to be specified")
  } else {
    a <- a %>%
      dplyr::relocate(all_of(endogenous), all_of(exogenous))
  }

  # Run netSEMp1 to get pairwise relationships
  p1.result <- netSEMp1(a, exogenous, endogenous)
  p1.df <- p1.result$bestModels
  p1.model <- p1.result$allModels

  p1.exogenous <- p1.result$bestModels %>% filter(Var == exogenous & Resp != endogenous) %>% arrange(adj.r.squared)
  p1.exogenous <- p1.exogenous %>% filter(row_number() == nrow(p1.exogenous))
  best_resp <- p1.exogenous$Resp #response variable (<S|R>) with the greatest adjR^2

  return(best_resp)

}