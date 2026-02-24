#' Aluminum Gradient Material for Metal's Design
#' 
#' Functional graded materials (FGM) are a class of materials with engineered continuous compositional gradients through the plate thickness. 
#' This work applies the netSEM approach on an aluminum FGM, produced via sequential alloy casting using planar solidification, to quantify the < Processing | Microstructure | Performance > relationships. 
#' The material has a continuous gradient in zinc (Zn) and magnesium (Mg) concentrations through the plate thickness. 
#' This subsequently produces a gradient in strengthening mechanisms from a dominant precipitate-strengthened aluminum alloy (AA) (Zn-based AA-7055) to a dominant strain-hardenable aluminum alloy (Mg-based AA-5456). 
#' Therefore, the material is simultaneously strengthened via solid solution strengthening and precipitation strengthening. 
#' 
#' @docType data
#' @usage data(metal)
#' @author Amit K. Verma, Roger H. French, Jennifer L. W. Carter 
#' 
#' @format A data frame with 72 rows and 6 variables:
#' \describe{
#'   \item{Hardness}{ Vickers hardness}
#'   \item{Z}{ the compositional gradient (z-direction)}
#'   \item{Mg}{ the element Zinc}
#'   \item{Zn}{ the element Magnesium }
#'   \item{MgZn2}{ the alloy }
#'   \item{Mgexcess}{ Mg-excess }
#' }
#' @source Solar Durability and Lifetime Extension (SDLE) Research Center, Case Western
#' Reserve University
"metal"