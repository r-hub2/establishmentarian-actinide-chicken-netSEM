#' A data frame of an degradation experiment of poly(ethylene-terephthalate) films
#'
#' The data set is a study of photolysis and hydrolysis of poly(ethylene-terephthalate) films that contain an ultraviolet stabilizer additive. 
#' In this work, polymeric samples were exposed to UV light and moisture according to ASTM G-154 Cycle 4 standard accelerated weathering conditions. 
#' Resulting optical chemical changes were determined through optical and infrared (IR) spectroscopy. *Time* is the main exogenous variable and *YI* (yellowness index) is the endogenous variable (response). 
#' The other columns in the data set (*AbsEdge*, *UV.Stab.Bleaching*, *Crystallization*, and *ChainScission*) are values extracted from optical and IR absorbance spectra as single metrics and used as intermediate unit level endogenous (response) variables in the netSEM analysis.
#'
#' @docType data
#' @usage data(pet)
#' @author Devin A. Gordon, Wei-Heng Huang, Roger H. French, Laura S. Bruckman
#' 
#'
#' @format A 37 by 6 data frame of continuous variables: 
#' \describe{
#'   \item{YI}{Yellowness index of PET film}
#'   \item{Time}{Time exposed to ASTM G-154 Cycle 4 conditions}
#'   \item{AbsEdge}{Fundamental Absorption Edge - Optical absorbance at 312 nm}
#'   \item{UVStabBleaching}{Ultraviolet Stabilizer Bleaching - Optical absorbance at 340 nm}
#'   \item{Crystallization}{IR spectropscopy peak of relative crystallinity - IR absorbance at 975 wavenumber}
#'   \item{ChainScission}{IR carbonyl peak  - IR absorbance at 1711 wavenumber}
#' }
#' 
#' @source Solar Durability and Lifetime Extension (SDLE) Research Center, Case Western
#' Reserve University
#' 
"pet"
