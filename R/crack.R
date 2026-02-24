#' Crack Quantification for Photovoltaic Backsheets
#' 
#' A dataset containing the average normalized crack depth for photovoltaic backsheets with inner layers of either ethylene-vinyl acetate or polyethylene exposed to 4,000 hours of continuous UVA irradiance with an intensity of 1.55 w/m2 at 340 nm and a chamber temperature of 70 deg C. 
#' Mechanistic variables from Fourier transform infrared spectroscopy are included to track chemical changes in the materials related to cracking. 
#' See the journal article titled 'A Non-Destructive Method for Crack Quantification in Photovoltaic Backsheets Under Accelerated and Real-World Exposures' in Polymer Degradation and Stability for more details. 
#' Doi: 10.1016/j.polymdegradstab.2018.05.008 
#' 
#' @docType data
#' @usage data(crack)
#' @author Addison G. Klinke, Abdulkerim Gok, Laura S. Bruckman, Roger H. French
#' 
#' @format A data frame with 97 rows and 5 variables:
#' \describe{
#'   \item{dAvgNorm}{ Average crack depth normalized by the backsheet's inner layer thickness}
#'   \item{uva360dose}{ Integrated, cumulative photodose for all wavelengths less than 360 nm}
#'   \item{crys730}{ Percent crystallinity calculated from the ratio of CH2 rocking peaks at 731 and 720 cm-1}
#'   \item{carb1715}{ Ketone carbonyl index calculated as the ratio of intensities at 1715 and 2851 cm-1}
#'   \item{carbPC1}{ 1st principal component score from the carbonyl region (1500-1800 cm-1)}
#' }
#' @source Solar Durability and Lifetime Extension (SDLE) Research Center, Case Western Reserve University
"crack"