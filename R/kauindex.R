#' Kau Index
#'
#' @description
#'
#' This function calculates the immunoglobulin A (IgA) Index as defined in Kau et al. (2015, \doi{10.1126/scitranslmed.aaa4877}) for a single taxon in a single sample.
#'
#' @param posabund The abundance of the bacteria in the IgA positive/high fraction (abundances should sum to 1 not as a \%).
#' @param negabund The abundance of the bacteria in the IgA negative/low fraction (abundances should sum to 1 not as a \%).
#' @param pseudo Pseudo count added to both the IgA positive and negative fraction values prior to calculation. Defaults to 1e-5. Recommend setting to minimum observed abundance in whole dataset.
#' @param nazeros Return NA if the pos and neg abundances are both zero. Default is TRUE.
#' @return A numeric value for the Kau index as defined in Kau et al. (2015, \doi{10.1126/scitranslmed.aaa4877}).
#' @keywords iga index Kau iga-seq
#' @export
#' @examples
#' kauindex(posabund=0.1,negabund=0.2,pseudo=0.0002)

kauindex <- function(posabund,negabund,pseudo=1e-5,nazeros=TRUE){
  if(posabund<0|negabund<0|pseudo<=0){
    stop("Postive and negative abundances must be greater than or equal to zero. Pseudo count must be greater than zero.")
  }
  if(posabund>1|negabund>1){
    stop("Abundance values should be less than 1. Function expects abundances relative to 1 not 100 (i.e. not a percentage).")
  }
  if(posabund==0&negabund==0&nazeros==TRUE){
    return(NA)
  }
  posabund <- posabund+pseudo
  negabund <- negabund+pseudo
  nume <- log(posabund)-log(negabund)
  denom <- log(posabund)+log(negabund)
  ii <- -(nume/denom)
  return(ii)
}
