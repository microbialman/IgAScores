#' Palm Index
#'
#' @description
#'
#' This function calculates the immunoglobulin A (IgA) Index as defined in Palm et al. (2014, \doi{10.1016/j.cell.2014.08.006}) for a single taxon in a single sample.
#'
#' @param posabund Abundance of the bacteria in the IgA positive/high fraction.
#' @param negabund Abundance of the bacteria in the IgA negative/low fraction.
#' @param pseudo Pseudo count added to the abundance of the IgA negative fraction if the bacteria is not in that fraction. Defaults to 1e-5. Recommend setting to minimum observed abundance in whole dataset.
#' @param nazeros Return NA if the pos and neg abundances are both zero. Default is TRUE.
#' @return A numeric value for the Palm index as defined in Palm et al. (2014, \doi{10.1016/j.cell.2014.08.006}).
#' @keywords iga coating index Palm iga-seq
#' @export
#' @examples
#' palmindex(posabund=0.1,negabund=0.2,pseudo=0.0002)

palmindex <- function(posabund,negabund,pseudo=1e-5,nazeros=TRUE){
  if(posabund<0|negabund<0|pseudo<=0){
    stop("Postive and negative abundances must be greater than or equal to zero. Pseudo count must be greater than zero.")
  }
  if(posabund==0&negabund==0&nazeros==TRUE){
    return(NA)
  }
  if(negabund==0){
    negabund <- negabund+pseudo
  }
  ici <- posabund/negabund
  return(ici)
}
