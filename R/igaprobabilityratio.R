#' IgA Probability Ratio
#'
#' @description
#'
#' Calculate the IgA Probability Ratio score as described in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#'
#' @details
#'
#' This function calculates the ratio of the immunoglobulin A (IgA) positive fraction probability relative to the IgA negative fraction probability for a single taxa in a single sample.
#' These probabilities can individually be calculated using the igaprobability() function. As both calculations have the whole fraction taxon abundance as a denominator it cancels.
#' This means the IgA probability ratio can be calculated without this information.
#' Further details can be found in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#'
#' @param posabund Abundance of the bacteria in the IgA positive/high fraction (abundances should sum to 1 not as a \%).
#' @param negabund Abundance of the bacteria in the IgA negative/low fraction (abundances should sum to 1 not as a \%).
#' @param possize The fraction of events in the flow cytometer classed as IgA positive/high (as a decimal fraction not a \%).
#' @param negsize The fraction of events in the flow cytometer classed as IgA negative/low (as a decimal fraction not a \%).
#' @param pseudo Pseudo count added to both the IgA positive and negative abundance values prior to calculation. Defaults to 1e-5. Recommend setting to minimum observed abundance in whole dataset.
#' @param scaleratio Should probratio scores be scaled to the pseudo count. Default is TRUE.
#' @param nazeros Return NA if the pos and neg abundances are both zero. Default is TRUE.
#' @return A numeric value for the IgA Probability Ratio as defined in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#' @keywords iga probability ratio iga-seq
#' @export
#' @examples
#' igaprobabilityratio(posabund=0.2,negabund=0.05,possize=0.05,negsize=0.6,pseudo=0.0002)

igaprobabilityratio <- function(posabund,negabund,possize,negsize,pseudo=1e-5,scaleratio=TRUE,nazeros=TRUE){
  if(posabund<0|negabund<0){
    stop("Positive and negative abundances must be greater than or equal to zero.")
  }
  if(possize<=0|negsize<=0|pseudo<=0){
    stop("Positive and negative gate fractions and pseudocount must be greater than zero.")
  }
  if(posabund>1|negabund>1|possize>1|negsize>1){
    stop("Abundance and fraction values should be less than 1. Function expects values relative to 1 not 100 (i.e. not a percentage).")
  }
  if(posabund==0&negabund==0&nazeros==TRUE){
    return(NA)
  }
  vals <- c(posabund,negabund)
  minval <- min(vals[vals!=0])
  if(pseudo>minval){
    stop("Pseudo count must be smaller than lowest non-zero abundance.")
  }
  nume <- (posabund*possize)+pseudo
  denom <- (negabund*negsize)+pseudo
  ipr <- log2((nume/denom))
  scaler <- 1
  if(scaleratio==TRUE){
  scaler <- log2((1+pseudo)/pseudo)}
  return(ipr/scaler)
}
