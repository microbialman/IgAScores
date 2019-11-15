#' IgA Probability Ratio
#'
#' This function calculates the ratio of the IgA positive fraction probability relative to the IgA negative fraction probability for a single taxa in a single sample.
#' These probabilities can individually be calculated using the igaprobability() function. As both calculations have the whole fraction taxon abundance as a denominator it cancels.
#' This means the IgA probiability ratio can be caluclated without this information.
#'
#' @param posabund Abundance of the bacteria in the IgA positive/high fraction (abundances should sum to 1 not as a \%).
#' @param negabund Abundance of the bacteria in the IgA negative/low fraction (abundances should sum to 1 not as a \%).
#' @param possize The fraction of events in the flow cytometer classed as IgA postive/high (as a decimal fraction not a \%).
#' @param negsize The fraction of events in the flow cytometer classed as IgA negative/low (as a decimal fraction not a \%).
#' @param psuedo Pseudo count added to both the IgA positive and negative abundance values prior to calculation. Defaults to 1e-5. Recommend setting to minimum observed abundance in whole dataset.
#' @keywords iga, probability, ratio, Jackson, iga-seq
#' @export
#' @examples
#' igaprobabilityratio(posabund=0.2,negabund=0.05,possize=0.05,negsize=0.6,pseudo=0.0002)

igaprobabilityratio <- function(posabund,negabund,possize,negsize,pseudo=1e-5){
  if(posabund<0|negabund<0){
    stop("Postive and negative abundances must be greater than or equal to zero.")
  }
  if(possize<=0|negsize<=0|pseudo<=0){
    stop("Positive and negative gate fractions and pseudocount must be greater than zero.")
  }
  if(posabund>1|negabund>1|possize>1|negsize>1){
    stop("Abundance and fraction values should be less than 1. Function expects values relative to 1 not 100 (i.e. not a percentage).")
  }
  posabund <- posabund+pseudo
  negabund <- negabund+pseudo
  nume <- posabund*possize
  denom <- negabund*negsize
  ipr <- log(nume/denom)
  return(ipr)
}
