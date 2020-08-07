#' IgA Probability
#' @description
#' This function calculates the conditional probability that at batceria will be sufficently bound/not bound to IgA to end up in a given IgA gate based on its taxonomy. Calculated on one taxa for one sample.
#'
#' This uses Bayes' theorem assuming:
#' \itemize{
#' \item That the relative abundance of a given taxon in the IgA gate under question represents the probability of being that taxa given that it is within the IgA gate (either high or low).
#' \item That the percentage of flow cytometery events binned into the IgA gate represents the probability of any bacteria being within the gate.
#' \item That the abundance of the given taxon in the input sample (or whole fraction) represent the probability that any bacteria is assigned to the taxon.
#' If there is insufficient levels of a taxa in the whole fraction to account for its abundance in the IgA gate, the function assumes all of the taxa fall within this gate (i.e. a probability of 1).
#' }
#' @param withinabund Abundance of the bacteria in the IgA gate under investigation (can be calculated for either the pos/high or neg/low gating) (abundances should sum to 1 not as a \%).
#' @param gatesize The fraction of events in the flow cytometer within the gate under investigation (as a decimal fraction not a \%).
#' @param presortabund Abundance of the bacteria in whole sample before sorting by IgA (abundances should sum to 1 not as a \%).
#' @param nazeros Return NA if the within and tot abundances are both zero. Default is TRUE.
#' @keywords iga, probability, Jackson, iga-seq
#' @export
#' @examples
#' igaprobability(withinabund=0.5,gatesize=0.05,presortabund=0.5)

igaprobability <- function(withinabund,gatesize,presortabund,nazeros=TRUE){
  if(withinabund<0|gatesize<0|presortabund<0){
    stop("Abundances and gate size must be greater than or equal to zero.")
  }
  if(withinabund>1|gatesize>1|presortabund>1){
    stop("Abundances and gate size should be less than 1. Function expects values relative to 1 not 100 (i.e. not a percentage).")
  }
  if(withinabund==0&presortabund==0&nazeros==TRUE){
    return(NA)
  }
  nume <- withinabund*gatesize
  denom <- presortabund
  if(nume>denom){
    ip <- 1
  }else if(nume==0&&denom==0){
    ip <- NA
  }else{
    ip <- nume/denom
  }
 return(ip)
}
