#' Relative Abundance from Counts
#'
#' This function converts values in a dataframe to a fraction/percentage of the sum of their column.
#'
#' @param counttable Dataframe with rows as observations and columns as samples.
#' @param percentage Should values be returned as a percentage? i.e multiplied by 100. Default is FALSE (as required for most IgA scoring approaches).
#' @keywords abundance, normalisation, microbiome
#' @export

relabund <- function(counttable,percentage=FALSE){
  reltab <- data.frame(t(t(counttable)/colSums(counttable)))
  if(percentage==TRUE){
    reltab <- reltab*100
  }
  return(reltab)
}
