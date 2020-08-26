#' Relative Abundance from Counts
#'
#' @description
#'
#' This function converts values in a dataframe to a fraction/percentage of the sum of their column.
#'
#' @param counttable Data frame of numeric values with rows as observations and columns as samples.
#' @param percentage Should values be returned as a percentage? i.e multiplied by 100. Default is FALSE (as required for most IgA scoring approaches).
#' @return A data frame of the input data normalised by column (to sum to either 1 or 100).
#' @keywords abundance normalisation microbiome relative
#' @export
#' @examples
#' taxcounts <- data.frame(Sample1=c(1,2,10,10),Sample2=c(3,10,5,1))
#' rownames(taxcounts) <- c("Taxon1","Taxon2","Taxon3","Taxon4")
#' relabund(taxcounts)

relabund <- function(counttable,percentage=FALSE){
  reltab <- data.frame(t(t(counttable)/colSums(counttable)))
  if(percentage==TRUE){
    reltab <- reltab*100
  }
  return(reltab)
}
