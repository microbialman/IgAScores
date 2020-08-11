#' Extract name at a given taxonomic level from a full name in the "p__;c__;o__;f__;g__;s__" format
#'
#' This function splits a full taxonomic lineage as a given level and returns the latter half.
#'
#' @param names Name string/ vector of name strings
#' @param level taxonomic level to split at must be in range phylum to species (default is genus).
#' @keywords microbiome taxonomy name split
#' @export
#' @examples
#' taxnamesplit("p__Bacteroidetes;c__Bacteroidia","class")


taxnamesplit <- function(names,level="genus"){
  delim <- paste0(substr(level,1,1),"__")
  spl <- strsplit(names,delim)
  nm <- unlist(lapply(spl,"[[",2))
  return(nm)
}
