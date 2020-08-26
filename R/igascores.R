#' Generate various scores for IgA binding in IgA-Seq experiments
#'
#' @description
#'
#' Calculate various different IgA-Seq scores across all the taxa and samples in an experiment.
#'
#' @details
#'
#' This function enables calculation of a variety of different indices for scoring immunoglobulin A (IgA) binding to taxa in IgA sequencing (IgA-Seq) experiments.
#' It is designed to be called on dataframes of abundance values, allowing easy calculation of scores across multiple taxa and samples.
#' The igaprobabilityratio(), igaprobability(), kauindex() and palmindex() functions can be used to calculate scores for one taxa and one sample.
#'
#' Scoring method can be chosen by specifying the method parameter as one of: "probratio", "prob", "kau", "palm" (Defaults to "probratio").
#' Each method requires different inputs as detailed below:
#'
#' \itemize{
#' \item probratio - equivalent to igaprobabilityratio() - requires two separate dataframes with iga positive abundances and iga negative abundances, two vectors with the sizes of the iga positive and negative gates per sample, and a pseudo count
#' \item prob - equivalent to igaprobability() - requires a dataframe with iga pos or neg fraction abundances, a vector of iga pos or neg gate size per sample, and a dataframe of taxa abundances in the presort samples
#' \item kau - equivalent to kauindex() - requires two separate dataframes with iga positive abundances and iga negative abundances, and a pseudo count
#' \item palm - equivalent to palmindex() - requires two separate dataframes with iga positive abundances and iga negative abundances, and a pseudo count
#' }
#'
#' @param posabunds A dataframe of taxa abundances in the positive/high IgA gate samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param negabunds A dataframe of taxa abundances in the negative/low IgA gate samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param pseudo The pseudo count to be used in scores. Default is 1e-5. Recommend setting to minimum observed abundance.
#' @param possizes A named vector containing the fraction of events in the IgA positive gate for each sample, with sample names matching abundance dataframes.
#' @param negsizes A named vector containing the fraction of events in the IgA negative gate for each sample, with sample names matching abundance dataframes.
#' @param presortabunds A dataframe of taxa abundances in the whole/initial samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param method Method to use to score IgA binding. One of: "probratio","prob","kau","palm". Default is "probratio".
#' @param scaleratio Should probratio scores be scaled to the pseudo count. Default is TRUE.
#' @param nazeros Should taxa with zero abundance in both the posabunds and negabunds (posabunds and presortabunds for prob method)  be scored as NA. Default is TRUE.
#' @return A data frame of IgA binding scores for all taxa and samples in the input data frame, generated using the scoring appraoch specified in 'method'.
#' @keywords iga score Kau Palm index ratio probability experiment iga-seq
#' @export
#' @examples
#' pab <- data.frame(Samp1=c(0.01,0.02,0.03),Samp2=c(0.05,0.02,0.04))
#' rownames(pab) <- c("Taxon1","Taxon2","Taxon3")
#' nab <- data.frame(Samp1=c(0.08,0.2,0.11),Samp2=c(0.05,0.0,0.07))
#' rownames(nab) <- c("Taxon1","Taxon2","Taxon3")
#' ps <- c(0.04,0.1)
#' ns <- c(0.08,0.4)
#' preab <- data.frame(Samp1=c(0.1,0.3,0.2),Samp2=c(0.15,0.05,0.2))
#' rownames(preab) <- c("Taxon1","Taxon2","Taxon3")
#' igascores(posabunds=pab,negabunds=nab, possizes=ps, negsizes=ns,pseudo=0.009)
#' igascores(posabunds=pab, possizes=ps, presortabunds=preab, method="prob")
#' igascores(posabunds=pab, negabunds=nab, pseudo=0.009, method="palm")
#' igascores(posabunds=pab, negabunds=nab, pseudo=0.009, method="kau")

igascores <- function(posabunds=NULL,negabunds=NULL,possizes=NULL,negsizes=NULL,pseudo=NULL,presortabunds=NULL, method="probratio", scaleratio=TRUE, nazeros=TRUE){
  methods <- c("probratio","prob","kau","palm")
  if(!method%in%methods){
    stop(paste0(method," is not a valid method."))
  }
  #data checks
  ##palm kau
  if(method=="palm" | method=="kau"){
    if(is.null(posabunds)|is.null(negabunds)|is.null(pseudo)){
      stop(paste0(method," method requires positive and negative abundance tables and a pseudocount. Specify posabunds, negabunds and pseudo."))
    }
    checkpseudo(pseudo)
    checktabs(posabunds,negabunds)
  }
  ##probratio method
  else if(method=="probratio"){
    if(is.null(posabunds)|is.null(negabunds)|is.null(pseudo)|is.null(possizes)|is.null(negsizes)){
      stop("probratio method requires positive and negative abundance tables, positive and negative gate sizes, and a pseudocount. Specify posabunds, negabunds, possizes, negsizes and pseudo.")
    }
    minval <- min(c(min(posabunds[posabunds!=0]),min(negabunds[negabunds!=0])))
    if(pseudo>minval){
      stop("Pseudo count must be lower than smallest non-zero abundance value.")
    }
    checkpseudo(pseudo)
    checktabs(posabunds,negabunds)
    checkvec(posabunds,possizes)
    checkvec(posabunds,negsizes)
  }
  ##prob method
  else if(method=="prob"){
    if(is.null(posabunds)&&is.null(negabunds)){
      stop("prob method requires either a positive or negative abundance table.")
    }
    if((!is.null(posabunds))&&(!is.null(possizes))){
      checkvec(posabunds,possizes)
      withinabunds <- posabunds
      gsizes <- possizes
      if(!is.null(negabunds)){
        message("Positive and negative abundances supplied to prob method. Will return positive gate probabilities. Remove positive variables to calculate negative gate probabilities.")}
    }else if((!is.null(negabunds))&&(!is.null(negsizes))){
      checkvec(negabunds,negsizes)
      withinabunds <- negabunds
      gsizes <- negsizes
    }
    else{
      stop("Must specify either posabunds and possizes or negabunds and negsizes for prob method.")
    }
    if(is.null(presortabunds)){
      stop("prob methods requires a total abundance table, specify presortabunds.")
    }
    checktabs(withinabunds,presortabunds)
  }

  #carry out calculations - methods specific for multiple samples/taxa (faster than calling the individual calculations within for loops)
  ##palm
  if(method=="palm"){
    nabunds <- negabunds
    nabunds[nabunds==0] <- nabunds[nabunds==0]+pseudo
    scores <- posabunds/nabunds
  }

  ##kau
  if(method=="kau"){
    pabunds <- posabunds+pseudo
    nabunds <- negabunds+pseudo
    nume <- log(pabunds)-log(nabunds)
    denom <- log(pabunds)+log(nabunds)
    scores <- -(nume/denom)
  }

  ##prob
  if(method=="prob"){
    nume <- rowprod(withinabunds,gsizes)
    denom <- presortabunds
    scores <- nume/denom
    #make NaN to NA
    scores[is.na(scores)] <- NA
    scores[scores>1] <- 1
  }

  ##probratio - default
  if(method=="probratio"){
    nume <- rowprod(posabunds,possizes)+pseudo
    denom <- rowprod(negabunds,negsizes)+pseudo
    scores <- log2(nume/denom)
    if(scaleratio==TRUE){
    scaler <- log2((1+pseudo)/pseudo)
    scores <- scores/scaler}
  }

  ##Convert scores to NA where no abundance detected
  if(nazeros==TRUE){
    if(method=="prob"){
     scores[presortabunds==0&withinabunds==0]=NA
    }else{
      scores[posabunds==0&negabunds==0]=NA
    }
  }

  return(scores)
}

#function to check pseudo and return erro
checkpseudo <- function(val){
  if((!is.numeric(val))|val<=0){
    stop("Pseudo count must be numeric and greater than zero.")
  }
}

#function to check table column and row names match
checktabs <- function(tab1,tab2){
  if((!is.data.frame(tab1))|(!is.data.frame(tab2))){
    stop("Abundance tables must be dataframes.")
  }
  if(any(!colnames(tab1)==colnames(tab2))){
    stop("Table column names do not match.")
  }
  if(any(!rownames(tab1)==rownames(tab2))){
    stop("Table row names do not match.")
  }
  if(any(tab1<0)|any(tab2<0)|any(tab1>1)|any(tab2>1)){
    stop("Abundance table values must be between 0 and 1.")
  }
}

#function to check the vector names match table column names
checkvec <- function(tab,vec){
  if(!is.data.frame(tab)){
    stop("Abundance tables must be dataframes.")
  }
  if(class(vec)!="numeric"|length(vec)!=ncol(tab)){
    stop("Gates sizes must be in numeric vectors matching abundance table columns.")
  }
  if(any(!colnames(tab)==names(vec))){
    stop("Size vector names do not match abundance table column names.")
  }
  if(any(tab<0)|any(vec<=0)|any(tab>1)|any(vec>1)){
    stop("Abundance table and gate size values must be between 0 and 1 (with gate sizes greater than 0).")
  }
}

#function to multiply rows of a dataframe by a vetor
rowprod <- function(tab,vec){
  return(data.frame(t(t(tab)*vec)))
}






