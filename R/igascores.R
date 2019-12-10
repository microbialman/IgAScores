#' Generate various scores for IgA binding in IgA-Seq experiments
#'
#' @description
#' This function enables calculation of a variety of different indicies for scoring IgA binding to taxa in IgA-Seq experiments.
#' It is designed to be called on dataframes of abundance values, allowing easy calculation of scores across multiple taxa and samples.
#' The igaprobabilityratio(), igaprobability(), kauindex() and palmindex() functions can be used to calculate scores for one taxa and one sample.
#'
#' Scoring method can be chosen by specifying the method parameter as one of: "probratio", "prob", "kau", "palm" (Defaults to "probratio").
#' Each method requires different inputs as detailed below:
#'
#' \itemize{
#' \item probratio - equivalent to igaprobabilityratio() - requires two seperate dataframes with iga positive abundances and iga negative abundances, two vectors with the sizes of the iga postive and negative gates per sample, and a pseudo count
#' \item prob - equivalent to igaprobability() - requires a dataframe with iga pos or neg fraction abundances, a vector of iga pos or neg gate size per sample, and a dataframe of taxa abudnances in the whole/inital samples
#' \item kau - equivalent to kauindex() - requires two seperate dataframes with iga positive abundances and iga negative abundances, and a pseudo count
#' \item palm - equivalent to palmindex() - requires two seperate dataframes with iga positive abundances and iga negative abundances, and a pseudo count
#' }
#'
#' @param posabunds A dataframe of taxa abundances in the positive/high IgA gate samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param negabunds A dataframe of taxa abundances in the negative/low IgA gate samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param pseudo The pseudo count to be used in scores. Default is 1e-5. Recommend setting to minimum observed abundance.
#' @param possizes A named vector containing the fraction of events in the IgA postive gate for each sample, with sample names matching abundance dataframes.
#' @param negsizes A named vector containing the fraction of events in the IgA negative gate for each sample, with sample names matching abundance dataframes.
#' @param totabunds A dataframe of taxa abundances in the whole/initial samples. Samples as columns and taxa as rows, column and row names must match across abundance tables.
#' @param method Method to use to score IgA binding. One of: "probratio","prob","kau","palm". Default is "probratio".
#' @keywords iga, score, Kau, Palm, Jackson, index, ratio, probability, experiment, iga-seq
#' @export

igascores <- function(posabunds=NULL,negabunds=NULL,possizes=NULL,negsizes=NULL,pseudo=NULL,totabunds=NULL,method="probratio"){
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
        message("Positive and negative abundances supplied to prob method. Will return postive gate probabilities. Remove postive variables to calculate negative gate probabilities.")}
    }else if((!is.null(negabunds))&&(!is.null(negsizes))){
      checkvec(negabunds,negsizes)
      withinabunds <- negabunds
      gsizes <- negsizes
    }
    else{
      stop("Must specify either posabunds and possizes or negabunds and negsizes for prob method.")
    }
    if(is.null(totabunds)){
      stop("prob methods requires a total abundance table, specify totabunds.")
    }
    checktabs(withinabunds,totabunds)
  }

  #carry out calculations - methods specific for multiple samples/taxa (faster than calling the individual calculations within for loops)
  ##palm
  if(method=="palm"){
    negabunds[negabunds==0] <- negabunds[negabunds==0]+pseudo
    scores <- posabunds/negabunds
  }
  ##kau
  else if(method=="kau"){
    posabunds <- posabunds+pseudo
    negabunds <- negabunds+pseudo
    nume <- log(posabunds)-log(negabunds)
    denom <- log(posabunds)+log(negabunds)
    scores <- -(nume/denom)
  }
  ##prob
  else if(method=="prob"){
    nume <- rowprod(withinabunds,gsizes)
    denom <- totabunds
    scores <- nume/denom
    #make NaN to NA
    scores[is.na(scores)] <- NA
    scores[scores>1] <- 1
  }
  ##probratio - default
  else if(method=="probratio"){
    nume <- rowprod(posabunds,possizes)+pseudo
    denom <- rowprod(negabunds,negsizes)+pseudo
    scores <- log2(nume/denom)
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






