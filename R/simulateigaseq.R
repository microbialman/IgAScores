#' Simulate an IgA-Seq dataset from a pre-defined set of IgA-binding probabilities
#'
#' @description
#'
#' Simulates IgA-Seq to create datasets with a defined binding distribution that can be used to test scoring method performance
#'
#' @details
#'
#' This function will generate a simulated immunoglobulin A sequencing (IgA-Seq) data set starting from a list containing the mean (and standard deviations) of IgA binding values expected for each species and cut-offs for defining the IgA positive and negative gates.
#' The input is a vector giving the average IgA value of each species (any arbitrary value that will represent the relative level of IgA binding between the species, ensure standard deviation and cut-offs are in the same magnitude).
#' These values are treated as the means of a normal distribution of IgA binding values for each species.
#' Species counts are generated on a log distribution for a given number of samples at an even depth.
#' For each bacteria in each sample, an IgA binding value is then assigned by sampling from its species IgA value distribution.
#' The value thresholds defining the positive and negative gates are then used to generate positive and negative counts tables of the bacteria whose values fall into these groups.
#' A second mode can also be used (by toggling betweengroups) that will introduce a consistent abundance change in half the samples by increasing one species in them. This can be used to simulate case-control experiments where, as an example,  one taxa has bloomed.
#' Further details can be found in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#'
#'Note: IgA values are simulated for each bacteria in each sample, setting the combination of the samplingdepth, number of species, and number of samples too high will slow the data generation.
#'
#' @param igavalmeans A vector of mean IgA values for as many species as you wish to simulate. Will default to an exponentially distributed vector of 10 species.
#' @param igavalsds A vector of standard deviations that will be used to generate IgA value distributions alongside the means. Defaults to 1 for all values.
#' @param nosamples The number of samples to generate simulated data from. Defaults to 10.
#' @param samplingdepth The number of bacteria to simulate in each sample. Defaults to 100000.
#' @param posthresh The IgA value threshold above which a bacteria will be considered IgA positive. Defaults to 4 (which is reasonable with the other defaults). It is recommended to run a simulation twice to determine reasonable thresholds on the first go.
#' @param negthresh The IgA value threshold below which a bacteria will be considered IgA negative. Defaults to 2 (which is reasonable with the other defaults). It is recommended to run a simulation twice to determine reasonable thresholds on the first go.
#' @param seed Seed for random number generation. Has a default so must be changed to rerun simulations.
#' @param betweengroups If TRUE this will modify starting abundances of half of the samples similarly (by adding betweenper\% of total counts to a single species) to simulate the case where there is an abundance shift without a change in IgA binding affinity. Defaults to FALSE.
#' @param betweenper Percentage of total counts to add to a species in the second group in the betweengroups mode.
#' @param betweensp Species (by index) to increased in between groups simulation. Chosen at random if NULL (default).
#' @return A list containing the simulated data set and relevant input parameters.
#' \itemize{
#'   \item presortcounts - A data frame containing simulated species counts for each sample in the pre-sort sample.
#'   \item presortabunds - presortcounts as relative abundances.
#'   \item poscounts - A data frame containing simulated species counts for each sample in the IgA positive fraction.
#'   \item posabunds - poscounts as relative abundances.
#'   \item negcounts - A data frame containing simulated species counts for each sample in the IgA negative fraction.
#'   \item negabunds - negcounts as relative abundances.
#'   \item possizes - A vector of the IgA positive fraction sizes for each sample.
#'   \item negsizes - A vector of the IgA negative fraction sizes for each sample.
#'   \item igabinding - A long format data frame containing the simulated IgA binding values for all simulated bacteria used to generate the count tables.
#'   \item igavalmeans - A vector of the mean IgA values for each species used in the simulation.
#'   \item igavalsds - A vector of the standard deviations of the IgA values for each species used in the simulation.
#'   \item posthresh - Numeric, the lower threshold used to determine a bacteria is IgA postive in the simulation.
#'   \item negthresh - Numeric, the upper threshold used to determine a bacteria is IgA negative in the simulation.
#'   \item expgroup - A vector showing class labels for the experimental group of each sample in the experiment. Will be uniform unless doing between group simulations.
#'   \item expspecies - Numeric, showing which species was modelled as differentially abundant between experimental groups when carryingout between group simulations.
#' }
#' @keywords iga iga-seq simulation benchmarking
#' @importFrom stats rlnorm rnorm rexp
#' @export
#' @examples
#' dat <- simulateigaseq(c(0.1,1,10,15),rep(1,4),posthresh=8,negthresh=4,samplingdepth=100)

simulateigaseq <- function(igavalmeans=NULL,igavalsds=NULL,nosamples=10,samplingdepth=100000,posthresh=4,negthresh=2,seed=66,betweengroups=FALSE, betweenper=10, betweensp=NULL){
  #use seed as there is a lot of random generation
  set.seed(seed)

  #set to default IgA mean value per species if none given - this value is just an arbitrary number relating to relative binding between species
  if(is.null(igavalmeans)){
  #by default use an exponential distribution of 10 values
  igavalmeans <- 2^rexp(10,rate=1)
  }
  #if no SD is given for each species IgA value distribution just use 1 for all
  if(is.null(igavalsds)){
  igavalsds=rep(1,length(igavalmeans))
  }
  #species names 1 to n
  speciesnames <- paste0("Species",c(1:length(igavalmeans)))
  names(igavalmeans) <- speciesnames
  nospecies <- length(speciesnames)
  #sample names 1 to n
  samplenames <- paste0("Sample",c(1:nosamples))

  #generate log distributed counts across the species within each sample
  #set of seeds to change the distribution for each sample
  aseeds <- sample(1:1e5,nosamples)
  #data frame to hold the counts
  speciesabunds <- matrix(ncol=nosamples,nrow=nospecies)
  for(i in 1:ncol(speciesabunds)){
    set.seed(aseeds[i])
    #generate a log normal distribution of species abundances in the samples
    speciesabunds[,i] <- rlnorm(nospecies,sdlog=2)
    #convert these relative scales to relative abundances and multiply by the total number of bacteria to get counts of each species in each sample
    speciesabunds[,i] <- round(speciesabunds[,i]/sum(speciesabunds[,i])*samplingdepth)
  }
  #level off to exactly the same depth (some samples will be slightly higher or lower due to rounding, just remove or add one to the first row of bacteria)
  speciesabunds[1,] <- speciesabunds[1,]+(samplingdepth-colSums(speciesabunds))
  #relabel
  rownames(speciesabunds) <- speciesnames
  colnames(speciesabunds) <- samplenames

  #if between groups is true we want to induce a similar change in half the samples starting abundances, for this we will add 0.1 to one of the species
  #simulating outgrowth of one taxa under a given condition
  expgroup <- rep(1,ncol(speciesabunds))
  names(expgroup) <- samplenames
  changesp <- NULL
  if(betweengroups==TRUE){
    #set half of the samples to be group 2
    g2size <- round(ncol(speciesabunds)/2)
    #choose a species to add if NULL
    if(is.null(betweensp)){
      changesp <- sample(c(1:nrow(speciesabunds)),1)
    }else{
      changesp <- betweensp
    }
    expgroup[(g2size+1):ncol(speciesabunds)] <- 2
    speciesabunds[changesp,(g2size+1):ncol(speciesabunds)] <- speciesabunds[changesp,(g2size+1):ncol(speciesabunds)] + (samplingdepth*(betweenper/100))
    #convert back to even depth
    speciesabunds <- round(t(t(speciesabunds)/colSums(speciesabunds))*samplingdepth)
    #level off again
    speciesabunds[1,] <- speciesabunds[1,]+(samplingdepth-colSums(speciesabunds))
    changesp <- speciesnames[changesp]
  }

  #for each individual bacteria in the counts generate above, assign it an IgA binding value that is sampled from the a normal distribution based on the mean and SD IgA value for its species
  sample <- c()
  species <- c()
  igabinding <- c()
  #for each species
  for(i in 1:nrow(speciesabunds)){
    #for each sample
    for(j in 1:ncol(speciesabunds)){
      #samples IgA binding values for the number of bacteria we see of this species in this sample, from the appropriate distribution
      igavals <- rnorm(speciesabunds[i,j],mean = igavalmeans[i], sd=igavalsds[i])
      #store in vecotrs of IgA value, Sample, and Species
      igabinding <- c(igabinding,igavals)
      sample <- c(sample,rep(colnames(speciesabunds)[j],length(igavals)))
      species <- c(species,rep(rownames(speciesabunds)[i],length(igavals)))
    }
  }
  #long format dataframe that can be used to plot the distributions later
  longigabinding <- data.frame(Sample=sample,Species=species,IgAValue=igabinding)

  #generate positive and negative fractions by counting the number of each species with values over and under the given thresholds
  poscounts <- filterigascorelong(longigabinding,posthresh,"over",speciesabunds)
  negcounts <- filterigascorelong(longigabinding,negthresh,"under",speciesabunds)

  #generate the neg and pos fraction sizes relative to total counts
  possizes <- colSums(poscounts)/colSums(speciesabunds)
  negsizes <- colSums(negcounts)/colSums(speciesabunds)

  #return: the presort counts table, the pos counts, the neg counts (and relative abundance equivalents), the pos fraction sizes, the neg fraction sizes, and the IgA binding table (for plotting density), as well as the input variables
  returnlist <- list(presortcounts=speciesabunds, presortabunds=relabund(speciesabunds), poscounts=poscounts, posabunds=relabund(poscounts), negcounts=negcounts, negabunds=relabund(negcounts), possizes=possizes, negsizes=negsizes, igabinding=longigabinding, igavalmeans=igavalmeans, igavalsds=igavalsds, posthresh=posthresh, negthresh=negthresh, expgroup=expgroup, expspecies=changesp)
  return(returnlist)
}



#function to filter IgA counts on intensity threshold
filterigascorelong <- function(longtab,threshold,direction,presortab){
  if(direction=="over"){filtlong <- longtab[longtab$IgAValue>threshold,]}
  if(direction=="under"){filtlong <- longtab[longtab$IgAValue<threshold,]}
  filtwide <- table(filtlong$Species,filtlong$Sample)
  filtwide <- as.data.frame.matrix(filtwide[match(rownames(presortab),rownames(filtwide)),colnames(presortab)])
  filtwide[is.na(filtwide)] <- 0
  rownames(filtwide) <- rownames(presortab)
  return(filtwide)
}
