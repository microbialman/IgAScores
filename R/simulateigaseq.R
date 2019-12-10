#' Simulate an IgA-Seq dataset from a pre-defined set of IgA-binding probabilities
#'
#' @description
#'
#' This function will generate a simulated IgA-Seq data set starting from a list containing the mean (and standard devations) of IgA binding values expected for each species and cut-offs for defining the IgA postive and negative gates.
#' The input is a vector giving the average IgA value of each species (any arbitrary value that will represent the relative level of IgA binding between the species, ensure SD and cut-offs are in the same mangitude).
#' These values are treated as the means of a normal distribution of IgA binding values for each species.
#' Species counts are generated on a log distribution for a given number of samples at an even depth.
#' For each bacteria in each sample, an IgA binding value is then assigned by sampling from its species IgA value distribution.
#' The value thresholds defining the positive and negative gates are then used to generate positive and negative counts tables of the bacteria whose values fall into these groups.
#'
#'Note: IgA values are simulated for each bacteria in each sample, setting the combination of the samplingdepth, number of species, and number of samples too high will slow the data generation.
#'
#' @param igavalmeans A vector of mean IgA values for as many species as you wish to simulate. Will default to log normal distributed vector of 10 species.
#' @param igavalsds A vector of standard deviations that will be used to generate IgA value distributions alongside the means. Defaults to 10 for all values.
#' @param nosamples The number of samples to generate simulated data from. Defaults to 10.
#' @param samplingdepth The number of bacteria to simulate in each sample. Defaults to 100000.
#' @param posthresh The IgA value threshold above which a bacteria will be considered IgA positive. Defaults to 25 (which is reasonable only with the other defaults). It is recommended to run a simulation twice to determine reasonble thresholds on the first go.
#' @param negthresh The IgA value threshold below which a bacteria will be considered IgA negative. Defaults to 0 (which is reasonable only with the other defaults). It is recommended to run a simulation twice to determine reasonble thresholds on the first go.
#' @param seed Seed for random number generation. Has a default so must be changed to rerun simulations.
#' @keywords iga, iga-seq, simulation, benchmarking
#' @importFrom stats rlnorm rnorm
#' @export
#' @examples
#' dat <- simulateigaseq()
#' dat <- simulateigaseq(c(0.1,1,10,15),rep(1,4),nosamples=10,posthresh=8,negthresh=4)

simulateigaseq <- function(igavalmeans=NULL,igavalsds=NULL,nosamples=10,samplingdepth=100000,posthresh=50,negthresh=10,seed=808){
  #use seed as there is a lot of random generation
  set.seed(seed)

  #set to default IgA mean value per species if none given - this value is just an arbitrary number relating to relative binding between species
  if(is.null(igavalmeans)){
  #by default use a log normal distribution of 10 values
  igavalmeans <- rlnorm(10, sdlog = 2)
  }
  #if no SD is given for each species IgA value distribution just use 10 for all
  if(is.null(igavalsds)){
  igavalsds=rep(20,length(igavalmeans))
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

  #return: the whole counts table, the pos counts, the neg counts (and relative abundance equivalents), the pos fraction sizes, the neg fraction sizes, and the IgA binding table (for plotting density), as well as the input variables
  returnlist <- list(wholecounts=speciesabunds, wholeabunds=relabund(speciesabunds), poscounts=poscounts, posabunds=relabund(poscounts), negcounts=negcounts, negabunds=relabund(negcounts), possizes=possizes, negsizes=negsizes, igabinding=longigabinding, igavalmeans=igavalmeans, igavalsds=igavalsds, posthresh=posthresh, negthresh=negthresh)
  return(returnlist)
}



#function to filter IgA counts on intensity threshold
filterigascorelong <- function(longtab,threshold,direction,wholeab){
  if(direction=="over"){filtlong <- longtab[longtab$IgAValue>threshold,]}
  if(direction=="under"){filtlong <- longtab[longtab$IgAValue<threshold,]}
  filtwide <- table(filtlong$Species,filtlong$Sample)
  filtwide <- as.data.frame.matrix(filtwide[rownames(wholeab),colnames(wholeab)])
  return(filtwide)
}
