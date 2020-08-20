## IgAScores

The IgAScores package is used to calculate taxon-level IgA binding scores from IgA-Seq data.
It also includes several helper functions for managing microbiome data and simulating IgA-Seq datasets for method testing. 

For a detailed consideration of the methods that are used for scoring IgA-Seq data, see the [associated paper](https://biorxiv.org/cgi/content/short/2020.08.19.257501v1).

### Install

This package can be installed from GitHub using [devtools](https://github.com/r-lib/devtools).

```r
if (!require(devtools)){install.packages('devtools')}
install_github("microbialman/IgAScores", build_vignettes = TRUE)
```

### Running

A detailed overview of how to run the various functions within IgAScores can be found in the R vignette: `vignette("IgAScores")`

A brief summary of the main *igascores()* function is given below:

```r
#load in IgAScores
library(IgAScores)

#dataframes with counts for the bacterial taxa in the IgA+ and IgA- fractions, as would be produced by 16S rRNA appraoches such as DADA2
igapos <- data.frame(Sample1=c(100,0,1,2,10),Sample2=c(110,0,11,42,50),Sample3=c(140,60,10,3,0))
iganeg <- data.frame(Sample1=c(200,0,40,20,4),Sample2=c(10,30,110,2,5),Sample3=c(30,20,0,123,20))

taxnames <- c("Taxon1","Taxon2","Taxon3","Taxon4","Taxon5")
rownames(igapos) <- taxnames
rownames(iganeg) <- taxnames

#convert the counts to relative abundances using the included helper function
igapos <- relabund(igapos)
iganeg <- relabund(iganeg)

#iga+ and iga- fraction sizes per sample (fraction, if a percentage divide by 100)
possize <- c(Sample1=0.04,Sample2=0.05,Sample3=0.03)
negsize <- c(Sample1=0.54,Sample2=0.47,Sample3=0.33)

#set a pseudo count for handling zero values in some scoring methods
#this should be of a similar value of the minimum non-zero observed value (e.g. if minum values is 0.007 use 0.001)
pseudo <- 0.001

#default method is the probability ratio "probratio"", additional methods available are "prob"", "kau"" and "palm".
prscores <- igascores(posabunds = igapos, negabunds = iganeg, 
                      possizes = possize, negsizes = negsize, 
                      pseudo = pseudo)

print(prscores)

```
