#' Species level counts for the OligoMM12-Colitis experiment used as an example in the IgAScores package
#'
#' @description
#'
#' Data from the colitis model described in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#'
#' @details
#'
#' Species level counts for an experiment where mice with a defined gut microbiota (OligoMM12) were either given *Helicobacter hepaticus* and IL10R antibody or the antibody alone (the first developing colitis).
#' Metadata for this experiment can be found in oligoMeta. Counts were generated from ASVs from the V4 region of the 16S rRNA gene processed using DADA2 and aligned to the *"RefSeq-RDP16S_v2_May2018.fa.gz"* database.
#' Further details can be found in Jackson et al. (2020, \doi{10.1101/2020.08.19.257501}).
#'
#' @docType data
#'
#' @usage data(oligoSpecies)
#'
#' @format An object of class \code{"tibble"}.
#'
#' @keywords dataset species oligoMM12 colitis ail10r helicobacter iga-seq counts
#'
#' @references To come...
#'
#' @examples
#' data(oligoSpecies)
"oligoSpecies"
