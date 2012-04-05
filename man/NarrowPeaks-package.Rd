\name{NarrowPeaks-package}
\alias{NarrowPeaks-package}
\alias{NarrowPeaks}
\docType{package}

\title{ Functional Principal Component Analysis to Narrow Down Transcription 
Factor Binding Site Candidates }
\description{ The package processes data in wiggle track format (WIG) commonly
 produced by several ChIP-seq data analysis tools by applying functional 
version of principal component analysis (FPCA) over a set of selected candidate
 enriched regions. This is done in order to shorten the genomic locations 
accounting for a given proportion of variation among the enrichment-score
 profiles. It allows the user to discriminate between binding regions in close
 proximity to each other and to narrow down the length of the putative 
transcription factor binding sites while preserving the information present in
 the variability of the dataset and capturing major sources of variation. }
\details{
\tabular{ll}{
Package: \tab NarrowPeaks\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-03-15\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
}

\author{
Pedro Madrigal, with contributions from Pawel Krajewski 
\email{pkra@igr.poznan.pl}

Maintainer: Pedro Madrigal \email{pm@engineering.com}

}
\references{ Madrigal, P. et al. (submitted) NarrowPeaks: an R/Bioconductor 
package to narrow down transcription factor binding site candidates using 
functional PCA.}

\examples{
owd <- setwd(tempdir())

##For this example we will use a subset of the AP1 ChIP-seq data (Kaufmann et
##al., 2010)
##The data is obtained after analysis using the CSAR package available in 
##Bioconductor 
data("NarrowPeaks-dataset")
writeLines(wigfile_test, con="wigfile.wig")

##Write binary files with the WIG signal values for each chromosome 
##independently and obtain regions of read-enrichment with score values greater
##than 't', allowing a gap of 'g'. Data correspond to enriched regions found up
##to 105Kb in the Arabidopsis thaliana genome
wigScores <- wig2CSARScore(wigfilename="wigfile.wig", nbchr = 1, 
chrle=c(30427671))
gc(reset=TRUE) 
library(CSAR)
candidates <- sigWin(experiment=wigScores$infoscores, t=1.0, g=30)

##Narrow down ChIPSeq enriched regions by functional PCA
shortpeaks <- narrowpeaks(inputReg=candidates, 
scoresInfo=wigScores$infoscores, lmin=0, nbf=150, rpenalty=0, 
nderiv=0, npcomp=2, pv=80, pmaxscor=3.0, ms=0)

###Export GRanges object with the peaks to annotation tracks in various 
##formats. E.g.:
library(GenomicRanges)
names(elementMetadata(shortpeaks$broadPeaks))[3] <- "score"
names(elementMetadata(shortpeaks$narrowPeaks))[2] <- "score"
library(rtracklayer)
export.bedGraph(object=candidates, con="CSAR.bed")
export.bedGraph(object=shortpeaks$broadPeaks, con="broadPeaks.bed")
export.bedGraph(object=shortpeaks$narrowPeaks, con="narrowpeaks.bed")

setwd(owd)
}
