\name{NarrowPeaks-package}
\alias{NarrowPeaks-package}
\alias{NarrowPeaks}
\docType{package}

\title{ Shape-based Analysis of Variation in ChIP-Seq using Functional PCA } 

\description{ The package applies a functional version of principal component 
analysis (FPCA) to: (1) Process data in wiggle track format (WIG) commonly produced 
by ChIP-Seq peak callers by applying FPCA over a set of read-enriched regions 
(ChIP-Seq peaks). This is done in order to shorten the genomic locations accounting 
for a given proportion of variation among the enrichment-score profiles. The function
 'narrowpeaks' allows splitting and trimming binding sites in close proximity to each 
 other, narrowing down the length of the putative transcription factor binding sites 
 while preserving the information present in the variability of the dataset and capturing 
 major sources of variation. (2) Analyse differential variation between multiple ChIP-Seq 
 samples with replicates. The function 'narrowpeaksDiff' quantifies differences between the 
 shapes, and uses Hotelling's T2 tests on the functional principal component scores to 
 identify significant differences across conditions. }
\details{
\tabular{ll}{
Package: \tab NarrowPeaks\cr
Type: \tab Package\cr
Version: \tab 1.11.4\cr
Date: \tab 2015-02-01\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
}

\author{
Pedro Madrigal, with contributions from Pawel Krajewski 
\email{pkra@igr.poznan.pl}

Maintainer: Pedro Madrigal \email{pm59@cam.ac.uk}

}
\references{ Madrigal P, Krajewski P (in preparation) Shape-based Dimensionality Reduction 
	Analyses by Functional PCA Reveal Associations between First and Higher Order Components 
	in Next-Generation Sequencing Coverage Profiles.}

\examples{
owd <- setwd(tempdir())

##For this example we will use a subset of the AP1 ChIP-Seq data (Kaufmann et
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
