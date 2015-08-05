\name{NarrowPeaks-package}
\alias{NarrowPeaks-package}
\alias{NarrowPeaks}
\docType{package}

\title{ Shape-based Analysis of Variation in ChIP-seq using Functional PCA } 

\description{ The package applies a functional
        version of principal component analysis (FPCA) to: (1) Postprocess
        data in wiggle track format, commonly produced by generic ChIP-seq
        peak callers, by applying FPCA over a set of read-enriched
      regions (ChIP-seq peaks). This is done to study variability of the 
	the peaks, or to shorten their genomic 
        locations accounting for a given proportion of variation among
        the enrichment-score profiles. (2) Analyse
        differential variation between multiple ChIP-seq samples with 
	replicates. The function 'narrowpeaksDiff' quantifies differences
        between the shapes, and uses Hotelling's T2 tests on
        the functional principal component scores to identify 
	significant differences across conditions. An application of the package 
	for Arabidopsis datasets is described in Mateos, Madrigal, et al. 
	(2015) Genome Biology: 16:31. }
\details{
\tabular{ll}{
Package: \tab NarrowPeaks\cr
Type: \tab Package\cr
Version: \tab 1.13.3\cr
Date: \tab 2015-06-01\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
}

\author{
Pedro Madrigal, with contributions from Pawel Krajewski 
\email{pkra@igr.poznan.pl}

Maintainer: Pedro Madrigal \email{pm12@sanger.ac.uk}

}
\references{ Mateos JL, Madrigal P, et al. (2015) Combinatorial activities of SHORT VEGETATIVE PHASE and FLOWERING LOCUS C define distinct modes of flowering regulation in Arabidopsis. Genome Biology 16: 31.\cr
Bailey T, Krajewski P, Ladunga I, Lefebvre C, Li Q, Liu T, Madrigal P, Taslim C, Zhang J (2013) Practical Guidelines for the Comprehensive Analysis of ChIP-seq data. PLOS Comput Biol. 9 (11): e1003326. }

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
