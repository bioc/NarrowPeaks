\name{narrowpeaks}
\alias{narrowpeaks}

\title{ Calculate Narrow Peaks from Enrichment-Score Profiles forming Broad
 Peaks }
\description{ Calculate narrow peaks from enrichment-score profiles forming
 broad peaks. }

\usage{
narrowpeaks(inputReg, scoresInfo, lmin = 0, nbf = 50, rpenalty= 0,  
nderiv= 0, npcomp = 5, pv = 80, pmaxscor = 0.0, ms = 0)
}

\arguments{
  \item{inputReg}{ Output of the function \code{sigWin} in package
 \code{CSAR}. }
  \item{scoresInfo}{ Output \code{infoscores} in the function 
\code{wig2CSARScore}, or the function \code{ChIPseqScore} after data analysis
 with package \pkg{CSAR}. }
  \item{lmin}{ Minimum length of an enriched region from the WIG file to be 
processed. Integer value. }
  \item{nbf}{ Number of order 4 B-spline basis functions that will represent 
the shape of each candidate site. Integer value.  }
  \item{rpenalty}{ Smoothing parameter for derivative penalization. Positive 
numeric value. }
  \item{nderiv}{ Order of derivative penalization, if \code{rpenalty}>0. 
Integer value. }
  \item{npcomp}{ Number of functional principal components. Integer value 
greater than or equal to \code{nbf}. }
  \item{pv}{ Minimum percentage of variation to take into account during the 
analysis. Numeric value in the range 0-100. }
  \item{pmaxscor}{ Cutoff for trimming of scoring function. Numeric value in 
the range 0-100. }
  \item{ms}{ Peaks closer to each other than \code{ms} nucleotides are to be 
merged in the final list. Integer value. }

}

\details{ This function produces shortened sites from a list of candidate 
transcription factor binding sites of arbitrary extension and shape. First, 
the enrichment signal from each candidate site is represented by a smoothed 
function constructed using a linear combination of order 4 B-spline basis 
functions. The data values are fitted using either least squares (if 
\eqn{rpenalty = 0}), or penalized residuals sum of squares (spline smoothing 
if \eqn{rpenalty > 0}). \cr Then, a functional principal component analysis 
for \code{npcomp} eigenfunctions is performed (Ramsay and Silverman, 2005), 
giving as a result a set of probe scores (principal component scores) which 
sum of squares is reported in \code{elementMetadata(broadPeaks)[,"fpcaScore"]
}.The higher the value of \code{fpcaScore}, the higher the variance that 
candidate peak accounts for within the original data. Details on the usage 
of semi-metrics in functional PCA is described in Ferraty and Vieu, 2006.\cr 
After that, we impose the condition that total scoring function for each 
reported narrow peak must be at least \code{pmaxscor} per cent of the maximum 
value. Max value is calculated from a set of scoring functions using only the 
eigenfunctions required to achieve \code{pv} percent of variance. A new set 
of scores is computed using trimmed versions of the eigenfunctions (Madrigal 
et al., submitted), and the root square is stored in 
\code{elementMetadata(narrowPeaks)[,"trimmedScore"]}.

}

\value{
  A list containing the following elements:
  \item{fdaprofiles}{ A functional data object encapsulating the enrichment 
profiles (see \pkg{fda} package. To plot the data use 
\code{plot.fd(fdaprofiles)}). }
  \item{broadPeaks}{ Description of the peaks prior to trimming. A 
\code{GRanges} object (see \pkg{GenomicRanges} package) with the information: 
\code{seqnames} (chromosome), \code{ranges} (start and end of the candidate 
site), \code{strand} (not used), \code{max} (maximum signal value for 
candidate site), \code{average} (mean signal value for candidate site), 
\code{fpcaScore} (sum of squares of the first \code{reqcomp} principal 
component scores for candidate site). }
  \item{narrowPeaks}{ Description of the peaks after trimming. A \code{GRanges}
 object (see \pkg{GenomicRanges} package) with the information: \code{seqnames}
 (chromosome), \code{ranges} (start and end after trimming), \code{strand} 
(not used), \code{broadPeak.subpeak}, \code{trimmedScore} (see details), 
\code{narrowedDownTo} (length reduction relative to the candidate), 
\code{merged} (logical value). }
  \item{reqcomp}{ Number of functional principal components used. Integer value.
}
  \item{pvar}{ Total proportion of variance accounted for by the \code{reqcomp}
 components used. Numeric value in the range 0-100 (always greater than or 
equal to argument \code{pv}). }
  
}

\references{ Madrigal, P. et al. (submitted) NarrowPeaks: an R/Bioconductor 
package for quantitative analysis of variation in ChIP-seq datasets. \cr
Ramsay, J.O. and Silverman, B.W. (2005) Functional Data Analysis. New York:
 Springer. \cr
Ferraty, F. and Vieu, P. (2006) Nonparametric Functional Data Analysis. 
New York: Springer. 
}

\author{ Pedro Madrigal, \email{pm@engineering.com} }

\seealso{ \code{\link{wig2CSARScore}}, \code{\link{NarrowPeaks-package}} }

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
