\name{narrowpeaksDiff}
\alias{narrowpeaksDiff}

\title{ Differential Analysis of Transcription Factor Binding using FPCA }
\description{ Shape-based differential binding analysis and testing for ChIP-Seq datasets 
	using Functional Principal Component Analysis and Hotelling's T2 tests. }

\usage{
narrowpeaksDiff(bedFile,  headerBed= TRUE, flank=100, bigwigs , conditions , 
nbasis=50, pcs = 10, bigWigSummaryPath=getwd(), variation = 0.6)
}

\arguments{

  \item{bedFile}{ Text file in BED format. It should contain at least 3 columns (chr, start, end), 
  	in which case the reference point for the FD analysis is calculated as the central point. A 4-th 
  	column can be provided containing the reference point.  }
  \item{headerBed}{ 'TRUE' if the first row in the BED file contain the name of the columns. 
  	'FALSE' otherwise.}
  \item{flank}{ Length (in bp.) that is considered upstream and downstream the reference point 
  	(or central point, if reference point is not given) for functional principal component analysis. }
  \item{bigwigs}{ Vector contaning the name of the bigWig files to be used in the analysis. }
  \item{conditions}{ Vector of characters with the Labels for the bigWig files. Biological replicates
  	must have the same label.}
  \item{nbasis}{ Number of order-4 B-Spline basis functions for functional data analysis.}
  \item{pcs}{ Number of principal components to be computed (default is 10). }
  \item{bigWigSummaryPath}{Path to the UCSC utility bigWigSummary (in case it is differente from the current 
  	directory). The tool can be downloaded for Linux and macOSX from the UCSC website: 
  	http://hgdownload.cse.ucsc.edu/admin/exe/}
  \item{variation}{Minimum proportion of varition that is considered to select the number of functional 
  	principal component scores used in the Hotelling's T2 tests (0-1, defualt is 0.6). }

}

\details{ Detailed information can be found in Madrigal and Krajewski (in preparation).}


\value{
  A list containing the following elements:
  \item{fdaprofiles}{ A list of matrices corresponding to the data of regions of interest (BED file) 
  	in the bigWig files. }
  \item{p.values}{ A list of pairwise comparisons between experimental conditions (taking into account 
  	replicates) for each region in the BED file. P-values are computed using the Hotelling's T2 test.}
  
}

\references{ Madrigal P, Krajewski P (in preparation) Shape-based Dimensionality Reduction 
	Analyses by Functional PCA Reveal Associations between First and Higher Order Components 
	in Next-Generation Sequencing Coverage Profiles. \cr
Bailey T, Krajewski P, Ladunga I, Lefebvre C, Li Q, Liu T, Madrigal P, Taslim C, Zhang J (2013)
 Practical guidelines for the comprehensive analysis of ChIP-seq data. PLoS Comput Biol
 9(11):e1003326. \cr
Ramsay, J.O. and Silverman, B.W. (2005) Functional Data Analysis. New York:
 Springer. 
}

\author{ Pedro Madrigal, \email{pm59@cam.ac.uk} }

\seealso{ \code{\link{narrowpeaks}}, \code{\link{NarrowPeaks-package}} }

\examples{

##library(NarrowPeaks)
##bigwigs <- c("A_rep1.bw","A_rep2.bw","A_rep3.bw","B_rep1.bw","B_rep2.bw","B_rep3.bw")
##conds   <- c("A","A","A","B","B","B")
##x <- narrowpeaksDiff(bedFile="regions.bed", bigwigs=bigwigs, conditions=conds, variation = 0.8)
##x$p.values

}
