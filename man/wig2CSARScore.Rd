\name{wig2CSARScore}
\alias{wig2CSARScore}

\title{ Convert Data from a Wiggle Track (WIG) File to CSAR Binary Format }
\description{ Convert data from a wiggle track (WIG) file to CSAR binary 
format and extract read-enriched regions. }

\usage{ wig2CSARScore(wigfilename, nbchr, chrle) }

\arguments{
  \item{wigfilename}{ WIG file containing the enrichment-score signal of a 
transcription factor binding experiment. }
  \item{nbchr}{ Number of chromosomes. }
  \item{chrle}{ Vector of lengths of the chromosomes (in base pairs). }
}

\details{ The Wiggle format (WIG/bigWig) is described on the UCSC Genome 
Bioinformatics web site: \url{http://genome.ucsc.edu/FAQ/FAQformat}. It allows
 the display of continuous value data in the genome browser. Although 
specifically designed for post-processing of WIG files, resulting from the 
analysis of ChIP-seq experiments (with Bioconductor packages \pkg{BayesPeak},
 \pkg{CSAR}, \pkg{PICS}, or other tools such as MACS, F-seq, etc.), 
\pkg{NarrowPeaks} can process other type of sequencing data encoded in WIG 
format in order to locate regions of high variability in the data. }

\value{
  A list of two elements:
   \item{infoscores}{A list with the same elements as reported by the function
 \code{ChIPseqScore} in the \pkg{CSAR} Bionductor package: chr (Chromosome
 names), chrL (Chromosome length (bp).), filenames (Name of the files where 
the score values are stored.), digits (Score values stored on the files need 
to be divided by 10^\code{digits}).}

}

\references{ Mateos JL, Madrigal P, et al. (2015) Combinatorial activities of SHORT VEGETATIVE PHASE and FLOWERING LOCUS C define distinct modes of flowering regulation in Arabidopsis. Genome Biology 16: 31.\cr
Bailey T, Krajewski P, Ladunga I, Lefebvre C, Li Q, Liu T, Madrigal P, Taslim C, Zhang J (2013) Practical Guidelines for the Comprehensive Analysis of ChIP-seq data. PLOS Comput Biol. 9 (11): e1003326. \cr
Muino JM, Kaufmann K, van Ham RC, Angenent GC, Krajewski P (2011) ChIP-seq analysis in 
R (CSAR): An R package for the statistical detection of protein-bound genomic regions. 
Plant Methods 7:11.
}

	

\author{ Pedro Madrigal, \email{dnaseiseq@gmail.com}}

\seealso{ \code{\link{narrowpeaks}}, \code{\link{NarrowPeaks-package}} }

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

setwd(owd)
}
