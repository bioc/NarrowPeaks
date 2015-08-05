\name{wigfile_test}
\docType{data}
\alias{wigfile_test}
\alias{NarrowPeaks-dataset}

\title{ Example Wiggle Track Produced After ChIP-seq Data Analysis}

\description{ Example of wiggle track produced after ChIP-seq data analysis. The 
data represents a small subset of a WIG file storing continuous-valued scores 
 based on a Poisson test for the chromosome 1 of \emph{Arabidopsis thaliana} 
(Kaufmann et al., 2010). It contains first 49515 lines of the WIG file for the
 complete experiment. }

\format{ Wiggle track format (WIG) data in a character vector. }

\source{ Gene Expression Omnibus GSE20176 
(\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20176}). Record from
 chromatin immunoprecipitation experiments wit AP1-specific antibodies followed
 by deep-sequencing in order to determine AP1 binding sites on a genome-wide 
scale in \emph{Arabidopsis thaliana}. }


\references{ Kaufmann et al. (2010) Orchestration of Floral Initiation by 
APETALA1. Science 328:85-89. }
\seealso{ \code{\link{NarrowPeaks-package}} }
\examples{ data(NarrowPeaks-dataset) }
\keyword{datasets}
