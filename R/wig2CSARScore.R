
wig2CSARScore <- function(wigfilename,nbchr,chrle, thr=1.0, gap=30) {

  Filename<-"a";
  length(Filename)<-nbchr
  Chrname<-"a";
  length(Chrname)<-nbchr
  
  # Interfacing C code
  returned_data = .C( 'wig2CSARScore_R_wrapper', 
        	   file_Name=as.character(wigfilename),
                   nbChr=as.integer(nbchr), 
                   chrL=as.integer(chrle),
                   filenames=as.character(Filename),
                   chr=as.character(Chrname),
                   digits=integer(1) );
 
  returned_data <- returned_data[c("chr","chrL","filenames","digits")]

  returned_data$chrL <- as.double(returned_data$chrL)     #CSAR specif
  returned_data$digits <- as.double(returned_data$digits) #CSAR specif

  #wrapper invoking CSAR::sigWin and converting the result to GenomicRanges::GRanges
  #Obtain regions of read-enrichment with score values greater than 't', allowing a gap of 'g'
  enrichedRegions <- sigWin(experiment=returned_data, t=thr, g=gap)
  #Genomic locations in GenomicRanges class GRanges
  CSARcandidates <- GRanges(seqnames=Rle(as.character(enrichedRegions$chr)),ranges=(IRanges(enrichedRegions$start, end=enrichedRegions$end)),posPeak=as.integer(enrichedRegions$posPeak), score=as.integer(enrichedRegions$score) )
  seqlengths(CSARcandidates) <- returned_data$chrL

  return(list(infoscores=returned_data,candidates=CSARcandidates))
}
