
wig2CSARScore <- function(wigfilename,nbchr,chrle) {

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


  return(list(infoscores=returned_data))
}
