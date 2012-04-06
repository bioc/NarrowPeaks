
wig2CSARScore <- function(wigfilename,nbchr,chrle) {
  # Interfacing C code
  returned_data = .Call( 'wig2CSARScore_R_wrapper', 
        	   file_Name=as.character(wigfilename),
                   nbChr=as.integer(nbchr), 
                   chrL=as.integer(chrle))

  names(returned_data) <- c("chr", "chrL", "filenames", "digits")

  returned_data$chrL <- as.double(returned_data$chrL)     #CSAR specif

  return(list(infoscores=returned_data))
}
