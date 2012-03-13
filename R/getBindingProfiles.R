
getBindingProfiles <- function(win,scoreFiles){
  # Extract binding site profile-values from the CSAR binary files

	numberEnrichedRegions <- nrow(win)
	numberOfChromosomes <- length(scoreFiles)
	profiles <- list()
	ChrScores <- list()
	
	for(i in (1:numberOfChromosomes)){
		ChrScores[[i]] <- LoadBinCSAR(scoreFiles[i])
	}
	
	for(j in (1:numberEnrichedRegions)){
		lengthEnrichedRegion <- win[j,6]  
		startEnrichedRegion <- win[j,2]
		endEnrichedRegion <- win[j,3] 
		chromosome <- win[j,1]
		chromosomeNumber <- as.integer(sub("Chr","000",chromosome)) 
		profiles[[j]] <- ChrScores[[chromosomeNumber]][startEnrichedRegion:endEnrichedRegion]
	}

  return(profiles)
}
