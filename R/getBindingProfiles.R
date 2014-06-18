
getBindingProfiles <- function(win,scoreFiles){
  # Extract binding site profile-values from the CSAR binary files

	numberEnrichedRegions <- nrow(win)        #Ok
	numberOfChromosomes <- length(scoreFiles) #Ok
	profiles <- list()
	ChrScores <- list()
	#solve bug_2014
	ChrNames <- c() #chr names
	
	for(i in (1:numberOfChromosomes)){
		ChrScores[[i]] <- LoadBinCSAR(scoreFiles[i])
		#solve bug_2014
		ChrNames[i] <- strsplit(x=scoreFiles[i], split="_", fixed=TRUE)[[1]][1]
	}
	
	for(j in (1:numberEnrichedRegions)){
		lengthEnrichedRegion <- win[j,6]  
		startEnrichedRegion <- win[j,2]
		endEnrichedRegion <- win[j,3] 
		chromosome <- win[j,1]
		# chromosome can start be start as "Chr", "chr", or just the number
		#chromosomeNumber <- as.integer(sub("Chr","000",chromosome)) 
		chromosomeNumber <- which( ChrNames == as.character(chromosome) )
		profiles[[j]] <- ChrScores[[chromosomeNumber]][startEnrichedRegion:endEnrichedRegion]
	}

  return(profiles)
}
