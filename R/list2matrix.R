
list2matrix <- function(coefList){

	L <- length(coefList)
	Lelements <- c()
	
	for (i in (1:L)){ Lelements[i] <- length(coefList[[i]]) }

	m <- max(Lelements)
	matrixData <- matrix(0,L,m)
	matrixMask <- matrix(0,L,m) 

	for (j in (1:L))
	{
		lj <- Lelements[j] 
		zerosAround <- m-lj
		startingIndex <- ceiling(zerosAround/2)
		matrixData[j,(startingIndex+1):(startingIndex+lj)] <- coefList[[j]]
		matrixMask[j,(startingIndex+1):(startingIndex+lj)] <- 1
	}

	return(list(dataMatrix = matrixData, dataMask = matrixMask))
}
