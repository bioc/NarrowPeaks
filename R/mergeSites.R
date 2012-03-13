
mergeSites <- function (filebs, profiles_for_analysis, candidates_position, iniLength, N=2){

	FILE <- filebs[order(filebs$chrom ,filebs$siteStart),]
	for (i in 1:nrow(FILE)){
		tps <- FILE[i,"broadPeak"]
		position <- candidates_position[tps,] 
		profile <- which(profiles_for_analysis==tps)[1]
	}
	j=1
	overwrite=0
	matrix <- data.frame(FILE[1,])

	for (o in levels(FILE$chrom)){
		FILEtps<-subset(FILE,FILE$chrom==o)
		matrix[j,] <- data.frame(FILEtps[1,])
		matrix$merged[j] <- "FALSE"
		
		for (i in 2:dim(FILEtps)[1]){
			if ((FILEtps$siteStart[i]<=(matrix$siteEnd[j]+N))){
				overwrite=overwrite+1
				matrix$mergged[1] <- "FALSE"
				if (matrix$siteEnd[j]>FILEtps$siteEnd[i]){
				    matrix$merged[j] <- "FALSE"
				}
				else{
				    matrix$siteEnd[j]<-FILEtps$siteEnd[i]
				    matrix$broadPeak.subpeak[j]<- paste( matrix$broadPeak.subpeak[j],FILEtps$broadPeak.subpeak[i], sep='-')
				    matrix$fpcaScore[j] <- ( matrix$fpcaScore[j] * matrix$length[j] + FILEtps$fpcaScore[i] * FILEtps$length[i] ) / ( 1 + matrix$siteEnd[j] - matrix$siteStart[j] )
				    matrix$trimmedScore[j] <- matrix$trimmedScore[j] + FILEtps$trimmedScore[i]     # If merged, final trimmed score is the sum
				    matrix$length[j]<-1+matrix$siteEnd[j]-matrix$siteStart[j]
				    matrix$merged[j] <- "TRUE"
				    if(matrix$broadPeak[j] == FILEtps$broadPeak[i]) { matrix$narrowedDownTo[j] <-  paste(round(100*matrix$length[j]/ iniLength[matrix$broadPeak[j]],2),"%",sep="") }
				    if(matrix$broadPeak[j] != FILEtps$broadPeak[i]) { matrix$narrowedDownTo[j] <-  paste(round(100*matrix$length[j]/ (iniLength[matrix$broadPeak[j]] + iniLength[FILEtps$broadPeak[i]]),2),"%",sep="") } 
				}
			}
			else{
				j<-j+1
				matrix[j,] <- data.frame(FILEtps[i,])
				matrix$merged[j] <- "FALSE"
			}
		}
		j<-j+1
	}
	return(matrix)
}
