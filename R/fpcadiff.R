

fpcadiff <- function(bedFile,  headerBed= TRUE, flank=100, bigwigs , conditions , nbasis=50, pcs = 10, bigWigSummaryPath=getwd(), variation = 0.6) {
 

	# Read BED file (chr start end  ref)
	bed <- read.table(bedFile, header = headerBed);
	if (dim(bed)[2] > 3){
		bed <- bed[,1:4]
		colnames(bed) <- c("chr", "start", "end", "ref")
		bed$start  <- bed$ref  - flank
		bed$end    <- bed$ref  + flank
	}	
	
	if (dim(bed)[2] == 3){   # If only chr,start,end; then we take the midpoint as reference
		bed <- bed[,1:3]
		colnames(bed) <- c("chr", "start", "end")
		ref <- round((bed$start + bed$end)/2)
		bed$start  <- ref  - flank
		bed$end    <- ref  + flank		
	}		


	fdamatrix <- list() #To store profiles of all datasets
	
	for (k in 1:length(bigwigs) ){

		fdamatrix[[k]] <- matrix(0.0,ncol=1+2*flank, nrow= length(bed$chr) )   

		for (i in 1:length(bed$chr)) {

			#owd <- setwd(tempdir())
			randString <- getRandString()
			UCSC <- paste(bigWigSummaryPath,"/bigWigSummary", sep="")
        	cmmd <- paste(UCSC, bigwigs[k], bed$chr[i], as.integer(bed$start[i]), as.integer(bed$end[i]), as.integer(1+2*flank), paste("> temp",".txt",sep=randString), sep=' ')

        	x <- system(cmmd, intern=FALSE);

        	if(length(readLines(paste("temp",".txt",sep=randString )))>0) {   #if no data the file is empty
			#if(length(x)==1) {

                x <- as.numeric( strsplit(x= readLines(paste("temp",".txt",sep=randString )) ,split="\t")[[1]]) 

				fdamatrix[[k]][i,] <- x
			#}
	        }

			system(paste("rm",paste("temp",".txt",sep=randString) ,sep=' '), intern=FALSE);
	        rm(x,randString)
		}
		
		# here we have the matrix with all the profiles
		# correct non-numeric values in case of NAs
		fdamatrix[[k]][which(is.na(fdamatrix[[k]])==TRUE)] <- 0.0
		fdamatrix[[k]][which(is.nan(fdamatrix[[k]])==TRUE)] <- 0.0
		fdamatrix[[k]][which(is.numeric(fdamatrix[[k]])==FALSE)] <- 0.0
		
	}
	
	#List of matrices with the profile data ready#

	#Create list of data.frames for the p-values
	PVALS <- list()
	uniqueCond <- unique(conditions)
	for(j in 1:length(bed$chr) ){
		Mtemp <- matrix(NA, ncol=length(uniqueCond), nrow=length(uniqueCond))
		colnames(Mtemp) <- uniqueCond
		rownames(Mtemp) <- uniqueCond
		PVALS[[j]] <- as.data.frame(Mtemp)
		rm(Mtemp)
	}
	#print(PVALS)	
	
	#FPCA:
	bspl <- create.bspline.basis(rangeval=c(-flank,flank),nbasis=nbasis, norder=4)
	argvalsBS <- -flank:flank
	
	
	for(j in 1:length(bed$chr) ){
		tempMatrix <- matrix(0.0, ncol=1+2*flank, nrow= length(bigwigs) )
		for (m in 1:length(bigwigs)) {
			tempMatrix[m,] <- fdamatrix[[m]][j,]
		}
		fdaData <- Data2fd(y=t(tempMatrix), argvals= argvalsBS, basisobj=bspl)		
		#plot.fd(fdaData)	
		pc <- pca.fd(fdobj=fdaData, nharm = pcs, harmfdPar=fdPar(fdaData),centerfns = TRUE)
		#plot.fd(pc$harmonics)
		
		# Select the PC scores for the components which amount >= 'variation'
		#print(pc$varprop)
		CS <- cumsum(pc$varprop)
		requiredPCs <- which(CS >= variation)[1]
		
		# Hotelling's T2 test between the conditions provided
		for (o in 1:length(uniqueCond)){
			for (p in 1:length(uniqueCond)){
				
				Xid <- which(conditions == uniqueCond[o] )
				Yid <- which(conditions == uniqueCond[p] )
				
				PVALS[[j]][o,p]  <- HotellingsT2(X=as.matrix(pc$scores[Xid ,1:requiredPCs],ncol=requiredPCs),Y=as.matrix(pc$scores[Yid,1:requiredPCs],ncol=requiredPCs), test="chi")$p.value
			}
		}
		# report the P-values in a list of data.frames (PVALS)		
	}
			
	return(list(fdaprofiles=fdamatrix, p.values = PVALS ) ) # report P-values!

}     

