
narrowpeaks <- function
(inputReg, scoresInfo, lmin = 0, nbf = 50, rpenalty= 0, nderiv= 0, npcomp = 5, pv = 80, pmaxscor = 0.0, ms = 0){
 
  propfilter <- 100;
  if( (propfilter < 0 | propfilter > 100) | (pv < 0 | pv >100) | (pmaxscor < 0 | pmaxscor > 100)) { stop(" 'propfilter', 'pv' and 'pmaxscor' must be in range 0-100 ") }
  if (npcomp > nbf){ stop(" The number of principal components 'ncomp' to be calculted must be less or equal to number of bicubic spline basis functions 'nbf'") }

  inputReg <- as.data.frame(inputReg)
  colnames(inputReg) <- c("chr", "start", "end", "length", "strand", "posPeak","score")
  inputReg <- inputReg[,c(1,2,3,6,7,4)]
 
  # Filter out regions smaller than Lmin base pairs and extract candidates from binary files
  filteredlmin <- which( inputReg$length >= lmin )
  enriched_regions <- inputReg[ filteredlmin, ]
  enriched_regions_filt <- enriched_regions
  binding_profiles <- list( )
  binding_profiles <- getBindingProfiles( win = enriched_regions, scoreFiles = scoresInfo$filenames )

  # Add vector of zeros for shift-to-zero correction of PCA scores
  # It means zero score when null binding (in the wig input file the estimated density of reads in the region is zero)
  # We want to move the origin of the euclidean space from (0,0) to PCs for no binding
  binding_profiles[[length(binding_profiles) + 1 ]] <- rep( 0, lmin + 1 )
  Chr0 = enriched_regions$chr[1]
  RowName0 <- as.numeric( row.names(enriched_regions)[length(row.names(enriched_regions))] ) + 1
  enriched_regions[ nrow( enriched_regions ) + 1, ] <- data.frame(chr = Chr0 , start = 1, end = lmin+1, posPeak = 1, score = 0, length = lmin, row.names=RowName0  )

  # List of profiles 2 matrix of profiles centered
  bp2Matrix <- list2matrix( binding_profiles )
  binding_profiles_extended <- bp2Matrix$dataMatrix
  mask_profiles <- bp2Matrix$dataMask  

  # Np is the number of profiles
  Np <- nrow( binding_profiles_extended )

  # L is the extension (nucleotides) of the profiles (extended) loaded
  L <- ncol( binding_profiles_extended )

  # Functional PCA
  pca_profiles <- fpca( profiles=binding_profiles_extended, numberOfComponents=npcomp, K=nbf, nderiv=nderiv, lambda=rpenalty ) 
  varpropFactor <- pca_profiles$fpcaout$varprop

  if (sum(100*varpropFactor) < pv){ stop(" Variance 'pv'% cannot be achieved by the number of components 'npcomp' specified. Please adjust parameters (i) Increasing the number of components 'npcomp'' or (ii) decreasing the proportion of variance 'pv'") }

  # Matrix W ( fda )
  BasisMatrix <- getbasismatrix(evalarg=pca_profiles$rangevalues, basisobj=pca_profiles$basisfunct, nderiv=0) 
  Wmat <- matrix(nrow = nbf, ncol = nbf)
  for (i in 1:nbf){ for(j in 1:nbf){ Wmat[i,j] <- sum( BasisMatrix[,i]*BasisMatrix[,j] ) } }

  # Matrix of scores on the Principal Components for every shape curve
  zScores <- pca_profiles$fpcaout$scores      

  # Scores to be shifted
  shiftFactor <- pca_profiles$fpcaout$scores[ nrow(pca_profiles$fpcaout$scores), ]

  # Eigenvalues and harmonics (principal components)
  eigenvalues <- pca_profiles$fpcaout$values
  harmonics <- pca_profiles$fpcaout$harmonics

  # New matrix of scores on the principal components 
  zShifted <- matrix(nrow = Np-1, ncol = npcomp)
  
  # New scores (no-binding on wig file will correspond to ZERO PCA scores)
  zShifted <- t(t(zScores[1: Np-1,]) - shiftFactor)  #Last row is deleted
  rm(enriched_regions)
  enriched_regions <- enriched_regions_filt

  # Global scores
  gz <- rowSums( zShifted^2 )  
  Np2 <- round( propfilter * (Np-1)/100 )
  sort_gz <- sort( gz, decreasing = TRUE, index.return =TRUE )
  #Profiles to be analyzed subsequently
  profiles_for_analysis = sort_gz[[2]][1:Np2]
  enriched_regions2 <- enriched_regions[ profiles_for_analysis, ]
  enriched_regions2[,7] <- data.frame(PCScore= sort_gz[[1]][1:Np2])

  # Calculation of minimum number of components for a given 'PV' (proportion of variation)
  mean_shape <- eval.fd( evalarg=1:L, fdobj=pca_profiles$fpcaout$mean )   #plot(mean_shape, type='l')
  ranked_components <-sort( 100*pca_profiles$fpcaout$var, decreasing = TRUE, index.return = TRUE )
  cum_var = 0.0
  comps_used = c()
  for (i in 1:npcomp){
	if (cum_var <= pv){ 
		cum_var <- cum_var + ranked_components[[1]][i]
		comps_used = c(comps_used, ranked_components[[2]][i])
	}
  }
  
  npcomp2 <- length(comps_used)  #Required number of components
  # 'cum_var' stores the variance accounted by the components
    
  # Matrix with NPCOMP2 vector-valued functional components for minimum variance PV
  functional_components <- matrix( 0, nrow = npcomp2, ncol = L )  
  for (i in 1:npcomp2)  { functional_components[i,] <- eval.fd( 1:L, pca_profiles$fpcaout$harmonics[comps_used[i]]) }

  # Functional score Function (measures global variation within ChIP-Seq wig signal values)
  functional_score <- matrix(0, nrow=Np2, ncol=L)
  for (i in 1:Np2) {
	a_s <- eval.fd(evalarg=1:L, fdobj=pca_profiles$fdaData[profiles_for_analysis[i]]) #approx. shape (" vector-valued function ")
	for (j in 1:npcomp2){
		functional_score[i,] <- functional_score[i,] + (functional_components[j,]*(a_s - mean_shape))^2  
	}

  }

  # signalValue is the average enrichment for the region
  meanSignal <- c()
  for (i in 1:(Np-1)) { meanSignal[i] <- mean(binding_profiles[[i]]) }

  # Selected candidate regions from the input wiggle track
  broadPeaks <- data.frame(n=1:(Np-1),chrom=enriched_regions[,1],length=enriched_regions[,6], siteStart=enriched_regions[,2], siteEnd=enriched_regions[,3],  max=enriched_regions[,5], average=round(meanSignal,2), fpcaScore=round(gz,2))

  # Thresholding loop, identify which candidate regions present variation contribution 
  # ...higher than PMAXSCOR. Then the region above the threshold is selected
  pmaxscor <- (pmaxscor/100) * max(functional_score,na.rm=TRUE)
  sites <- matrix(0,Np2,L)  # 'sites' indicate the genomic position above the cutoff 'PMAXSCOR'
  for (i in 1:Np2){
 	funcVar <- functional_score[i,]
 	for (j in 1:L) {
		if (is.na(funcVar[j])) { funcVar[j] = 0.0 }  
 		if (funcVar[j] >= pmaxscor) { sites[i,j]=1 }
 	}
  }

  # Number of candidates found above threshold Tsf
  sitesNumber <- 0  
  for (i in 1:Np2){   if (max(sites[i,]) != 0) sitesNumber = sitesNumber+1}
  rgz <- c() 
  for (i in 1:(Np-1)){ rgz[i] = which(sort_gz[[2]] == i) }


  bindingSites <- sitesCreator(s=sites*(mask_profiles[profiles_for_analysis,]), perfiles = binding_profiles_extended , componentes = functional_components, media=mean_shape,  enriched_regions=enriched_regions2,gz=gz,rgz=rgz, IDs=sort_gz[[2]][1:Np2], correction=shiftFactor, fdaprofiles=pca_profiles$fdaData) 
  bindingSites <- bindingSites[order(bindingSites$chrom ,bindingSites$siteStart),]
  bindingSites$broadPeak.subpeak <- as.character(bindingSites$broadPeak.subpeak)
  bindingSites$narrowedDownTo <- as.character(bindingSites$narrowedDownTo)
  candidates_position <- data.frame(start = as.numeric(enriched_regions[,2]), end = as.numeric(enriched_regions[,3]), length = as.numeric(enriched_regions[,6]) )

  if (nrow(bindingSites)>1){ narrowpeaksfinal <- mergeSites(filebs=bindingSites,profiles_for_analysis=profiles_for_analysis,candidates_position=candidates_position,iniLength=broadPeaks$siteEnd-broadPeaks$siteStart+1,N=ms) }
  else {
    narrowpeaksfinal <- bindingSites
    narrowpeaksfinal$merged <- FALSE
  }

  narrowpeaksfinal <- narrowpeaksfinal[order(narrowpeaksfinal$chrom ,narrowpeaksfinal$siteStart),]
  narrowpeaksfinal <- narrowpeaksfinal[ ,c(1:4,6:9)]

  # Genomic locations in GenomicRanges class GRanges
  broadPeaksGR <- GRanges(seqnames=Rle(as.character(broadPeaks$chrom)),ranges=(IRanges(broadPeaks$siteStart, end=broadPeaks$siteEnd)),max=as.integer(broadPeaks$max), average=broadPeaks$average, fpcaScore=broadPeaks$fpcaScore )
  narrowpeaksfinalGR <- GRanges(seqnames=Rle(as.character(narrowpeaksfinal$chrom)),ranges=(IRanges(narrowpeaksfinal$siteStart, end=narrowpeaksfinal$siteEnd)),broadPeak.subpeak=as.character(narrowpeaksfinal$broadPeak.subpeak), trimmedScore=narrowpeaksfinal$trimmedScore, narrowedDownTo=narrowpeaksfinal$narrowedDownTo, merged=as.logical(narrowpeaksfinal$merged) )
  seqlengths(broadPeaksGR) <- scoresInfo$chrL
  seqlengths(narrowpeaksfinalGR) <- scoresInfo$chrL

  return(list(fdaprofiles=pca_profiles$fdaData, broadPeaks = broadPeaksGR, narrowPeaks = narrowpeaksfinalGR, reqcomp = length(comps_used), pvar=cum_var  ))

}
