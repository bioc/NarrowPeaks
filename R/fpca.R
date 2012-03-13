
fpca <- function(profiles,numberOfComponents,K,nderiv,lambda){
  # Functional Principal Component Analysis of the read-enriched regions
#require(splines)

	numberOfProfiles <- length(profiles)
	profMatrix <- t(profiles)
	argvalsBS <- 1:nrow(profMatrix)
	bspl <- create.bspline.basis(rangeval=c(1,nrow(profMatrix)),nbasis=K, norder=4)
	fdamatrix <- Data2fd(y=profMatrix, argvals= argvalsBS, basisobj=bspl,nderiv=nderiv, lambda=lambda)
	fpcamatrix <- pca.fd(fdobj = fdamatrix, nharm = numberOfComponents, harmfdPar=fdPar(fdamatrix), centerfns = TRUE)
	pcOfprofMatrix <- eval.fd(argvalsBS,fpcamatrix$harmonics)
	return(list(eigenfunctions=pcOfprofMatrix,fpcaout=fpcamatrix,fdaData=fdamatrix,rangevalues=argvalsBS,basisfunct=bspl))

}
