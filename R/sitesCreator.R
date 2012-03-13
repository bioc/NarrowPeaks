
sitesCreator <- function(s, perfiles, componentes, media, enriched_regions, gz,rgz, IDs, correction, fdaprofiles){

  Nrows = nrow(s)
  Ncols = ncol(s)
  Nbs = 0  
  chromosome <- c()
  start_bs <- c()
  end_bs <- c()
  prof_id <- c()
  bs_id <- c()
  length_bs <- c()
  startBS <-c()
  scre <-c()
  reducted <-c()
  
  for (i in 1:Nrows){
    bs <- s[i,]

    if( max(bs) == 1 ){
      bs <- c(1000,bs,1000)
      profile_id <- IDs[i]
      chr <- as.character(enriched_regions[i,1])
      lengthEnr <- as.numeric(enriched_regions[i,6])
      over_bp <- Ncols - lengthEnr +1 
      genome_position <- seq(from= floor(as.numeric(enriched_regions[i,2])-over_bp/2), to=floor(as.numeric(enriched_regions[i,3])+over_bp/2)-1 , by=1)

      if (over_bp %% 2 == 1){
	genome_position <- seq(from= floor(as.numeric(enriched_regions[i,2])-over_bp/2)+1, to=floor(as.numeric(enriched_regions[i,3])+over_bp/2) , by=1)
      }

      bs_index <- which(bs==1)
      bsID <- 0
      for(j in 1:length(bs_index)){

	  if (bs[bs_index[j]-1]==0 || bs[bs_index[j]-1]==1000)  {
		  
		      startBS <- genome_position[bs_index[j]-1]  
		      startFiltro = bs_index[j]-1
	  }

	  if (bs[bs_index[j]+1]==0 || bs[bs_index[j]+1]==1000) {
	      endFiltro = bs_index[j]-1
	      Nbs = Nbs+1
	      bsID = bsID +1
	      endBS <- genome_position[bs_index[j]-1]  
	      chromosome[Nbs] <- chr
	      start_bs[Nbs] <- startBS
	      end_bs[Nbs] <- endBS
	      prof_id[Nbs] <- profile_id
	      bs_id[Nbs] <- paste(profile_id,bsID,sep=".")
	      length_bs[Nbs] <- endBS-startBS+1
	      reducted[Nbs] <- 100*( length_bs[Nbs]/(lengthEnr + 0) ) 
	      filtroSubpeak <- rep(0,Ncols)
	      filtroSubpeak[startFiltro:endFiltro] <-1
	      scre[Nbs] <- 0.0
	      scoreK <- 0

	      for (k in 1:nrow(componentes)) { 
		  scoreK <- sum( ( filtroSubpeak * eval.fd(1:Ncols,fdaprofiles[profile_id]) -media)*componentes[k,] )
		  scre[Nbs] <- scre[Nbs] + ( scoreK -correction[k] )^2 }  
	      }
	  }
      }
  }

return(data.frame(chrom=chromosome, length= length_bs,siteStart=start_bs,siteEnd=end_bs, broadPeak=prof_id, broadPeak.subpeak=bs_id, trimmedScore = round(sqrt(scre),2), narrowedDownTo = paste(round(reducted,2),'%',sep=''))) 

}
