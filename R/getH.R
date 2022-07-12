`getH` <- function (sts) 
	{
	##from a vector of envtl states, and returns the envtl transition matrix
	newState<-as.numeric(as.factor(sts))
	HH<-matrix(0,nrow=max(newState),ncol=max(newState))

	for (i in 2:length(newState))
		{HH[newState[i],newState[i-1]]<-HH[newState[i],newState[i-1]]+1}
	for (i in 1:dim(HH)[1]) {HH[,i]<-HH[,i]/colSums(HH)[i]}
	return(HH)
	}

