`projectPop` <-
function (demogT, envtlT, sts=NULL, tS=100000, tT=10000)
	{
	##if !is.null(sts) then there's within and between-state variation
	##get habitat state vector: choose initial state at random
    totLength <- ((tS + (2 * tT)) + 1) 
    states <- numeric(totLength)
    ns <- dim(envtlT)[1]
    states[1] <- ceiling(runif(1, 0, 2))

  for (i in 2:totLength)
    {
    nxt <- runif(1)
    Acs <- apply(envtlT, 2, cumsum)
    states[i] <- which(nxt<Acs[,states[i-1]])[1]
    }

	nClass <- dim(demogT[[1]])[1]
	growth <- matrix(0,nrow=tS,ncol=1) 		#to receive the sequence of one-time growth rates
	uvecs  <- matrix(0,nrow=nClass,ncol=tS)		#to receive the sequence of population structures
	vvecs  <- matrix(0,nrow=nClass,ncol=tS)		#to receive the sequence of reproductive value vectors
	rates <- vector("list",tS)

	v1 <- rep(1,nClass)	##initialize population vector with uniform proportions in each age-class
	v2 <- v1 <- v1/sum(v1)
	mgrowth <- 0
	maxTime <- tS + (2*tT)
	trun <- tS + tT

	for (i in 1:trun)
		{
			
		specMat <- ifelse(is.null(sts), states[i], sample(which(sts==states[i]))[1])
		Tf <- as.matrix(demogT[[specMat]])
	
		v1 <- Tf %*% v1
		growth1 <- sum(v1)
		v1 <- v1/growth1
		mgrowth <- mgrowth + log(growth1)

		if( i > tT)      				#only store the ones after the transients
			{
			i1          <- (i-tT)
			uvecs[,i1]  <- v1			#store right vectors (population structure)
			growth[i1]  <- growth1
			rates[[i1]] <- Tf
			}
	
		j  <- (maxTime-(i+1))
		v2 <- v2 %*% Tf					#makes 'left' vectors (reproductive value)
		v2 <- v2/(sum(v2))
				
		if (i > tT)	        			#only store the ones after the transients
			{
			vvecs[,(j-tT)] <- v2
			}
		}
	
	popn<-list(growth=growth, popStruct=uvecs, repVal=vvecs, rates=rates, nState=dim(envtlT)[1], states=states[(tT+1):(tS+tT)])
	class(popn) <- "projectedPop"
	return(popn)
	}

