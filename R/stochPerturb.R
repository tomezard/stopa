`stochPerturb` <-
function (simltn)
	{
	if (class(simltn) != "projectedPop") {stop("Object class is not projectedPop.")}
	
	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	#CALCULATE GROWTH RATES AND AVERAGE QUANTITIES FROM THE PROJECTION
	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	nClass <- dim(simltn$rates[[1]])[1]
	avmat  <- matrix(rowMeans(sapply(simltn$rates, c)), nClass)
	lambda0 <- eigen.analysis(avmat)$lambda1
	a <- mean(log(simltn$growth))
	
	statemats <- vector("list", simltn$nState)
	stLambda0 <- rep(0, simltn$nState)
	usedStates <- as.numeric(names(tapply(simltn$states,simltn$states,length)))
	for (i in 1:length(usedStates)) 
		{
		statemats[[usedStates[i]]] <- matrix(rowMeans(sapply(simltn$rates[simltn$states==usedStates[i]], c)), nClass)
		stLambda0[i] <- Re(eigen(statemats[[usedStates[i]]])$values[1])
		}
	averageStuff <- list(a=a, lambdaS=exp(a), r=log(lambda0), lambda0=lambda0, avmat=avmat, statemats=statemats, stLambda0=stLambda0)

	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	#INITIALIZATION FOR PERTURBATION ANALYSIS, FOR DESCRIPTIONS OF EACH SEE HELP FILE
	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	aa<-dim(avmat)[1]
	sensE <- sensA <- sensV <- sensW <- sensB <- matrix(0,nrow=aa,ncol=aa)
	sensHS<-vector("list",simltn$nState)
	for(i in 1:simltn$nState) {sensHS[[i]]<-matrix(0,aa,aa)}

	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	#PERTURB ALONG THE SIMULATION
	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	ts<-length(simltn$growth)
	for (t in 2:ts) 
		{ 	
		#only store relevant information, therefore start at 2, as dependent upon conditions at time t-1
		Tf <- simltn$rates[[t]]

		mat2S <- Tf
		mat2A <- avmat
		mat2V <- Tf - avmat

		stateav <- statemats[[simltn$states[t]]]
		mat2W <- Tf - stateav
		mat2B <- stateav - avmat
		
		mat2S<-diag(as.numeric(simltn$repVal[,t])) %*% mat2S %*% diag(as.numeric(simltn$popStruct[,(t-1)]))
		mat2A<-diag(as.numeric(simltn$repVal[,t])) %*% mat2A %*% diag(as.numeric(simltn$popStruct[,(t-1)]))
		mat2V<-diag(as.numeric(simltn$repVal[,t])) %*% mat2V %*% diag(as.numeric(simltn$popStruct[,(t-1)]))
		mat2W<-diag(as.numeric(simltn$repVal[,t])) %*% mat2W %*% diag(as.numeric(simltn$popStruct[,(t-1)]))
		mat2B<-diag(as.numeric(simltn$repVal[,t])) %*% mat2B %*% diag(as.numeric(simltn$popStruct[,(t-1)]))

		scale1<-drop(as.numeric(simltn$repVal[,t]) %*% as.numeric(simltn$popStruct[,t]))
		
		mat2S<-mat2S/(scale1*simltn$growth[t])	
		sensE<-sensE + mat2S; 					
		
		mat2A<-mat2A/(scale1*simltn$growth[t])	
		sensA<-sensA + mat2A; 					

		mat2V<-mat2V/(scale1*simltn$growth[t]) 	
		sensV<-sensV + mat2V; 					

		mat2W<-mat2W/(scale1*simltn$growth[t]) 	
		sensW<-sensW + mat2W; 					

		mat2B<-mat2B/(scale1*simltn$growth[t]) 	
		sensB<-sensB + mat2B;

		sensHS[[simltn$states[t]]]<-sensHS[[simltn$states[t]]]+mat2S 
		}

	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	#RESCALING AND "CONVIENIENT CHECKS ON CALCULATION"
	#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
	sensE<-sensE/(ts-1);
	sensA<-sensA/(ts-1);
	sensV<-sensV/(ts-1);
	sensW<-sensW/(ts-1);
	sensB<-sensB/(ts-1);
	for(i in 1:length(sensHS)) {sensHS[[i]]<-sensHS[[i]]/(ts-1)}

	if (round(sum(sensE),12) != 1) {warning("sum of eS is not 1:\n",paste("sum(eS) = ",sum(sensE),sep=""))}
	hsSum <- sum(sapply(sensHS, sum))
	#if (round(hsSum != 1) {warning("Habitat-stage elasticities do not sum to 1:\n",paste("sum(hsSum) = ",sum(hsSum),sep=""))}
	if (round(sum(sensA+sensV-sensE),12) != 0) 
			{warning("eSmu and eSsigma do not sum to eS:\n",paste("sum(eS-(eSmu+eSsigma)) = ",sum(sensA+sensV-sensE),sep=""))}
	if (round(sum(sensV-(sensB+sensW)),12) != 0) 
			{warning("eSigmaW and eSigmaB do not sum to eSsigma:\n",paste("sum(eSsigma-(eSigmaB+eSigmaW)) = ",sum(sensV-(sensB+sensW)),sep=""))}

	elastcts<-list(eD=eigen.analysis(avmat)$elasticities,eS=sensE,eSmu=sensA, eSsigma=sensV, eSigmaW=sensW, eSigmaB=sensB, eSBeta=sensHS)
	return(elastcts)
	}
