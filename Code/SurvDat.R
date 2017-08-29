#mortDat 	<- read.csv(file="../Data/PumaMortKitten.csv")
mortDat 	<- read.csv(file="../Data/PumaMort.csv")
#mortDat 	<- read.csv(file="../Data/PumaMort_noDisperse.csv")
#mortDat 	<- read.csv(file="../Data/PumaMort_JMF.csv")
#remYoung 	<- which(mortDat$InitialAge <= 18)
#mortDat	<- mortDat[-remYoung,]

survDat 	<- mortDat[,-c(1:6, 71:73)]
survival.mat <- as.matrix(survDat)
survival.mat <- survival.mat[,-c(1:3, 64)]
#survival.mat <- survival.mat[,-c(1:3)]
survival.mat <- survival.mat[-which(apply(survival.mat, 1, sum, na.rm=T) ==0),]

#survival.mat <- survival.mat[-unique(which(is.na(survival.mat), arr.ind=T)[,1]),]
#stop()
N <- dim(survival.mat)[1]
tmax <- dim(survival.mat)[2]

age.mat 	<- matrix(NA, N, tmax)
state.mat 	<- age.mat  #1 is kitten, 2 is subadult, 3 is adult
firstObs 	<- vector('numeric', N)
lastObs 	<- vector('numeric', N)

for(i in 1:dim(age.mat)[1]) {

	firstObs[i] <- which(survival.mat[i,] >= 1)[1]
	temp 		<- which(survival.mat[i,] >= 1)
	
	if(any(survival.mat[i,] == 2, na.rm=T)) {
		lastObs[i] <- which(survival.mat[i,] == 2)
		if(lastObs[i] < tmax) {
			survival.mat[i, (lastObs[i]+2):tmax] <- NA
		}
	} else { lastObs[i] <- temp[length(temp)] }

	age.mat[i,firstObs[i]:tmax] 	<- seq(from=mortDat$InitialAge[i], by=3, length.out=tmax-firstObs[i]+1)
	temp 							<- age.mat[i,firstObs[i]:tmax]
	state.mat[i, which(age.mat[i,] <= 6)] <- 1
	state.mat[i, which(age.mat[i,] > 6 & age.mat[i,] <= 18)] <- 2
	state.mat[i, which(age.mat[i,] > 18)] <- 3 
	#state.mat[i, which(age.mat[i,] <= 6)] <- 1
	#state.mat[i, which(age.mat[i,] > 6)]  <- 2
	#state.mat[i, which(age.mat[i,] > 18)] <- 3 
		
}

survival.mat[which(survival.mat == 2)] <- 0.0