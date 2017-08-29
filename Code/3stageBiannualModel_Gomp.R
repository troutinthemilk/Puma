library(runjags)
source('BiannualIntegrated_Gomp.jags')
source('JagsFuncs.R')

#read in survival data
source('SurvDat.R')
N 			<- dim(state.mat)[1]
tmax 		<- dim(state.mat)[2]

abund.dat 	<- read.csv(file="../Data/Puma3stage.csv")
Adult.dat 	<- read.csv(file="../Data/PumaAdult.csv", header=T)
Cov.dat 	<- read.csv(file="../Data/AbundancesUpdated.csv", header=T)
effort.dat 	<- read.csv(file="../Data/Effort.csv")[,2]
effort.scale <- effort.dat/50
tmax 		<- dim(state.mat)[2]

A1 <- data.frame(State=rep("Adult", length(Adult.dat[,4])), Season=rep("NoHunt", length(Adult.dat[,4])), Abundance=Adult.dat[,4]+1)
A2 <- data.frame(State=rep("Adult", length(Adult.dat[,3])), Season=rep("Hunt", length(Adult.dat[,4])), Abundance=Adult.dat[,3]+1)
S1 <- data.frame(State=rep("Subadult", length(abund.dat[,5])), Season=rep("NoHunt", length(Adult.dat[,4])), Abundance=abund.dat[,5]+1)
S2 <- data.frame(State=rep("Subadult", length(abund.dat[,4])), Season=rep("Hunt", length(Adult.dat[,4])), Abundance=abund.dat[,4]+1)
K1 <- data.frame(State=rep("Kitten", length(abund.dat[,3])), Season=rep("NoHunt", length(Adult.dat[,4])), Abundance=abund.dat[,3]+1)
K2 <- data.frame(State=rep("Kitten", length(abund.dat[,2])), Season=rep("Hunt", length(Adult.dat[,4])), Abundance=Adult.dat[,4]+1)

abund.frame <- rbind(A1, A2, S1, S2, K1, K2)
abund.frame <- cbind(abund.frame, effort=rep(effort.dat, 6))

library(AICcmodavg)
lm0 <- lm(log(Abundance) ~ 1, data=abund.frame)
lm1 <- lm(log(Abundance) ~ effort, data=abund.frame)
lm1.5 <- lm(log(Abundance) ~ effort + I(effort^2), data=abund.frame)
lm2 <- lm(log(Abundance) ~ effort+State, data=abund.frame)
lm2.5 <- lm(log(Abundance) ~ State, data=abund.frame)
lm3 <- lm(log(Abundance) ~ effort+Season, data=abund.frame)
lm3.5 <- lm(log(Abundance) ~ effort+Season+State, data=abund.frame)
lm4 <- lm(log(Abundance) ~ Season*State + effort, data=abund.frame)
print(summary(lm3.5)	)
print(aictab(list(lm0, lm1, lm1.5, lm2, lm2.5, lm3, lm3.5, lm4)))


par(mfrow=c(1,2))
library(RColorBrewer)
col.vec <- brewer.pal(3, "Set1")
plot(effort.dat, Adult.dat[,4], pch=19, col=col.vec[1], cex=1.3, xlab="Effort", ylab="Abundance", ylim=c(0, 20))
points(effort.dat, Adult.dat[,3], pch=18, col=col.vec[1], cex=1.3)
points(effort.dat, abund.dat[,3], pch=19, col=col.vec[2], cex=1.3)
points(effort.dat, abund.dat[,2], pch=18, col=col.vec[2], cex=1.3)
points(effort.dat, abund.dat[,5], pch=19, col=col.vec[3], cex=1.3)
points(effort.dat, abund.dat[,4], pch=18, col=col.vec[3], cex=1.3)

#b <- coef(lm3.5)
b <- coef(lm3.5)
b[3:5] <- 0
effort.vec <- seq(0, 50, length.out=100)

lines(effort.vec, exp(b[1] + b[2]*effort.vec), lwd=2, col=col.vec[1])
lines(effort.vec, exp(b[1] + b[2]*effort.vec+b[4]), lwd=2, col=col.vec[2])
lines(effort.vec, exp(b[1] + b[2]*effort.vec+b[5]), lwd=2, col=col.vec[3])


legend('topleft', legend=c("Adult", "Subadult", "Kitten", "Hunt", "No hunt"), col=c(col.vec, 'black', 'black'), pch=c(19,19,19, 18, 19))

library(lme4)

z <- matrix(NA, N, tmax)
for(i in 1:N) {
	z[i, 1:(firstObs[i])] <- NA
	z[i, (firstObs[i]):lastObs[i]] <- 1
}

##read in abundance data
PumaHunt		<- (Adult.dat$Hunt+1)
PumaNoHunt		<- (Adult.dat$NoHunt+1)
AdultTot		<- c(matrix(c(PumaNoHunt, PumaHunt), 2, byrow = T)) 
KittenHunt		<- (abund.dat$Kittens.Hunt)
KittenNoHunt	<- (abund.dat$Kittens.NoHunt)
KittenTot		<- c(matrix(c(KittenNoHunt, KittenHunt), 2, byrow = T)) 
SubHunt			<- (abund.dat$Subadult.Hunt) - KittenHunt + 1
SubNoHunt		<- (abund.dat$Subadult.NoHunt) - KittenNoHunt +1
SubTot			<- c(matrix(c(SubNoHunt, SubHunt), 2, byrow = T)) 
KittenHunt		<- (abund.dat$Kittens.Hunt)+1
KittenNoHunt	<- (abund.dat$Kittens.NoHunt)+1

plot(PumaHunt, pch=19, ylim=c(0,20))
lines(Adult.dat$Hunt*exp(-b[2]*effort.dat)+1)
points(PumaNoHunt, pch=19, col='red')
lines(Adult.dat$NoHunt*exp(-b[2]*effort.dat), col='red')
lines(effort.dat/max(effort.dat), col=rgb(0,0,0,0.1), lwd=2)

Cov 			<- abund.dat$Wolf
nYears 			<- length(Cov)

Cov.surv <- YearVec	<- NULL
for(i in 1:length(Cov)) { 
	YearVec  <- c(YearVec, rep(i,4))
	Cov.surv <- c(Cov.surv, rep(Cov[i], 4))
}

#prep data for JAGS
inits 		<- list(z=z, a=0, b=0)
hunt.state	<- rep(c(1,1,2,2), length.out=dim(state.mat)[2]) #is 1 or 2 hunting?
data.list 	<- list(nObs=tmax, nYears=nYears, nMarked=N, firstObs=firstObs, CapHis=survival.mat, stateMat=state.mat, huntState=hunt.state, YearVec=YearVec, PumaHunt=log((PumaHunt)*exp(-b[2]*effort.dat - b[3])-1), PumaNoHunt=log((PumaNoHunt)*exp(-b[2]*effort.dat)-1), KittenHunt=log((KittenHunt)*exp(-b[2]*effort.dat- b[3]-b[5])-1), KittenNoHunt=log((KittenNoHunt)*exp(-b[2]*effort.dat-b[5])-1), effort=effort.dat, SubHunt=log((SubHunt)*exp(-b[2]*effort.dat-b[3]-b[4])-1), SubNoHunt=log((SubNoHunt)*exp(-b[2]*effort.dat-b[4])-1), Cov=scale(Cov.dat$Wolves)[,1], PumaCov=scale(Adult.dat$Total)[,1])

params <- c('f0', 'fA', 'fC', 'aK1', 'aK2', 'aS1', 'aS2', 'aA1', 'aA2', 'sAC', 'sSC', 'sKC', 'AS', 'SS', 'KS', 'S1Htrue[1]', 'S1NHtrue[1]', 'S2Htrue[1]', 'S2NHtrue[1]', 'tauK', 'tauS', 'tauA', 'adetect', 'deviance')
par(mfrow=c(1,1))

	jags.Null <- run.jags(fullIntBiannualMod0, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.Null$DIC 	<- DIC.extract(jags.Null)
print(jags.Null$DIC)

	jags.DD1 <- run.jags(fullIntBiannualMod1, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.DD1$DIC 	<- DIC.extract(jags.DD1)

	jags.Wolf <- run.jags(fullIntBiannualMod2, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.Wolf$DIC 	<- DIC.extract(jags.Wolf)

	jags.fullWolf <- run.jags(fullIntBiannualMod3, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.fullWolf$DIC 	<- DIC.extract(jags.fullWolf)

	data.list 	<- list(nObs=tmax, nYears=nYears, nMarked=N, firstObs=firstObs, CapHis=survival.mat, stateMat=state.mat, huntState=hunt.state, YearVec=YearVec, PumaHunt=log((PumaHunt)*exp(-b[2]*effort.dat - b[3])-1), PumaNoHunt=log((PumaNoHunt)*exp(-b[2]*effort.dat)-1), KittenHunt=log((KittenHunt)*exp(-b[2]*effort.dat- b[3]-b[5])-1), KittenNoHunt=log((KittenNoHunt)*exp(-b[2]*effort.dat-b[5])-1), effort=effort.dat, SubHunt=log((SubHunt)*exp(-b[2]*effort.dat-b[3]-b[4])-1), SubNoHunt=log((SubNoHunt)*exp(-b[2]*effort.dat-b[4])-1), Cov=scale(Cov.dat$Elk.density)[,1], PumaCov=scale(Adult.dat$Total)[,1])

	jags.Elk <- run.jags(fullIntBiannualMod2, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.Elk$DIC 	<- DIC.extract(jags.Elk)
	
	jags.fullElk <- run.jags(fullIntBiannualMod3, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.fullElk$DIC 	<- DIC.extract(jags.fullElk)
	
	
	data.list 	<- list(nObs=tmax, nYears=nYears, nMarked=N, firstObs=firstObs, CapHis=survival.mat, stateMat=state.mat, huntState=hunt.state, YearVec=YearVec, PumaHunt=log((PumaHunt)*exp(-b[2]*effort.dat - b[3])-1), PumaNoHunt=log((PumaNoHunt)*exp(-b[2]*effort.dat)-1), KittenHunt=log((KittenHunt)*exp(-b[2]*effort.dat- b[3]-b[5])-1), KittenNoHunt=log((KittenNoHunt)*exp(-b[2]*effort.dat-b[5])-1), effort=effort.dat, SubHunt=log((SubHunt)*exp(-b[2]*effort.dat-b[3]-b[4])-1), SubNoHunt=log((SubNoHunt)*exp(-b[2]*effort.dat-b[4])-1), Cov=scale(Cov.dat$Elk.off)[,1], PumaCov=scale(Adult.dat$Total)[,1])
	
	jags.fullElkoff <- run.jags(fullIntBiannualMod3, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.fullElkoff$DIC 	<- DIC.extract(jags.fullElkoff)

	jags.postHoc <- run.jags(fullIntBiannualMod4, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)
jags.postHoc$DIC 	<- DIC.extract(jags.postHoc)

	data.list 	<- list(nObs=tmax, nYears=nYears, nMarked=N, firstObs=firstObs, CapHis=survival.mat, stateMat=state.mat, huntState=hunt.state, YearVec=YearVec, PumaHunt=log((PumaHunt)*exp(-b[2]*effort.dat - b[3])-1), PumaNoHunt=log((PumaNoHunt)*exp(-b[2]*effort.dat)-1), KittenHunt=log((KittenHunt)*exp(-b[2]*effort.dat- b[3]-b[5])-1), KittenNoHunt=log((KittenNoHunt)*exp(-b[2]*effort.dat-b[5])-1), effort=effort.dat, SubHunt=log((SubHunt)*exp(-b[2]*effort.dat-b[3]-b[4])-1), SubNoHunt=log((SubNoHunt)*exp(-b[2]*effort.dat-b[4])-1), Elk=scale(Cov.dat$Elk.off)[,1], Wolves=scale(Cov.dat$Wolves)[,1], PumaCov=scale(Adult.dat$Total)[,1])

	params <- c('f0', 'fA', 'fC', 'aK1', 'aK2', 'aS1', 'aS2', 'aA1', 'aA2', 'sKC', 'sACW', 'sACE', 'AS', 'SS', 'KS', 'adetect', 'deviance')
	jags.apriori <- run.jags(fullIntBiannualModapriori, burnin=1e4, sample=1e4, n.chains=4, thin=10, adapt=1e4, monitor=params, data=data.list, inits=inits, method='parallel', silent.jags=T)

jags.apriori$DIC 	<- DIC.extract(jags.apriori)


	DIC.vec <- c("No DD"=jags.Null$DIC[3], "DD1"=jags.DD1$DIC[3], "Wolf"=jags.Wolf$DIC[3], "Elk"=jags.Elk$DIC[3], "Puma + Wolf"=jags.fullWolf$DIC[3], "Puma + Elk"=jags.fullElk$DIC[3], "Puma + Elk.on" =jags.fullElkon$DIC[3], "Puma + Elk.off"=jags.fullElkoff$DIC[3], 'apriori'=jags.apriori$DIC[3], 'postHoc'=jags.postHoc$DIC[3])

	print(DIC.vec-min(DIC.vec))
		
	save.image('PumaModelFits.Rdata')

par(mfrow=c(1,3))
MCMCplot(jags.DD1$mcmc, params=c("fA", "KS", "SS", "AS"), labels=c("Density-dependent fecundity", "Density-dependent kitten survival", "Density-dependent subadult survival", "Density-dependent adult survival"))
MCMCplot(jags.Wolf$mcmc, params=c("fC", "sKC", "sSC", "sAC"), labels=c("Wolves on fecundity", "Wolves on kitten survival", "Wolves on subadult survival", "Wolves on adult survival"))
MCMCplot(jags.Elk$mcmc, params=c("fC", "sKC", "sSC", "sAC"), labels=c("Elk on fecundity", "Elk on kitten survival", "Elk on subadult survival", "Elkon adult survival"))

par(mfrow=c(1,2))
MCMCplot(jags.fullWolf$mcmc, params=c("fA", "fC", "KS", "SS", "AS", "sKC", "sSC", "sAC"), labels=c("Puma on fecundity", "Wolves on fecundity", "Puma on kitten survival", "Puma on subadult survival", "Puma on adult survival", "Wolves on kitten survival", "Wolves on subadult survival", "Wolves on adult survival"), main="Model with wolf covariate", xlim=c(-5,5))
MCMCplot(jags.fullElk$mcmc, params=c("fA", "fC", "KS", "SS", "AS", "sKC", "sSC", "sAC"), labels=
	c("Puma on fecundity", "Elk on fecundity", "Puma on kitten survival", "Puma on subadult survival", "Puma on adult survival", "Elk on kitten survival", "Elk on subadult survival", "Elk on adult survival"), main="Model with elk covariate", xlim=c(-5,5))
