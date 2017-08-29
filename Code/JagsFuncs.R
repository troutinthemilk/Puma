logit <- function(x) { exp(x)/(1 + exp(x)) }

DIC.extract <- function(jags.out) {
	if(is.numeric(jags.out$DIC[1])) {return(jags.out$DIC)}

	dev <- extract(jags.out, 'mean.deviance', force.resample=F, silent=TRUE, progress.bar="none")
	pd.temp  <- extract(jags.out, 'mean.pd', force.resample=F, silent=TRUE, progress.bar="none")
	print(sum(is.na(pd.temp)))	
	pd <- sum(pd.temp, na.rm=T)

	return(c(Deviance=sum(dev), Pd=pd, DIC=sum(dev) + pd))
}

AIC.extract <- function(jags.out) {
	if(is.numeric(jags.out$DIC[1])) {return(jags.out$DIC)}

	dev <- extract(jags.out, 'mean.deviance', silent=TRUE, progress.bar="none")
	pd.temp  <- extract(jags.out, 'mean.pd', silent=TRUE, progress.bar="none")
	print(sum(is.na(pd.temp)))	
	pd <- sum(pd.temp, na.rm=T)

	return(c(Deviance=sum(dev), Pd=pd, DIC=sum(dev) + pd))
}