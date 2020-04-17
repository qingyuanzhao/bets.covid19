library(coda)

for (alg in c("simple_alpha_10.000","simple_alpha_1.000","twoR_alpha_10.000", "twoR_alpha_1.000")) {
  d = read.csv(paste(paste("simulations_",alg,sep=""),".csv",sep=""),check.names = F) 
  names(d) = c("r","w_mean","copy")
  x = list()
  for (copy in 1:8) {
    this.sampler = (d$copy==copy)
    x=c(x,list(mcmc(d[this.sampler,c("w_mean","r")])))
  }
  x = mcmc.list(x)
  print(paste("Gelman-Rubin diagnostic for simulation: ",alg,sep=""))
  print(gelman.diag(x, confidence = 0.95, transform=F, autoburnin=T, multivariate=F))
}


