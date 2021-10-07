rm(list=ls())
library(coda)
#Function

target<-function(mu){
  #These can be defined inside or outside the function
  n<-10
  sumy<-17.34
  sumy2<-32
  eta<-1.8
  tau2<-0.25
  
  t1<-(sumy2-2*mu*sumy+n*mu^2)^(-n/2)
  t2<-exp(-((mu-eta)^2)/(2*tau2))
  t<-t1*t2
  return(t)
}


#Random Walk Metropolis (More complicated)
mu0<-1.8
bi<-5000 #Burn In
M<-25000 #Sample Size
accrate<-0 #Counter for acceptance rate
mu<-rep(NA,M) #Initialise vector to store x
nu2<-0.02 #Proposal Variance

for (j in 1:M){
  muold<-ifelse((j==1),mu0,mu[j-1]) #Set muold
  munew<-rnorm(1,muold,sqrt(nu2)) #Propose a value mu new.  ATTENTION: sqrt(nu2)
  a<-min(1,target(munew)/target(muold)) #Compute acceptance probability
  accept<-rbinom(1,1,a) #Accept YES/NO
  if(accept){mu[j]<-munew} #If accept set munew
  if(accept&&(j>bi)){accrate<-accrate+1} #Also if accept add to 
  if(!accept){mu[j]<-muold}
}

accrate<-accrate/(M-bi)
plot(1:M,mu,"l",xlab='j',ylab=expression(mu))

#Now throw away burn in (nobi='no burn in')
mu_nobi<-mu[(bi+1):M]

#Posterior Mean
mean(mu_nobi)

#Posterior Median
median(mu_nobi)

#Credible interval
quantile(mu_nobi,0.025) #2.5% percentile
quantile(mu_nobi,0.975) #97.5% percentile

#Probability that mu>1.7

sum(mu_nobi>1.7)/(M-bi)  #The numerator counts the number of elements where mu_nobi>1.7 is true


#Effective Sample Size
effectiveSize(mu_nobi)

