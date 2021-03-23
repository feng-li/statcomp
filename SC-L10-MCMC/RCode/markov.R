#Binary Markov chain
rm(list=ls())
#Legth of chain
n<-100
#Initial value
x0<-0
#Vector to store values
x<-rep(NA,n)
#Set first value
if (x0==0){x[1]<-rbinom(1,1,0.2)}
if (x0==1){x[1]<-rbinom(1,1,0.7)}
#We use rbinom(1,1,p) to simulate
#p is the  probability X=1

#Put this in a loop
for  (i in 2:n){
  if (x[i-1]==0){x[i]<-rbinom(1,1,0.2)}
  if (x[i-1]==1){x[i]<-rbinom(1,1,0.7)}
}
plot(1:n,x,"l",main='Sequence of X',xlab = 'i')


#Transition Matrix
P<-matrix(c(0.8,0.2,0.3,0.7),2,2,byrow = TRUE)
#P<-matrix(c(0.9,0.1,0.15,0.85),2,2,byrow = TRUE)
#Set pi_1
pi<-c(0.8,0.2)
#pi<-c(0.15,0.85)
#Set n
n<-100

#Repeat Recursion
for (i in 1:(n-1)){
  pi<-pi%*%P
}

print(pi)

#Code to keep pi and plot it
#Initialise all values of pi in a matrix
pi<-matrix(NA,n,2)

#Set first value of pi
pi[1,]<-c(0.8,0.2)

#Recursion
for (i in 2:n){
  pi[i,]<-pi[i-1,]%*%P
}


#Two Plots
par(mfrow=c(2,1))
#First plot
plot(1:n,pi[,1],'l',main=expression(X[0]==0),xlab = 'i',ylab=expression(Pr(X[i]==0)))

#Repeat with different initial condition
pi[1,]<-c(0.3,0.7)
for (i in 2:n){
  pi[i,]<-pi[i-1,]%*%P
}
plot(1:n,pi[,1],'l',main=expression(X[0]==1),xlab = 'i',ylab=expression(Pr(X[i]==0)))

#Random Walk
x0<-0
n<-100
x<-rep(NA,n)
x[1]<-rnorm(1,x0,1)
for (i in 2:n){
  x[i]<-rnorm(1,x[i-1],1)
}
plot(1:n,x,"l")


#Metropolis Algorithm (Simple)
set.seed(3)
x0<-0
n<-110000
x<-rep(NA,n)

#Metropolis  Algorithm First Step
xold<-x0 #Set x old
xnew<-rnorm(1,xold,100) #Generate proposal
a<-min(1,dt(xnew,5)/dt(xold,5)) #Compute acceptance probability
accept<-rbinom(1,1,a) #Accept yes/no
if(accept){x[1]<-xnew} #If accept x1=xnew
if(!accept){x[1]<-xold} #Otherwise accept x1=xold
sd = 100
#Now inside loop
for (i in 2:n){
  xold<-x[i-1]
  xnew<-rnorm(1,xold,sd)
  a<-min(1,dt(xnew,5)/dt(xold,5))
  accept<-rbinom(1,1,a)
  if(accept){x[i]<-xnew}
  if(!accept){x[i]<-xold}
}
#pdf('../Pics/sp5.pdf')
plot(1:1000,x[1:1000],"l",ylab = 'x',xlab='',main=paste('Initial Value =0, Proposal S.D. =', sd))
#dev.off()


#Random Walk Metropolis (More complicated)
x0<-0
bi<-10000 #Burn In
n<-110000 #Sample Size
accrate<-0 #Counter for acceptance rate
x<-rep(NA,n) #Initialise vector to store x
nu<-2 #Proposal Variance

for (i in 1:n){
  xold<-ifelse((i==1),x0,x[i-1])
  xnew<-rnorm(1,xold,nu)
  a<-min(1,dt(xnew,5)/dt(xold,5))
  accept<-rbinom(1,1,a)
  if(accept){x[i]<-xnew}
  if(accept&&(i>bi)){accrate<-accrate+1}
  if(!accept){x[i]<-xold}
}
accrate<-accrate/(n-bi)
plot(1:n,x,"l")
library("coda")
geweke.diag(x[10001:110000])
effectiveSize(x[10001:110000])
