rm(list=ls())
n<-100 #Set sample size
p<-3 # Set number of regressors

#### Steps to randomly generate data from a Logistic Regression model 

#Randomly generate X from a U(0,1) distrbution, and put in an nxp matrix
X<-matrix(runif(n*p),n,p) 

#Set true values of the coefficients
betatrue=as.vector(c(-1,0,1)) 

#Compute vector Xb
Xbetatrue<-X%*%betatrue 

#Compute probabilities Y_i=1
prob<-exp(Xbetatrue)/(1+exp(Xbetatrue)) 

#Geneate y by making n draws from the uniform distribution
u<-runif(n) 
y<-as.numeric(u<prob) 


#Log Likelihood Function
logliklogit<-function(beta){
  Xb<-X%*%beta
  log1Xb<-log(1+exp(Xb))
  loglik<-(t(y)%*%Xb)-sum(log1Xb)
  return(-loglik)
}

gradlogit<-function(beta){
  Xb<-X%*%beta
  lamXb<-(exp(Xb))/(exp(Xb)+1)
  grad<-t(X)%*%(y-lamXb)
  return(-grad)
}

hesslogit<-function(beta){
  Xb<-X%*%beta
  weight<-(1+exp(Xb))^(-2)
  weightm<-diag(as.vector(weight))
  hess<--t(X)%*%weightm%*%X
  return(-hess)
}

#Initialise beta0 at least squares estimate
beta0<-solve(t(X)%*%X,t(X)%*%y)


#Set Tolerance
epsilon<-1E-10
#Set Maximum Iterations
MaxIter<-500
#Initialise n
n<-1 
#Initialise xnew
betanew<-beta0-solve(hesslogit(beta0),gradlogit(beta0))
#Initialise change in function value
change<-abs(logliklogit(betanew)-logliklogit(beta0)) 
#While Loop
while((n<=MaxIter)&&(change>epsilon)){ #NOTE: while loop breaks if either MaxIter reached of change in function value is less than tolerance
  #Update n  
  n<-n+1
  
  #Print Iteration Number
  print(paste('Iteration:',n))
  
  #New becomes old
  betanew->betaold 
  #xnew is updated
  betanew<-betaold-solve(hesslogit(betaold),gradlogit(betaold)) 
  
  #Print betaold and betanew at each step
  print(('Betaold:'))
  print((betaold))
  print(('Betanew:'))
  print((betanew))
  
  #Change in function value updated
  change<-abs(logliklogit(betanew)-logliklogit(betaold)) 
  
  #Print Change at each step
  print(paste('Change:',change))
  
}

betahat1<-betanew

betaoptim<-optim(beta0,fn = logliklogit,gr = gradlogit,hessian=TRUE,method='L-BFGS-B')

betahat2<-betaoptim$par
















