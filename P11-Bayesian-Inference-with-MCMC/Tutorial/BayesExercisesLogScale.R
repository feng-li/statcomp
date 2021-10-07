#Aim is to find x that maximises this function
fun<-function(x){exp(-1000*x^2)}

#Let's evaluate x at a few points
a<-fun(0.5)
b<-fun(1)
c<-fun(1.01)

#Let's look at the ratio of two of these points?
ratio<-fun(1)/fun(1.01)

#b and c are both 0 and d is NaN.  What is going on?

#Let's try to optimise
oo<-optim(5,fun,lower=-10,upper=10,method='L-BFGS-B',control=list(fnscale=-1))
oo$par
#The maximum is not 5! What's going on?

##Let's try to use the log of this function
#Aim is find x that maximises this function
logfun<-function(x){(-1000*x^2)}

#Let's evaluate x at a few points
loga<-logfun(0.5)
logb<-logfun(1)
logc<-logfun(1.01)

#Let's look at the ratio of two of these points?
logratio<-logfun(1)-logfun(1.01)
ratio<-exp(logratio)

#b and c are both 0 and d is NaN.  What is going on?

#Let's try to optimise
oologscale<-optim(5,logfun,lower=-10,upper=10,method='L-BFGS-B',control=list(fnscale=-1))
#The maximum is at x=0
oologscale$par