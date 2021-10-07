## Direct Method##

##Generate from Weibull (a=1,b=2)
n<-1000
invWeib<-function(p,a,b){b*((-log(-1*(p-1)))^(1/a))}
u<-runif(n)
x1Weib<-invWeib(u,0.5,2)
#Use rweibull to check answer
checkWeib<-rweibull(n,0.5,2)
summary(x1Weib)
summary(checkWeib)

##Generate from Standard Logistic 

invLogi<-function(p){-log((1/p)-1)}
u<-runif(n)
x1Logi<-invLogi(u)
#Use rlogis to check answer

checkLogi<-rlogis(n)
summary(x1Logi)
summary(checkLogi)


##Indirect Method ##

##Q1: LogNormal

#Target Density
p<-function(x){((2*pi)^(-0.5))*(1/x)*exp(-(log(x)^2)/2)}

#Proposal Density (Use Weibull)
q<-function(x,a,b){(a/b)*((x/b)^(a-1))*(exp(-(x/b)^a))}

#Find sup(p(x)/q(x))
rat<-function(x)(-p(x)/q(x,0.5,2)) #NEGATIVE Ratio
oo<-optim(0.1,rat,method = 'L-BFGS-B')#Optim
M<--oo$value #Optimal M


#Use our Weibull sample from previous answer
n<-length(x1Weib)

#Generate y
y<-runif(n)


#Only keep accepted values
accept<-(y<(p(x1Weib)/(M*q(x1Weib,0.5,2))))
xLogN<-x1Weib[accept]

#Use lnorm to check answer
check<-rlnorm(length(xLogN))


summary(xLogN)
summary(check)

##Question 2: Indirect Normal

aralg<-function(n,p,q,rq,...){  #Note use of ellipsis(...).  This allows us to pass arguments for dq and rq 
  #Find M
  rat<-function(x)(-p(x)/q(x,...))#Ratio of negative
  oo<-optim(0.1,rat,method = 'L-BFGS-B') #Use Optim
  M<--oo$value #Optimal M
  #Simulate from proposal
  x<-rq(n,...)
  #Simulate Uniform
  y<-runif(n)
  #Check Accept Condition
  accept<-(y<(p(x)/(q(x,...)*M)))
  #Keep x that satisfy accept condition
  xkeep<-x[accept]
  return(xkeep)
}

n<-10000
#Logistic as proposal
x1<-aralg(n,dnorm,dlogis,rlogis)
length(x1)/n
#Cauchy as proposal
x2<-aralg(n,dnorm,dcauchy,rcauchy)
length(x2)/n

#T with 5 df as proposal
x3<-aralg(n,dnorm,dt,rt,df=5)
length(x3)/n

#T with 5 df as proposal
x4<-aralg(n,dnorm,dt,rt,df=20)
length(x4)/n

#N(1,2) df as proposal
x5<-aralg(n,dnorm,dnorm,rnorm,mean=1,sd=2)
length(x5)/n
