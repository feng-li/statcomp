##Simulating from an exponential
set.seed(1)
n<-1000 #Sample Size
u<-runif(n) #Uniform 
invexp<-function(u,lambda){-log(1-u)/lambda} #Inverse cdf
y1<-invexp(u,1) #Using my own inverse cdf
y2<-qexp(u,1) #Using R inverse cdf
y3<-rexp(n,1) #Using R random number generator

#Histograms
par(mfrow=c(1,3))
hist(y1,breaks=0:10,xlim = c(0,10)) #breaks and xlim allows all three  
hist(y2,breaks=0:10,xlim = c(0,10)) #histograms to use the same
hist(y3,breaks=0:10,xlim = c(0,10)) #bins and scale on x axis

#Summary of all three methods
summary(y1)
summary(y2)
summary(y3)

##Estimating pi
n<-1000
u1<-runif(n) #Uniform Variable on x coordinate
u2<-runif(n) #Uniform Variable on y coordinate
inside<-(u1^2+u2^2<1) #Checks whether point lies in quarter circle
pi4<-sum(inside)/n 

print(paste('Estimate of pi/4',pi4))
print(paste('True value of pi/4',pi/4))

##Accept Reject algorithm (Beta distribution)
n<-1000
M<-1.5 #Scaling factor
x<-runif(n) #Generate x coordinate
y<-runif(n) #Generate y coordinate
accept<-(y<dbeta(x,2,2)/M) #Does point lie inside density?
xkeep<-x[accept] #Keep these x's

#Plot histogram
par(mfrow=c(1,1))
hist(xkeep)

#Check number of iterates we keep
length(xkeep)

##Accept Reject algorithm (Logistic distribution)
n<-1000
M<-2
x<-rlogis(n)#Generate x coordinate from proposal
y<-runif(n) #Generate y coordinate from U(0,1)
accept<-(y<dnorm(x)/(dlogis(x)*M)) #Do we accept (TRUE/FALSE)
xkeep<-x[accept] #Keep only TRUE
#Plot histogram
par(mfrow=c(1,1))
hist(xkeep)
#Check number of iterates we keep
length(xkeep)

