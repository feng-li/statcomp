rm(list=ls())

eta<-1.6
tau2<-0.001
sumy<-1734
n<-1000
sigma2<-0.01
b<-1/((1/tau2)+(n/sigma2))  
a<-((sumy/sigma2)+(eta/tau2))*b
print(a)
print(b)
x<-seq(1.55,1.85,0.001)
posterior1<-dnorm(x,a,sqrt(b))
plot(x,posterior1,'l')
eta<-1.6
tau2<-0.25
b<-1/((1/tau2)+(n/sigma2))  
a<-((sumy/sigma2)+(eta/tau2))*b
print(a)
print(b)
posterior2<-dnorm(x,a,sqrt(b))
pdf('../Pics/compare1.pdf')
plot(x,posterior1,"l",lwd=3,xlab=expression(mu),ylab=expression(paste('p(',mu,'|y)')),main='Effect of Prior')
lines(x,posterior2,"l",lwd=3,col='red')
text(1.6,10,labels = expression(tau^2==0.001))
text(1.8,8,labels = expression(tau^2==0.25),col='red')
dev.off()

sumy<-1734
n<-1000
b<-1/((1/tau2)+(n/sigma2))  
a<-((sumy/sigma2)+(eta/tau2))*b
print(a)
print(b)
posterior3<-dnorm(x,a,sqrt(b))
pdf('../Pics/compare2.pdf')
plot(x,posterior3,"l",lwd=3,xlab=expression(mu),ylab=expression(paste('p(',mu,'|y)')),main='Effect of Prior',col='red')
lines(x,posterior1,"l",lwd=3)
text(1.76,60,labels = expression(n==1000),col='red')
text(1.6,10,labels = expression(n==10))
dev.off()
