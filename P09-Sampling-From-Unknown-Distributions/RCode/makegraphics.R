rm(list=ls())

need1x<-c(0.25,0.55)
need1y<-c(0.4,0.8)
need2x<-c(0.7,0.7)
need2y<-c(0.1,0.6)
pdf('../Pics/buffonneedle.pdf')
plot(c(0,0),c(0,1),"l",col='red',lty=2,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',main='Buffon Needle')
for (x in seq(0.2,1,0.2)){
  lines(c(x,x),c(0,1),"l",col='red',lty=2,xlim=c(0,1),ylim=c(0,1))
}

lines(need1x,need1y,lwd=3)
lines(need2x,need2y,lwd=3,col='blue')
dev.off()


cdfd1x<-c(0,0,1,1)
cdfd1y<-c(0,0.3,0.3,1)
set.seed(1)
u<-runif(5)

pdf('../Pics/d1p1.pdf')
plot(cdfd1x,cdfd1y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')
points(rep(-0.02,5),u,pch=19)
dev.off()

pdf('../Pics/d1p2.pdf')
plot(cdfd1x,cdfd1y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')
for (i in 1:5){
  if (u[i]<0.3){
    lines(c(-1,0),c(u[i],u[i]),col='red',lwd=2,lty=1)
  }
  else{
    lines(c(-1,1),c(u[i],u[i]),col='blue',lwd=2,lty=1)
  }
}
dev.off()

pdf('../Pics/d1p3.pdf')
plot(cdfd1x,cdfd1y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')

for (i in 1:5){
  if (u[i]<0.3){
    lines(c(-1,0),c(u[i],u[i]),col='red',lwd=2,lty=1)
    lines(c(0,0),c(0,u[i]),col='red',lwd=2,lty=1)  
  }
  else{
    lines(c(-1,1),c(u[i],u[i]),col='blue',lwd=2,lty=1)
    lines(c(1,1),c(0,u[i]),col='blue',lwd=2,lty=1)
  }
}
dev.off()


cdfd2x<-c(0,0,1,1,2,2)
cdfd2y<-c(0,0.3,0.3,0.55,0.55,1)

pdf('../Pics/d2p1.pdf')
plot(cdfd2x,cdfd2y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')
points(rep(-0.02,5),u,pch=19)
dev.off()

pdf('../Pics/d2p2.pdf')
plot(cdfd2x,cdfd2y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')
for (i in 1:5){
  if (u[i]<0.3){
    lines(c(-1,0),c(u[i],u[i]),col='red',lwd=2,lty=1)
  }
  if (u[i]>0.55){
    lines(c(-1,2),c(u[i],u[i]),col='magenta',lwd=2,lty=1)
  }
  if ((u[i]<0.55)&&(u[i]>0.3)){
    lines(c(-1,1),c(u[i],u[i]),col='blue',lwd=2,lty=1)
  }
}
dev.off()

pdf('../Pics/d2p3.pdf')
plot(cdfd2x,cdfd2y,"l",lwd=3,main="Example 1",xlab='X',ylab='U')

for (i in 1:5){
  if (u[i]<0.3){
    lines(c(-1,0),c(u[i],u[i]),col='red',lwd=2,lty=1)
    lines(c(0,0),c(0,u[i]),col='red',lwd=2,lty=1)  
  }
  if (u[i]>0.55){
    lines(c(-1,2),c(u[i],u[i]),col='magenta',lwd=2,lty=1)
    lines(c(2,2),c(0,u[i]),col='magenta',lwd=2,lty=1)  
  }
  if ((u[i]<0.55)&&(u[i]>0.3)){
    lines(c(-1,1),c(u[i],u[i]),col='blue',lwd=2,lty=1)
    lines(c(1,1),c(0,u[i]),col='blue',lwd=2,lty=1)
  }
}
dev.off()

gridx<-seq(-2,2,0.01)
gridy<-pnorm(gridx)
y<-qnorm(u)

pdf('../Pics/cp1.pdf')

plot(gridx,gridy,"l",main='Continuous Case',xlab='X',ylab='U',lwd=3)
points(rep(-2.01,5),u,pch=19)
dev.off()

pdf('../Pics/cp2.pdf')

plot(gridx,gridy,"l",main='Continuous Case',xlab='X',ylab='U',lwd=3)
points(rep(-2.01,5),u,pch=19)
for (i in 1:5){
  lines(c(-2.01,y[i]),c(u[i],u[i]),lwd=2,lty=1)
}
dev.off()

pdf('../Pics/cp3.pdf')

plot(gridx,gridy,"l",main='Continuous Case',xlab='X',ylab='U',lwd=3)
points(rep(-2.01,5),u,pch=19)
for (i in 1:5){
  lines(c(-2.01,y[i]),c(u[i],u[i]),lwd=2,lty=1)
  lines(c(y[i],y[i]),c(0,u[i]),lwd=2,lty=1)
}
points(y,rep(0,5),pch=19)
dev.off()

gridx<-seq(0,1,0.01)
gridy<-sqrt(1-gridx^2)
pdf('../Pics/appi1.pdf')
plot(gridx,gridy,"l",main='Approximate pi',xlab='',ylab='')
polygon(c(0,gridx),c(0,gridy),col='red')
dev.off()
n<-100
u1<-runif(n)
u2<-runif(n)

pdf('../Pics/appi2.pdf')
plot(gridx,gridy,"l",main='Approximate pi',xlab='',ylab='')
points(u1,u2,pch=16)
dev.off()

pdf('../Pics/appi3.pdf')
plot(gridx,gridy,"l",main='Approximate pi',xlab='',ylab='')
points(u1[(u1^2+u2^2<1)],u2[(u1^2+u2^2<1)],pch=16,col='red')
points(u1[(u1^2+u2^2>1)],u2[(u1^2+u2^2>1)],pch=16,col='blue')
dev.off()

pdf('../Pics/appi4.pdf')
plot(gridx,gridy,"l",main='Approximate pi',xlab='',ylab='')
points(u1[(u1^2+u2^2<1)],u2[(u1^2+u2^2<1)],pch=16,col='red')
dev.off()

x<-seq(0,1,0.01)
y<-dbeta(x,2,2)

pdf('../Pics/bet1.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
dev.off()
n<-100
u1<-runif(n)
u2<-2*runif(n)

pdf('../Pics/bet2.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
points(u1,u2,pch=16)
dev.off()

pdf('../Pics/bet3.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
points(u1[(u2<dbeta(u1,2,2))],u2[(u2<dbeta(u1,2,2))],pch=16,col='red')
points(u1[(u2>dbeta(u1,2,2))],u2[(u2>dbeta(u1,2,2))],pch=16,col='blue')
dev.off()

pdf('../Pics/bet4.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
points(u1[(u2<dbeta(u1,2,2))],u2[(u2<dbeta(u1,2,2))],pch=16,col='red')

dev.off()

pdf('../Pics/bet5.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
lines(c(0.5,0.5),c(0,dbeta(0.5,2,2)),col='magenta',lty=2)
points(u1[(u2<dbeta(u1,2,2))],u2[(u2<dbeta(u1,2,2))],pch=16,col='red')
points(u1[(u2<dbeta(u1,2,2))&(u1>0.5)],u2[(u2<dbeta(u1,2,2))&(u1>0.5)],pch=16,col='magenta')

dev.off()

a<-qbeta(0.6,2,2)
b<-qbeta(0.8,2,2)
pdf('../Pics/bet6.pdf')
plot(x,y,"l",lwd=2,main='Beta Distribution',ylim=c(0,2))
lines(c(a,a),c(0,dbeta(a,2,2)),col='magenta',lty=2)
lines(c(b,b),c(0,dbeta(b,2,2)),col='magenta',lty=2)
points(u1[(u2<dbeta(u1,2,2))],u2[(u2<dbeta(u1,2,2))],pch=16,col='red')
points(u1[(u2<dbeta(u1,2,2))&(u1>a)&(u1<b)],u2[(u2<dbeta(u1,2,2))&(u1>a)&(u1<b)],pch=16,col='magenta')

dev.off()

xgrid<-seq(-5,5,0.01)
nm<-dnorm(xgrid)
lg<-dlogis(xgrid)
M<-dnorm(0)/dlogis(0)
Mlg<-M*lg
pdf('../Pics/nmlg1.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y')
lines(xgrid,lg,lwd=3,col='blue')
dev.off()
pdf('../Pics/nmlg2.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y')
lines(xgrid,Mlg,lwd=3,col='blue')
dev.off()
n<-500
xp<-rlogis(n)
yp<-M*dlogis(xp)*runif(n)

pdf('../Pics/nmlg3.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y')
lines(xgrid,Mlg,lwd=3,col='blue')
points(xp,yp,col='blue',pch=16,cex=0.5)
points(xp[(yp<dnorm(xp))],yp[(yp<dnorm(xp))],cex=0.5,pch=16,col='red')
dev.off()

pdf('../Pics/nmlg4.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y')
lines(xgrid,Mlg,lwd=3,col='blue')
lines(c(xp[1],xp[1]),c(0,M*dlogis(xp[1])),lty=2)
text(xp[1],M*dlogis(xp[1]),'y=Mq(0.96)',pos=4)
text(xp[1],dnorm(xp[1]),'y=p(0.96)',pos=4)
text(xp[1],0.15,paste('x=',round(xp[1],2)),pos=2)
dev.off()

pdf('../Pics/nmlg5.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y',xlim=c(-1,2))
lines(xgrid,Mlg,lwd=3,col='blue')
lines(c(xp[1],xp[1]),c(dnorm(xp[1]),M*dlogis(xp[1])),lty=1,col='blue')
lines(c(xp[1],xp[1]),c(0,dnorm(xp[1])),lty=1,col='red')
text(xp[1],0,'y=0',pos=4)
text(xp[1],M*dlogis(xp[1]),'y=Mq(0.96)',pos=4)
text(xp[1],dnorm(xp[1]),'y=p(0.96)',pos=4)
dev.off()

pdf('../Pics/nmlg6.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y',xlim=c(-1,2))
lines(xgrid,Mlg,lwd=3,col='blue')
lines(c(xp[1],xp[1]),c(dnorm(xp[1]),M*dlogis(xp[1])),lty=1,col='blue')
lines(c(xp[1],xp[1]),c(0,dnorm(xp[1])),lty=1,col='red')
points(xp[1],0.2,pch=18,col='red',cex=3)
text(xp[1],0,'y=0',pos=4)
text(xp[1],M*dlogis(xp[1]),'y=Mq(0.96)',pos=4)
text(xp[1],dnorm(xp[1]),'y=p(0.96)',pos=4)
dev.off()

pdf('../Pics/nmlg7.pdf')
plot(xgrid,nm,"l",lwd=3,col='red',main='Simulate from Normal',xlab = 'x',ylab = 'y',xlim=c(-1,2))
lines(xgrid,Mlg,lwd=3,col='blue')
lines(c(xp[1],xp[1]),c(dnorm(xp[1]),M*dlogis(xp[1])),lty=1,col='blue')
lines(c(xp[1],xp[1]),c(0,dnorm(xp[1])),lty=1,col='red')
points(xp[1],0.28,pch=18,col='blue',cex=3)
text(xp[1],0,'y=0',pos=4)
text(xp[1],M*dlogis(xp[1]),'y=Mq(0.96)',pos=4)
text(xp[1],dnorm(xp[1]),'y=p(0.96)',pos=4)
dev.off()

pdf('../Pics/rat1.pdf')
plot(xgrid,dnorm(xgrid)/dlogis(xgrid),"l",lwd=3,main='Ratio Target/Proposal',xlab = 'x',ylab = 'y')
dev.off()

pdf('../Pics/rat2.pdf')
plot(xgrid,dlogis(xgrid)/dnorm(xgrid),"l",lwd=3,main='Ratio Target/Proposal',xlab = 'x',ylab = 'y')
dev.off()
