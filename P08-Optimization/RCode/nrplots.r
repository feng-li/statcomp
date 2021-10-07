rm(list=ls()) # Clear Workspace

gx<-function(x){
  y=(3*(x^5))-(4*(x^4))+(6*(x^3))+(4*x)-4
  return(y)  
}

dgx<-function(x){
  y=(15*(x^4))-(16*(x^3))+(18*(x^2))+4
  return(y)  
}

pdf('rfout.pdf')
gridx<-seq(-2,2,0.01)
gridy<-gx(gridx)
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4))
lines(c(-2,2),c(0,0))
dev.off()

x0<-1.4
gridx<-seq(0.2,1.5,0.01)
gridy<-gx(gridx)

pdf('rf1.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
dev.off()

pdf('rf2.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)
dev.off()

x1<-x0-(gx(x0)/dgx(x0))

pdf('rf3.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)
lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
dev.off()


pdf('rf4.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

dev.off()


x2<-x1-(gx(x1)/dgx(x1))

pdf('rf5.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

lines(c((x1),x2),c(gx(x1),0),lwd=2,col='blue')
dev.off()

pdf('rf6.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=20)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
lines(c((x1),x2),c(gx(x1),0),lwd=2,col='blue')

text(x2,0,bquote(paste('x'['2']*' = ', .(round(x2,3)))),pos = 1,col='blue')

points(x2,0,pch=18,cex=2,col='blue')
lines(c(x2,x2),c(0,gx(x2)),lwd=2,lty=2,col='blue')
lines(c(0,x2),c(gx(x2),gx(x2)),lwd=2,lty=2,col='blue')

dev.off()
x3<-x2-(gx(x2)/dgx(x2))

pdf('rf7.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding', xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

lines(c((x1),x2),c(gx(x1),0),lwd=2,col='blue')

text(x2,0,bquote(paste('x'['2']*' = ', .(round(x2,3)))),pos = 1,col='blue')

points(x2,0,pch=18,cex=2,col='blue')
lines(c(x2,x2),c(0,gx(x2)),lwd=2,lty=2,col='blue')
lines(c(0,x2),c(gx(x2),gx(x2)),lwd=2,lty=2,col='blue')

lines(c((x2),x3),c(gx(x2),0),lwd=2,col='magenta')

#text(x3,0,paste('x3=',round(x3,3)),pos = 3,col='magenta')
text(x3,0,bquote(paste('x'['3']*' = ', .(round(x3,3)))),pos = 3,col='magenta')

points(x3,0,pch=18,cex=2,col='magenta')
lines(c(x3,x3),c(0,gx(x3)),lwd=2,lty=2,col='magenta')
lines(c(0,x3),c(gx(x3),gx(x3)),lwd=2,lty=2,col='magenta')

dev.off()

pdf('findtangent1.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding',xlim=c(0.2,1.6), xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

#text(1.4,9,expression(g(x[0])),pos=4,col='blue',cex=2)
#text(1.2,-1,expression(x[1]-x[0]),pos=1,col='blue',cex=2)
text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
dev.off()


pdf('findtangent2.pdf')
plot(gridx,gridy,"l",lwd=3,main='Root Finding',xlim=c(0.2,1.6), xlab = 'x', ylab=expression(y=3*x^5-4*x^4+6*x^3+4*x-4)) #Plot
lines(c(-2,2),c(0,0)) #x axis
text(x0,0,bquote(paste('x'['0']*' = ', .(round(x0,3)))),pos = 1) 
points(x0,0,pch=18,cex=2) 
lines(c(x0,x0),c(0,gx(x0)),lwd=2,lty=2)
lines(c(0,x0),c(gx(x0),gx(x0)),lwd=2,lty=2)

lines(c((x0),x1),c(gx(x0),0),lwd=2,col='red')
lines(c(0,x1),c(gx(x1),gx(x1)),lwd=2,lty=2,col='red')

lines(c(x0,x0),c(0,gx(x0)),lwd=2.5,lty=1,col='blue')

text(1.4,9,expression(g(x[0])),pos=4,col='blue',cex=2)
lines(c(x1,x0),c(0,0),lwd=2.5,lty=1,col='magenta')

text(1.2,-1,expression(x[0]-x[1]),pos=1,col='magenta',cex=2)
text(x1,0,bquote(paste('x'['1']*' = ', .(round(x1,3)))),pos = 1,col='red')
points(x1,0,pch=18,cex=2,col='red')
lines(c(x1,x1),c(0,gx(x1)),lwd=2,lty=2,col='red')
dev.off()
