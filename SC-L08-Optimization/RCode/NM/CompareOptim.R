rm(list=ls())

#Define Sphere Function
f<-function(x){x[1]^2+x[2]^2}
g<-function(x){c(2*x[1],2*x[2])}



#Use Nelder Mead and optim for trace
x0<-c(2,2)
onm<-optim(x0,f)
#Compare to BFGS
obfgs<-optim(x0,f,method="BFGS") #Optimization by BFGS
olbfgsb<-optim(x0,f,method="L-BFGS-B") #Optimization by L-BFGS-B
obfgsg<-optim(x0,f,g,method="BFGS") #Optimization by BFGS (Derivative Known)
olbfgsbg<-optim(x0,f,g,method="L-BFGS-B") #Optimization by L-BFGS-B (Derivative Known)

booth<-function(x){(x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2}

boothnm<-function(x){
  y<-(x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2
  if ((abs(x[1])>10)||(abs(x[2])>10)){y<-1E30}
  return(y)
}

boothg<-function(x){c(2*(x[1]+2*x[2]-7)+4*(2*x[1]+x[2]-5),
                      4*(x[1]+2*x[2]-7)+2*(2*x[1]+x[2]-5))}

x0<-c(-10,1)
oboothnm<-optim(x0,boothnm) #Optimization by Nelder Mead
oboothqn<-optim(x0,booth,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,10)) #Optimization by L-BFGS-B (Derivative Known)


bukin<-function(x){100*sqrt(abs(x[2]-(x[1]^2/100)))+(abs(x[1]+10)/100)}
bukinnm<-function(x){
  y<-100*sqrt(abs(x[2]-(x[1]^2/100)))+(abs(x[1]+10)/100)
  if ((x[1]>5)||(x[1]<(-15))||(abs(x[2])>3)){y<-1E30}
  return(y)
}
x0<-c(-9.5,1.5)
obukinnm<-optim(x0,bukinnm) #Optimization by Nelder Mead
obukinqn<-optim(x0,bukin,method="L-BFGS-B",lower=c(-15,-3),upper=c(5,3)) #Optimization by L-BFGS-B (Derivative Known)
