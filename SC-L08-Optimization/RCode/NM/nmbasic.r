rm(list=ls())
#Set Function
f<-function(x){x[1]^2+x[2]^2}

#Set dimension
n<-2

#Set Parameters
alpha<-1 #Reflection
gamma<-2 #Expansion
rho<-0.5 #Contraction
sigma<-0.5 #Shrinkage

#Initialise Simplex
simplex<-matrix(c(1,1,2,
                  2,1,2),n,n+1,byrow = TRUE)

#Evaluate function at Simplex
fx<-rep(NA,n+1)
for(j in 1:(n+1)){fx[j]<-f(simplex[,j])}

#Sort by function values
sortf<-sort(fx,index.return=TRUE)
  
#Reorder points in simplex
simplex<-simplex[,sortf$ix]
fx<-sortf$x
  
#Compute Centroid
cent<-apply(simplex[,1:n],1,mean)
  
#Compute Reflected Point
refl<-cent+alpha*(cent-simplex[,(n+1)])
frefl<-f(refl)
  
###Three Cases###
  
###Case 1###
  
if((fx[1]<=frefl)&&(frefl<fx[n])){
  simplex[,n+1]<-refl #Replace Worst Point with reflected point
  fx[n+1]<-frefl #Update fx
}
  
###Case 2###

if(frefl<fx[1]){
  #Compute Expansion
  expa<-cent+gamma*(refl-cent)  
  fexpa<-f(expa)
  if(fexpa<frefl){
    simplex[,n+1]<-expa #Replace Worst Point with expanded point  
    fx[n+1]<-fexpa #Update fx
  }
  else {
    simplex[,n+1]<-refl #Replace Worst Point with reflected point
    fx[n+1]<-frefl #Update fx
  }
}
  
###Case 3###

if(frefl>=fx[n]){
  #Compute Contraction
  cont<-cent+rho*(simplex[,n+1]-cent)  
  fcont<-f(cont)
  if(fcont<fx[n+1]){
    simplex[,n+1]<-cont #Replace Worst Point with contracted point
    fx[n+1]<-fcont #Update fx
  }
  else {
    #Shrink
    for (j in 2:(n+1)){
      simplex[,j]<-simplex[,1]+sigma*(simplex[,j]-simplex[,1])
      fx[j]<-f(simplex[,j])
    }
  }
}
print('Simplex:')
print(simplex)
