rm(list=ls())
g1<-301
g2<-301
gridx1<-seq(-0.9,2.1,length.out = g1)
gridx2<-seq(-0.9,2.1,length.out = g2)
count<-0
ymat<-matrix(0,g1,g2)
for (i in 1:g1){
  for (j in 1:g2){
    ymat[i,j]<-gridx1[i]^2+gridx2[j]^2+3*exp(-(gridx1[i]+1)^2)*exp(-(gridx2[j]+1)^2)
  }
}



#Set Function
f<-function(x){x[1]^2+x[2]^2+3*exp(-(x[1]+1)^2)*exp(-(x[2]+1)^2)}

#Set dimension
n<-2

#Set Parameters
alpha<-1 #Reflection
gamma<-2 #Expansion
rho<-0.5 #Contraction
sigma<-0.5 #Shrinkage

#Initialise Simplex
simplex<-matrix(c(-0.7,-0.7,-0.4,
                  -0.4,-0.7,-0.7),n,n+1,byrow = TRUE)

png('out/out00.png')
image(gridx1,gridx2,ymat,main='Nelder-Mead',xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(-0.9,1.5),ylim=c(-0.9,1.5),zlim = c(0,2))
polygon(t(simplex),lwd=4)

dev.off()

#Volume of Simplex and tolerance(for stopping condition)
volume<-0.5*abs(det(simplex[,1:n]-simplex[,(n+1)]))
tol<-1E-10

#Maximum Iterations (and initialise iteration count)
maxiter<-100
iter<-0

#Evaluate function at Simplex
fx<-rep(NA,n+1)
for(i in 1:(n+1)){fx[i]<-f(simplex[,i])}

while ((volume>tol)&(iter<maxiter)){
  iter<-iter+1
  print(paste('Iteration:',iter))
  #Sort by function values
  sortf<-sort(fx,decreasing = FALSE,index.return=TRUE)
  
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
  volume<-0.5*abs(det(simplex[,1:n]-simplex[,(n+1)]))
  print('Simplex:')
  print(simplex)
  print(paste('Volume:',volume))
  fn<-ifelse((iter<10),paste('out/out0',iter,'.png',sep = ''),paste('out/out',iter,'.png',sep = ''))
  png(fn)
  image(gridx1,gridx2,ymat,main='Nelder-Mead',xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(-0.9,1.5),ylim=c(-0.9,1.5),zlim = c(0,2))
  polygon(t(simplex),lwd=4)
  dev.off()
}
