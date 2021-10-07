rm(list=ls())

#First define function


g<-function(x){(3*(x^5))-(4*(x^4))+(6*(x^3))+4*x-4}
#g<-function(x){abs(x)^0.5}
#g<-function(x){x*exp(-x^2)-0.4*(1/(exp(x)+1))+0.2}
#g<-function(x){2*(x^3)-3*(x^2)+5}
#g<-function(x){(x^3)-2*(x^2)-11*x+12}

#Then define first derivative


dg<-function(x){(15*(x^4))-(16*(x^3))+(18*(x^2))+4}
#dg<-function(x){ifelse((x>0),0.5*(x^(-0.5)),-0.5*((-x)^(-0.5)))}
#dg<-function(x){(1-2*(x^2))*exp(-x^2)+0.4*(exp(x)/(exp(x)+1)^2)}
#dg<-function(x){6*x*(x-1)}
#dg<-function(x){3*(x^2)-4*x-11}  



#Code 1 Newton-Raphson using "for" loop

nralgo1<-function(x0,func,deri){
  epsilon<-1E-10 #Set Tolerance Level
  MaxIter<-500 #Maximum Number of Iterations
  x<-rep(NA,MaxIter+1) #A vector to store sequence of x
  x[1]<-x0 # Initialise
  for(n in 1:MaxIter){
    x[n+1]<-x[n]-(func(x[n])/deri(x[n])) #Update Step
    if(abs(func(x[n+1]))<=epsilon){break} #Stopping Rule
  }
  if(n==MaxIter){warning('Maximum number of iterations reached without convergence')}
  
  return(x[n+1]) #Return value of x
}

#CODE 2: Newton-Raphson using "while"

nralgo2<-function(x0,func,deri){
  epsilon<-1E-10 #Tolerance Level
  x<-x0 #Initialise
  while(abs(func(x))>epsilon){
    x<-x-(func(x)/deri(x))  #Update step
    print(x)
  }
  return(x)
}

#With all functions defined we can now do optimisation:

x0<-1.4
#Use Code 1
root1=nralgo1(x0,g,dg)
print(paste('Root using Code 1:',round(root1,3)))
print(paste('Function Evaluated at this Root',round(g(root1),3)))

#Using Code 2
root2<-nralgo2(x0,g,dg)  
print(paste('Root using Code 2:',round(root2,3)))
print(paste('Function Evaluated at this Root',round(g(root2),3)))