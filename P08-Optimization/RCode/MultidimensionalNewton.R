rm(list=ls())
#There are two ways to code up the functions.  Either with two
#scalar inputs or one vector input

#func<-function(x1,x2){x1^2-(x1*x2)+(x2^2)+exp(x2)}
#grad<-function(x1,x2){c(2*x1-x2,-x1+2*x2+exp(x2))}
#hess<-function(x1,x2){matrix(c(2,-1,-1,2+exp(x2)),2,2)}



func<-function(x){x[1]^2-(x[1]*x[2])+(x[2]^2)+exp(x[2])} #Function
grad<-function(x){c(2*x[1]-x[2],-x[1]+2*x[2]+exp(x[2]))} #Gradient
hess<-function(x){matrix(c(2,-1,-1,2+exp(x[2])),2,2)} #Hessian

#Function for Newton's Method

nralgo<-function(x0,f,g,h){
  epsilon<-1E-10 #Tolerance
  MaxIter<-500 #Maximum Iterations
  n<-1 #Initialise n
  xnew<-x0-solve(h(x0),g(x0)) #Initialise xnew
  change<-abs(f(xnew)-f(x0)) #Initialise change in function value
  while((n<=MaxIter)&&(change>epsilon)){ #Note while loop breaks if either MaxIter reached of change in function value is less than tolerance
    n<-n+1 #Update n
    xnew->xold #New becomes old. Pay attention to the direction of the arrow
    xnew<-xold-solve(h(xold),g(xold)) #xnew is updated
    change<-abs(f(xnew)-f(xold)) #Change in function value updated
    }
    return(xnew)
}

x0<-c(2,1) #Initial value
x_Newton<-nralgo(x0,func,grad,hess) #Solution using Newton
o_BFGS<-optim(x0,func,grad,method="BFGS")
x_BFGS<-o_BFGS$par
o_LBFGSB<-optim(x0,func,grad,method="L-BFGS-B")
x_LBFGSB<-o_LBFGSB$par

#Print out all answers.  All give the same root
print(paste('Newton Method:',x_Newton))
print(paste(' BFGS Method:',x_BFGS))
print(paste('L-BFGS-B Method:',x_LBFGSB))

