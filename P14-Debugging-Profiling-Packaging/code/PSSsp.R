PSSsp=function(A1, A2, x0, lambda) #fx=t(x0)%*%A2%*%x0+t(x0)%*%A1+lambda*sum(abs(x0)) return x1
{
Lx=t(x0)%*%A2%*%x0+t(x0)%*%A1
fk=Lx+lambda*sum(abs(x0))
dLx=A1+2*A2%*%x0
gk=(dLx+lambda*sign(x0))*(abs(x0)>10^(-4))+(dLx+lambda)*(abs(x0)<10^(-4))*(dLx<(-lambda))+(dLx-lambda)*(abs(x0)<10^(-4))*(dLx>lambda)
d=-min(1,1/sum(abs(gk)))*gk
it=0
while ((it<10)*(max(abs(gk))>10^(-4)))
{
  a=1
  x1=(x0+a*d)*((x0*(x0+a*d))>=0)
  Lx1=t(x1)%*%A2%*%x1+t(x1)%*%A1
  fk1=Lx1+lambda*sum(abs(x1))
  dLx1=A1+2*A2%*%x1
  gk1=(dLx1+lambda*sign(x1))*(abs(x1)>10^(-4))+(dLx1+lambda)*(abs(x1)<10^(-4))*(dLx1<(-lambda))+(dLx1-lambda)*(abs(x1)<10^(-4))*(dLx1>lambda)
  
  itt=0
  while ((itt<10)*(fk1>(fk+0.1*t(gk)%*%(x1-x0))))
  {
    a=a-0.1
    x1=(x0+a*d)*((x0*(x0+a*d))>=0)
    Lx1=t(x1)%*%A2%*%x1+t(x1)%*%A1
    fk1=Lx1+lambda*sum(abs(x1))
    dLx1=A1+2*A2%*%x1
    gk1=(dLx1+lambda*sign(x1))*(abs(x1)>10^(-4))+(dLx1+lambda)*(abs(x1)<10^(-4))*(dLx1<(-lambda))+(dLx1-lambda)*(abs(x1)<10^(-4))*(dLx1>lambda)
    itt=itt+1
  } 
  
  sk=x1-x0
  yk=gk1-gk
  
  if ((sum(abs(sk))==0)|(sum(abs(yk))==0))
  {sigma=1}else 
  {sigma1=t(sk)%*%sk/(t(yk)%*%sk)
  sigma2=t(sk)%*%yk/(t(yk)%*%yk)
  sigma=min(100,sigma1,sigma2)
  }
  d=-sigma*gk1*((sigma*gk1*gk1)>0)
  gk=gk1
  fk=fk1
  x0=x1*(abs(x1)>10^(-4))
  
  it=it+1
}
return(x0)
}


 
 