library(quantreg)
library(SparseM)
library(KernSmooth)
library(MASS)
library(stats)
library(MultiRNG)
library(hqreg)
library(cqrReg)
library(readxl)
setwd('C:/Users/JR/Desktop')
source('PSSsp.R')
dataallg<-read_excel("YearPredictionMSD.xlsx")
dataallg<-na.omit(dataallg)
dataallg<- as.matrix(dataallg)
dataallg[dataallg==123456789]<-NA
dataallg<-na.omit(dataallg)
dataall<- as.matrix(dataallg) 

p=90
y1<-dataall[,1]
Nall<-length(y1)
yall<-dataall[1:500000,1]
xallg<-dataall[1:500000,2:91]
xall=cbind(1,xallg)
N<-length(yall)

ytest<-dataall[500001:Nall,1]
xtest<-cbind(1,dataall[500001:Nall,2:91])

##############################################
##############################################

tauall=c(0.1,0.3,0.7,0.9, 0.5)
Kall=c(50,100,200,300)
for (tauq in 1:4)
{
  tau=tauall[tauq]
  print(tau)
  
  lambda = cv.hqreg(xallg, yall,method = "quantile", tau = tau)$lambda.min
  t1=Sys.time()
  bglasso1=QR.lasso.admm(xallg,yall,tau)
  ball=c(bglasso1$b,bglasso1$beta)
  t2=Sys.time()
  t_ball=t2-t1
  print(t_ball)
  print(mean(abs(ytest-xtest%*%ball)))
  
  for (tk in 1:4)
  {
    K=Kall[tk]
    print(K)
    n=N/K
    ##################################
    ####### 2- b0-L1  ###################
    ##################################
    t3=Sys.time()
    y0=yall[1:n]
    x0=xall[1:n,]
    x0g=xallg[1:n,]
    b0g=hqreg(x0g, y0, method = "quantile", tau = tau, alpha=1,nlambda=10)$beta
    b0=b0g[,10]
    t4=Sys.time()
    t_b0=t4-t3
    ##################################
    #######end 2- b0-L1  ################
    ##################################
    

    ##################################
    #######3- bAV-QR-L1  #############
    ##################################
    b_AVg=matrix(0,nrow=p+1,ncol=K)
    b_AV=matrix(0,nrow=p+1,ncol=1)
    t5=Sys.time()
    for(ki in 1:K)
    {
      y=yall[((ki-1)*n+1):(ki*n)]
      x=xall[((ki-1)*n+1):(ki*n),] 
      xg=x[1:n,2:(p+1)]
      b0g=hqreg(xg, y, method = "quantile", tau = tau, alpha=1,nlambda=10)$beta
      b_AVg[,ki]=b0g[,10]
    }
    t6=Sys.time()
    t_AV=t6-t5
    b_AV=rowSums(b_AVg)/K
    ##################################
    #######end 3- bAV-L1  ############
    ##################################
    
    ##################################
    #######4- Wang  ##################
    ##################################
    t7=Sys.time()
    bhg=b0
    aa=matrix(0,nrow=K,ncol=(p+1))
    for(ki in 1:K)
    {
      y=yall[((ki-1)*n+1):(ki*n)]
      x=xall[((ki-1)*n+1):(ki*n),]  
      aa[ki,]=t(x)%*%(tau-(y<=(x%*%bhg)))
    }
    dL=t(x0)%*%(tau-(y0<=(x0%*%bhg)))/n-colSums(aa)/N
    itt=0
    cha=1
    th=matrix(0,nrow=n,ncol=1)
    while ((cha>0.0001)*(itt<10))
    {
      c=y0-x0%*%bhg-th/1.2*n
      rt1=(c-tau/1.2)*((c-tau/1.2)>0)-(-c-(1-tau)/1.2)*((-c-(1-tau)/1.2)>0)
      ynew=y0-rt1
      A1=dL-1.2*t(x0)%*%ynew/n+t(x0)%*%th
      A2=0.6*t(x0)%*%x0/n
      bhg1=PSSsp(A1, A2, bhg, lambda) 
      cha=mean(abs(bhg1-bhg))
      itt=itt+1
      bhg=bhg1
      th=th-1.2*(y0-rt1-x0%*%bhg)/n
    }
    b_W=bhg
    t8=Sys.time()
    t_W=t8-t7+t_b0
    ##################################
    #######4- end Wang  ##############
    ##################################
    
    
    n1=n
    s=p
    ##################################
    #######5- chen  ##################
    ##################################
    t9=Sys.time()
    bdsg=b0
    tit=max(1,ceiling(log(n1/N)/log((s+1)^2*log(N)/n1)))
    fk=matrix(0,nrow=1,ncol=K)
    zk=matrix(0,nrow=(p+1),ncol=K)
    sigmak=matrix(0,nrow=(p+1),ncol=K)
    for(qi in 1:tit)
    {
      h=max(((s+1)*log(N)/N)^{1/2},(s+1)^{qi-0.5}*(log(N)/n1)^{qi/2})
      for(ki in 1:K)
      {
        y=yall[((ki-1)*n+1):(ki*n)]
        x=xall[((ki-1)*n+1):(ki*n),]  
        u=(y-x%*%bdsg)/h
        Kh=(105/64-525/64*u^2+735/64*u^4-315/64*u^6)*(abs(u)<=1)
        fk[ki]=sum(Kh)/n/h
      }
      f=max(mean(fk),0.001)
      
      for(ki in 1:K)
      {
        y=yall[((ki-1)*n+1):(ki*n)]
        x=xall[((ki-1)*n+1):(ki*n),]  
        zk[,ki]=t(x)%*%(x%*%bdsg-((y<=x%*%bdsg)-tau)/f)
        sigmak[,ki]=t(x)%*%x%*%bdsg
      }
      znn=colSums(t(zk))/N
      snn=colSums(t(sigmak))/N
      
      A1=snn-t(x0)%*%(x0)%*%bdsg/n1-znn
      A2=0.5*t(x0)%*%(x0)/n1
      bdsg=PSSsp(A1, A2, bdsg, lambda) 
    }
    b_C=bdsg
    t10=Sys.time()
    t_C=t10-t9+t_b0
    
    ##################################
    #######end 5- Chen  ##############
    ##################################
    
    ##################################
    #######6- DSQRL1  ################
    ##################################
    t11=Sys.time()
    bdsg=b0
    h=((s+1)*log(N)/N)^{1/3}
    h1=((s+1)*log(n1)/n1)^{1/3}
    Q=max(2,ceiling(1+3*log(n1/N)/log((s+1)^5*(log(N))^2/n1^2)))
    S1g=matrix(0,nrow=(p+1),ncol=K)
    S2=matrix(0,nrow=(p+1),ncol=(p+1))
    for(qi in 1:Q)
    {
      S2=t(x0)%*%(x0*as.vector(dnorm((y0-x0%*%bdsg)/h1)))/n1/h1
      lambda=0.1*max(sqrt(log(N)/N),(s+1)^(qi*5/6)*(log(n1)/n1)^(qi/3+0.5))
      for(ki in 1:K)
      {
        y=yall[((ki-1)*n+1):(ki*n)]
        x=xall[((ki-1)*n+1):(ki*n),]  
        S1g[,ki]=t(x)%*%(pnorm(-(y-x%*%bdsg)/h)-tau)
      }
      S1=rowSums(S1g)/N
      A1=S1-S2%*%bdsg
      A2=0.5*S2
      bdsg=PSSsp(A1, A2, bdsg, lambda) 
    }
    b_DS=bdsg
    t12=Sys.time()
    t_DS=t12-t11+t_b0
    ##################################
    #######end 6- DSQRL1  ############
    ##################################
    
    B0=mean(abs(ytest-xtest%*%b0))
    BAV=mean(abs(ytest-xtest%*%b_AV))
    BW=mean(abs(ytest-xtest%*%b_W))
    BC=mean(abs(ytest-xtest%*%b_C))
    BDS=mean(abs(ytest-xtest%*%b_DS))
    
    print(round(cbind(B0,BAV,BW,BC,BDS),3))
    print(round(cbind(t_b0,t_AV,t_W,t_C,t_DS),2))
    print(sum(abs(b_DS[2:(p+1)])>10^(-3)))
  }
}