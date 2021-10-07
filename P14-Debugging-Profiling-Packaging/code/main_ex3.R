library(quantreg)
library(SparseM)
library(KernSmooth)
library(MASS)
library(stats)
library(MultiRNG)
library(hqreg)
library(cqrReg)
setwd('C:/Users/JR/Desktop')
source('PSSsp.R')
####### simple begin#########
N=10^6 ###全样本数
nite=100 ###模拟迭代数
p=100
s=3
b1=c(1,2,3)
b2=matrix(0,nrow=p-s,ncol=1)
beta_true=matrix(c(1,b1,b2),nrow=p+1,ncol=1)
Kall=c(100,200,500,800,1000)
KL=length(Kall)
tauall=c(0.5,0.1,0.9)

t_ball=matrix(0,nrow=nite,ncol=1)
t_DS=matrix(0,nrow=nite,ncol=KL)
t_AV=matrix(0,nrow=nite,ncol=KL)
t_C=matrix(0,nrow=nite,ncol=KL)
t_0=matrix(0,nrow=nite,ncol=KL)
t_W=matrix(0,nrow=nite,ncol=KL)

MSE_ball=matrix(0,nrow=nite,ncol=1)
MSE_b0=matrix(0,nrow=nite,ncol=KL)
MSE_C=matrix(0,nrow=nite,ncol=KL)
MSE_AV=matrix(0,nrow=nite,ncol=KL)
MSE_DS=matrix(0,nrow=nite,ncol=KL)
MSE_W=matrix(0,nrow=nite,ncol=KL)

C_ball=matrix(0,nrow=nite,ncol=1)
C_b0=matrix(0,nrow=nite,ncol=KL)
C_C=matrix(0,nrow=nite,ncol=KL)
C_AV=matrix(0,nrow=nite,ncol=KL)
C_DS=matrix(0,nrow=nite,ncol=KL)
C_W=matrix(0,nrow=nite,ncol=KL)

IC_b0=matrix(0,nrow=nite,ncol=KL)
IC_C=matrix(0,nrow=nite,ncol=KL)
IC_AV=matrix(0,nrow=nite,ncol=KL)
IC_DS=matrix(0,nrow=nite,ncol=KL)
IC_ball=matrix(0,nrow=nite,ncol=1)
IC_W=matrix(0,nrow=nite,ncol=KL)

ss=matrix(0,nrow=p,ncol=p)
for (si in 1:p)
{
  for (sj in 1:p)
  {
    ss[si,sj]=0.5^(abs(si-sj))
  }
}

for (tauq in 1:3)
{
  tau=tauall[tauq]
  print(tau)
    ##################################
    #############sample###############
    ##################################
    for (nii in 1:nite)
    {
      print(nii)
      
      
      xallg=draw.d.variate.uniform(no.row=N,d=p,cov.mat=ss)
      xall=cbind(1,xallg)
      ##################################
      #error=rnorm(N,0,1)
      #yall=xall%*%beta_true+(error-quantile(error,tau))
      ##################################
      error=rt(N,3)
      yall=xall%*%beta_true+0.5*(1+(xallg[,1])^2)*(error-quantile(error,tau))
      ##################################
      
      ##################################
      #######1- ball  ###################
      ##################################
      t1=Sys.time()
      bglasso1=QR.cd(xallg,yall,tau)
      ball=c(bglasso1$b,bglasso1$beta)
      t2=Sys.time()
      t_ball[nii]=t2-t1
      ##################################
      #######end 1- ball  ###############
      ##################################
      MSE_ball[nii]=sqrt(sum((ball-beta_true)^2))
      C_ball[nii]=mean(ball[1:(s+1)]>10^(-3))
      IC_ball[nii]=mean(ball[(s+2):(p+1)]>10^(-3))
      lambda = cv.hqreg(xallg, yall,method = "quantile", tau = tau)$lambda.min
      
      for (Kq in 1:KL)
      {
       K=Kall[Kq] ###分块数
        #print(K)
        n=N/K   
      ##################################
      ####### 2- b0  ###################
      ##################################
      n1=n
      t3=Sys.time()
      y0=yall[1:n]
      x0=xall[1:n,]
      x0g=xallg[1:n,]
      b0g=hqreg(x0g, y0, method = "quantile", tau = tau, alpha=1,nlambda=10)$beta
      b0=b0g[,10]
      t4=Sys.time()
      t_b0=t4-t3
      t_0[nii,Kq]=t_b0
      ##################################
      #######end 2- b0  ################
      ##################################
      
      ##################################
      #######3- bAV-QR  ###################
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
      t_AV[nii,Kq]=t6-t5
      b_AV=rowSums(b_AVg)/K
      ##################################
      #######end 3- bAV  ###############
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
      t_C[nii,Kq]=t10-t9+t_b0
      
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
      Q=max(1,ceiling(1+3*log(n1/N)/log((s+1)^5*(log(N))^2/n1^2)))
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
      t_DS[nii,Kq]=t12-t11+t_b0
      ##################################
      #######end 6- DSQRL1  ############
      ##################################
      
      MSE_b0[nii,Kq]=sqrt(sum((b0-beta_true)^2))
      MSE_C[nii,Kq]=sqrt(sum((b_C-beta_true)^2))
      MSE_DS[nii,Kq]=sqrt(sum((b_DS-beta_true)^2))
      MSE_AV[nii,Kq]=sqrt(sum((b_AV-beta_true)^2))
      MSE_W[nii,Kq]=sqrt(sum((b_W-beta_true)^2))
      
      C_b0[nii,Kq]=mean(b0[1:(s+1)]>10^(-3))
      C_C[nii,Kq]=mean(b_C[1:(s+1)]>10^(-3))
      C_DS[nii,Kq]=mean(b_DS[1:(s+1)]>10^(-3))
      C_AV[nii,Kq]=mean(b_AV[1:(s+1)]>10^(-3))
      C_W[nii,Kq]=mean(b_W[1:(s+1)]>10^(-3))
      
      IC_b0[nii,Kq]=mean(b0[(s+2):(p+1)]>10^(-3))
      IC_C[nii,Kq]=mean(b_C[(s+2):(p+1)]>10^(-3))
      IC_DS[nii,Kq]=mean(b_DS[(s+2):(p+1)]>10^(-3))
      IC_AV[nii,Kq]=mean(b_AV[(s+2):(p+1)]>10^(-3))
      IC_W[nii,Kq]=mean(b_W[(s+2):(p+1)]>10^(-3))
    }
    }

 
  for (Kq in 1:KL)
  {
    K=Kall[Kq] ###分块数
    print(K)
  MSE_mean=cbind(mean(MSE_ball),mean(MSE_b0[,Kq]),mean(MSE_AV[,Kq]),mean(MSE_W[,Kq]),mean(MSE_C[,Kq]),mean(MSE_DS[,Kq]))
  MSE_sd=cbind(sd(MSE_ball),sd(MSE_b0[,Kq]),sd(MSE_AV[,Kq]),sd(MSE_W[,Kq]),sd(MSE_C[,Kq]),sd(MSE_DS[,Kq]))
  
  C_mean=cbind(mean(C_ball),mean(C_b0[,Kq]),mean(C_AV[,Kq]),mean(C_W[,Kq]),mean(C_C[,Kq]),mean(C_DS[,Kq]))
  IC_mean=cbind(mean(IC_ball),mean(IC_b0[,Kq]),mean(IC_AV[,Kq]),mean(IC_W[,Kq]),mean(IC_C[,Kq]),mean(IC_DS[,Kq]))
  IC_sd=cbind(sd(IC_ball),sd(IC_b0[,Kq]),sd(IC_AV[,Kq]),sd(IC_W[,Kq]),sd(IC_C[,Kq]),sd(IC_DS[,Kq]))
  
  t=cbind(mean(t_ball),mean(t_0[,Kq]),mean(t_AV[,Kq]),mean(t_W[,Kq]),mean(t_C[,Kq]),mean(t_DS[,Kq]))
  print(round(MSE_mean,3))
  print(round(MSE_sd,3))
  print(round(IC_mean,3))
  print(round(IC_sd,3))
  print(round(C_mean,3))
  print(round(t,3))
  }
}
