library(quantreg)
library(SparseM)
library(KernSmooth)
library(MASS)
library(stats)
library(MultiRNG)
library(cqrReg)
####### simple begin#########
n=500 ###全样本数
nite=100###模拟迭代数
p=100
beta_trueg=matrix(1,nrow=p,ncol=1)
beta_true=c(1,beta_trueg)

t_ball=matrix(0,nrow=nite,ncol=1)
t_DS=matrix(0,nrow=nite,ncol=1)
t_AV=matrix(0,nrow=nite,ncol=1)
t_LE=matrix(0,nrow=nite,ncol=1)
t_0=matrix(0,nrow=nite,ncol=1)

MSE_ball=matrix(0,nrow=nite,ncol=1)
MSE_b0=matrix(0,nrow=nite,ncol=1)
MSE_LE=matrix(0,nrow=nite,ncol=1)
MSE_AV=matrix(0,nrow=nite,ncol=1)
MSE_DS=matrix(0,nrow=nite,ncol=1)



ss=matrix(0,nrow=p,ncol=p)
for (si in 1:p)
{
  for (sj in 1:p)
  {
    ss[si,sj]=0.5^(abs(si-sj))
  }
}

Kall=c(20,40,60,80,100,200,300,400,500,600,700,800)
tauall=c(0.5,0.3,0.7)
for (tauq in 1:3)
{
  tau=tauall[tauq]
  print(tau)
  for (Kq in 1:12)
  {
    K=Kall[Kq] ###分块数
    print(K)
    N=n*K  
    ##################################
    #############sample###############
    ##################################
    for (nii in 1:nite)
    {
      #print(nii)
      xallg=draw.d.variate.normal(no.row=N,d=p,mean.vec=as.vector(matrix(0,nrow=p,ncol=1)),cov.mat=ss)
      xall=cbind(1,xallg)
      error=rchisq(N,2)
      yall=xall%*%beta_true+(error-quantile(error,tau))
      ##################################
      
      ##################################
      #######1- ball  ###################
      ##################################
      t1=Sys.time()
      ball=rq(yall~0+xall, tau = tau)$coef
      t2=Sys.time()
      t_ball[nii]=t2-t1
      ##################################
      #######end 1- ball  ###############
      ##################################
      
      ##################################
      #######2- b0  ###################
      ##################################
      t3=Sys.time()
      y0=yall[1:n]
      x0=xall[1:n,]  
      b0=rq(y0~0+x0, tau = tau)$coef
      t4=Sys.time()
      t_b0=t4-t3
      t_0[nii]=t4-t3
      ##################################
      #######end 2-b0  ###############
      ##################################
      
      ############3-AVQR################
      ##################################
      b_AVg=matrix(0,nrow=p+1,ncol=K)
      b_AV=matrix(0,nrow=p+1,ncol=1)
      t4=Sys.time()
      for(ki in 1:K)
      {
        y=yall[((ki-1)*n+1):(ki*n)]
        x=xall[((ki-1)*n+1):(ki*n),]  
        b_AVg[,ki]=rq(y~0+x, tau = tau)$coef
      }
      t5=Sys.time()
      t_AV[nii]=t5-t4
      b_AV=rowSums(b_AVg)/K
      ############END 3-AVQR##############
      ##################################
      
      
      ##################################
      #######4- LEQR  ##################
      ##################################
      t6=Sys.time()
      M=ceiling(2+log(log(sqrt((p+1)/N))/log((p+1)/n))/log(2)) ###迭代数
      blqrg=b0
      for(mi in 1:M)
      {
        hg1=max(sqrt((p+1)/N),((p+1)/n)^(2^(mi-2)))
        Ul=matrix(0,nrow=(p+1),ncol=1)
        Vl=matrix(0,nrow=(p+1),ncol=(p+1))
        for(ki in 1:K)
        {
          y=yall[((ki-1)*n+1):(ki*n)]
          x=xall[((ki-1)*n+1):(ki*n),]  
          sg1=(y-x%*%blqrg)/hg1
          S1=(sg1>=1)+(0.5+15/16*(sg1-2/3*sg1^3+0.2*sg1^5))*(abs(sg1)<1)
          Sp1=15/16*(1-2*sg1^2+sg1^4)*(abs(sg1)<1)
          g11=S1+tau-1+y/hg1*Sp1
          txg11l=matrix(0,nrow=n,ncol=(p+1))
          txg21l=matrix(0,nrow=n,ncol=(p+1))
          for (pi in 1:(p+1)) 
          {
            txg11l[,pi]=x[,pi]*Sp1
            txg21l[,pi]=x[,pi]*g11
          }
          Ul=Ul+colSums(txg21l)
          Vl=Vl+t(x)%*%txg11l/hg1
        }
        blqrg=ginv(Vl)%*%Ul
      }
      b_LE=blqrg
      t7=Sys.time()
      t_LE[nii]=t7-t6+t_b0
      ##################################
      #######END- LEQR  ################
      ##################################  
      
      ##################################
      #######5- DSQR  ##################
      ##################################
      t8=Sys.time()
      h=((p+1)/N)^{0.4}
      h1=((p+1)*log(n)/n)^{1/3}
      S1g=matrix(0,nrow=(p+1),ncol=K)
      
        bt0=b0
        for(ki in 1:K)
        {
          y=yall[((ki-1)*n+1):(ki*n)]
          x=xall[((ki-1)*n+1):(ki*n),]  
          S1g[,ki]=t(x)%*%(pnorm(-(y-x%*%bt0)/h)-tau)
        }
        S1=rowSums(S1g)/N
        bt1=bt0-S1
        gt=t(x0)%*%(pnorm(-(y0-x0%*%bt1)/h1)-tau)/n-t(x0)%*%(pnorm(-(y0-x0%*%bt0)/h1)-tau)/n
        it=0
        while ((it<10)*(max(abs(gt))>10^(-3)))
        {
        bt=bt1-bt0
        gt=t(x0)%*%(pnorm(-(y0-x0%*%bt1)/h1)-tau)/n-t(x0)%*%(pnorm(-(y0-x0%*%bt0)/h1)-tau)/n
        if ((sum(abs(bt))==0)|(sum(abs(gt))==0))
        {a=1}else 
        {
        a1=t(bt)%*%bt/(t(bt)%*%gt)
        a2=t(bt)%*%gt/(t(gt)%*%gt)
        a=min(a1,a2,100)}
        bt0=bt1
        
        for(ki in 1:K)
        {
          y=yall[((ki-1)*n+1):(ki*n)]
          x=xall[((ki-1)*n+1):(ki*n),]  
          S1g[,ki]=t(x)%*%(pnorm(-(y-x%*%bt0)/h)-tau)
        }
        S1=rowSums(S1g)/N
        bt1=bt0-a*S1
        it=it+1
      }
      b_DS=bt1
      t9=Sys.time()
      t_DS[nii]=t9-t8+t_b0
      
      ##################################
      #######5- DSQR  ##################
      ##################################
      
      MSE_ball[nii]=sqrt(sum((ball-beta_true)^2))
      MSE_b0[nii]=sqrt(sum((b0-beta_true)^2))
      MSE_LE[nii]=sqrt(sum((b_LE-beta_true)^2))
      MSE_DS[nii]=sqrt(sum((b_DS-beta_true)^2))
      MSE_AV[nii]=sqrt(sum((b_AV-beta_true)^2))
      
      
    }
    
    t=cbind(mean(t_ball),mean(t_AV),mean(t_LE),mean(t_DS))
    MSE_mean=cbind(mean(MSE_ball),mean(MSE_AV),mean(MSE_LE),mean(MSE_DS))
    
    print(round(MSE_mean,3))
    print(round(t,3))
  }
}