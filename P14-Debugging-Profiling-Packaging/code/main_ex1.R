library(quantreg)
library(SparseM)
library(KernSmooth)
library(MASS)
library(stats)
library(MultiRNG)
library(cqrReg)
library(conquer)
library(hqreg)
####### simple begin#########
N=5*10^6 ###全样本数
nite=100###模拟迭代数
p=10
tau=0.5


beta_trueg=matrix(1,nrow=p,ncol=1)
beta_true=c(1,beta_trueg)

t_PK=matrix(0,nrow=nite,ncol=1)
t_H=matrix(0,nrow=nite,ncol=1)
t_P=matrix(0,nrow=nite,ncol=1)
t_DS=matrix(0,nrow=nite,ncol=1)


MSE_PK=matrix(0,nrow=nite,ncol=1)
MSE_H=matrix(0,nrow=nite,ncol=1)
MSE_P=matrix(0,nrow=nite,ncol=1)
MSE_DS=matrix(0,nrow=nite,ncol=1)


ss=matrix(0,nrow=p,ncol=p)
for (si in 1:p)
{
  for (sj in 1:p)
  {
    ss[si,sj]=0.5^(abs(si-sj))
  }
}


    ##################################
    #############sample###############
    ##################################
    for (nii in 1:nite)
    {
      print(nii)
      xallg=draw.d.variate.normal(no.row=N,d=p,mean.vec=as.vector(matrix(0,nrow=p,ncol=1)),cov.mat=ss)
      xall=cbind(1,xallg)
      error=rnorm(N,0,1)
      #error=rt(N,5)
      yall=xall%*%beta_true+error
      ##################################
      
      ##################################
      #######1- Portnoy and Koenker(1997)  ###################
      ##################################
      t1=Sys.time()
      bpk=rq(yall~0+xall, tau = tau, method="fn")$coef
      t2=Sys.time()
      t_PK[nii]=t2-t1
      ##################################
      #######end 1  ###############
      ##################################
      
      
      ##################################
      #######2- He, Pan, Tan, and Zhou (2020) 
      ##################################
      t3=Sys.time()
      bh=conquer(xallg, yall, tau = tau, kernel = "Gaussian")$coef
      t4=Sys.time()
      t_H[nii]=t4-t3
      ##################################
      #######end 2  ###############
      ##################################
      
      ##################################
      #######3- Pietrosanu et al. (2020)  
      ##################################
      t5=Sys.time()
      bg=QR.admm(xallg,yall,tau)
      bp=c(bg$b,bg$beta)
      t6=Sys.time()
      t_P[nii]=t6-t5
      ##################################
      #######end 3  ###############
      ##################################
      
      ##################################
      #######4- DSQR  ##################
      ##################################
      t7=Sys.time()
      h=((p+1)/N)^{0.4}
      b0g=lm(yall ~ xallg)$coef
      gg=IQR(yall)/10
      b0=ginv(t(xall)%*%(xall*(as.vector((abs(yall-xall%*%b0g)<gg)))))%*%(t(xall)%*%(yall*(abs(yall-xall%*%b0g)<gg))-gg*t(xall)%*%(sign(yall-xall%*%b0g)*(abs(yall-xall%*%b0g)>gg)))
      abs=1
      it=0
      while ((it<50)*(abs>10^(-6)))
      {
      S2=t(xall)%*%(xall*as.vector(dnorm((yall-xall%*%b0)/h)))/N/h  
      S1=t(xall)%*%(pnorm(-(yall-xall%*%b0)/h)-tau) /N
      b1=b0-ginv(S2)%*%S1
      abs=sqrt(sum((b1-b0)^2))
      b0=b1
      it=it+1
      }
      b_DS=b1
      t8=Sys.time()
      t_DS[nii]=t8-t7
      
      ##################################
      #######5- DSQR  ##################
      ##################################
   
    
      MSE_PK[nii]=sqrt(sum((bpk-beta_true)^2))
      MSE_H[nii]=sqrt(sum((bh-beta_true)^2))
      MSE_P[nii]=sqrt(sum((bp-beta_true)^2))
      MSE_DS[nii]=sqrt(sum((b_DS-beta_true)^2))
    }
    
    t=cbind(mean(t_PK),mean(t_P),mean(t_H),mean(t_DS))
    MSE_mean=cbind(mean(MSE_PK),mean(MSE_P),mean(MSE_H),mean(MSE_DS))*100
    MSE_sd=cbind(sd(MSE_PK),sd(MSE_P),sd(MSE_H),sd(MSE_DS))*100

    print(round(MSE_mean,3))
    print(round(MSE_sd,3))
    print(round(t,3))