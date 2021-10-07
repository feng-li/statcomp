library(quantreg)
library(SparseM)
library(KernSmooth)
library(MASS)
library(stats)
library(MultiRNG)
library(conquer)
library(cqrReg)
library(readxl)
setwd('C:/Users/JR/Desktop/2')
dataallg1<-read_excel("1.xlsx")
dataallg1<-na.omit(dataallg1)
dataallg1<- as.matrix(dataallg1)
dataallg1[dataallg1==123456789]<-NA
dataallg1<-na.omit(dataallg1)
dataallg1<- as.matrix(dataallg1) 
n1=length(dataallg1[,1])

dataallg2<-read_excel("2.xlsx")
dataallg2<-na.omit(dataallg2)
dataallg2<- as.matrix(dataallg2)
dataallg2[dataallg2==123456789]<-NA
dataallg2<-na.omit(dataallg2)
dataallg2<- as.matrix(dataallg2) 
n2=length(dataallg2[,1])

dataallg3<-read_excel("3.xlsx")
dataallg3<-na.omit(dataallg3)
dataallg3<- as.matrix(dataallg3)
dataallg3[dataallg3==123456789]<-NA
dataallg3<-na.omit(dataallg3)
dataallg3<- as.matrix(dataallg3) 
n3=length(dataallg3[,1])

dataallg4<-read_excel("4.xlsx")
dataallg4<-na.omit(dataallg4)
dataallg4<- as.matrix(dataallg4)
dataallg4[dataallg4==123456789]<-NA
dataallg4<-na.omit(dataallg4)
dataallg4<- as.matrix(dataallg4) 
n4=length(dataallg4[,1])

dataallg5<-read_excel("5.xlsx")
dataallg5<-na.omit(dataallg5)
dataallg5<- as.matrix(dataallg5)
dataallg5[dataallg5==123456789]<-NA
dataallg5<-na.omit(dataallg5)
dataallg5<- as.matrix(dataallg5) 
n5=length(dataallg5[,1])

dataallg6<-read_excel("6.xlsx")
dataallg6<-na.omit(dataallg6)
dataallg6<- as.matrix(dataallg6)
dataallg6[dataallg6==123456789]<-NA
dataallg6<-na.omit(dataallg6)
dataallg6<- as.matrix(dataallg6) 
n6=length(dataallg6[,1])

dataallg7<-read_excel("7.xlsx")
dataallg7<-na.omit(dataallg7)
dataallg7<- as.matrix(dataallg7)
dataallg7[dataallg7==123456789]<-NA
dataallg7<-na.omit(dataallg7)
dataallg7<- as.matrix(dataallg7) 
n7=length(dataallg7[,1])

dataallg8<-read_excel("8.xlsx")
dataallg8<-na.omit(dataallg8)
dataallg8<- as.matrix(dataallg8)
dataallg8[dataallg8==123456789]<-NA
dataallg8<-na.omit(dataallg8)
dataallg8<- as.matrix(dataallg8) 
n8=length(dataallg8[,1])

dataallg9<-read_excel("9.xlsx")
dataallg9<-na.omit(dataallg9)
dataallg9<- as.matrix(dataallg9)
dataallg9[dataallg9==123456789]<-NA
dataallg9<-na.omit(dataallg9)
dataallg9<- as.matrix(dataallg9) 
n9=length(dataallg9[,1])

dataallg10<-read_excel("10.xlsx")
dataallg10<-na.omit(dataallg10)
dataallg10<- as.matrix(dataallg10)
dataallg10[dataallg10==123456789]<-NA
dataallg10<-na.omit(dataallg10)
dataallg10<- as.matrix(dataallg10) 
n10=length(dataallg10[,1])

dataallg11<-read_excel("11.xlsx")
dataallg11<-na.omit(dataallg11)
dataallg11<- as.matrix(dataallg11)
dataallg11[dataallg11==123456789]<-NA
dataallg11<-na.omit(dataallg11)
dataallg11<- as.matrix(dataallg11) 
n11=length(dataallg11[,1])

dataallg12<-read_excel("12.xlsx")
dataallg12<-na.omit(dataallg12)
dataallg12<- as.matrix(dataallg12)
dataallg12[dataallg12==123456789]<-NA
dataallg12<-na.omit(dataallg12)
dataallg12<- as.matrix(dataallg12) 
n12=length(dataallg12[,1])

dataall<-rbind(dataallg1,dataallg2,dataallg3,dataallg4,dataallg5,dataallg6,dataallg7,dataallg8,dataallg9,dataallg10,dataallg11,dataallg12)

p=7
y1<-dataall[,1]
x2<-dataall[,2]
x3<-dataall[,3]
x4<-dataall[,4]
x5<-dataall[,5]
x6<-dataall[,6]
x7<-dataall[,7]
x8<-dataall[,8]


yall<-y1
xall2<-x2
xall3<-x3
xall4<-x4/100
xall5<-x5
xall6<-x6-mean(x6)
xall7<-x7-mean(x7)
xall8<-x8

xallg<-cbind(xall2,xall3,xall4,xall5,xall6,xall7,xall8)
xall=cbind(1,xallg)
N<-length(yall)
K=12


tall=c(0.1,0.3,0.5,0.7,0.9)
####################################################################
#######1- ball Chernozhukov et al. 2020. one-step ##################
####################################################################
#t1=Sys.time()
#b0=rq(yall~0+xall, tau = 0.1, method="fn")$coef
#taucall=seq(0.1,0.9,0.01)
#tl=length(taucall)
#cheab=matrix(0,nrow=tl,ncol=1)
#chemse=matrix(0,nrow=tl,ncol=1)
#b1=b0
#for (tii in 1:tl)
#{
#  ehat=yall-xall%*%b1
#  h=min(1.06*sd(ehat)*N^(-0.2),1)
# J2=t(xall)%*%(xall*as.vector(0.75*(1-ehat^2/h^2)*(abs(ehat/h)<=1))/N/h)
#  b2=b1-ginv(J2)%*%colSums(xall*as.vector(taucall[tii]-(ehat<=0)))/N
#  b1=b2
#  cheab[tii]=mean(abs(yall-xall%*%b1))
#  chemse[tii]=mean((yall-xall%*%b1)^2)
#}



t1=Sys.time()
b1=rq(yall~0+xall, tau = tall[1], method="fn")$coef
b2=rq(yall~0+xall, tau = tall[2], method="fn")$coef
b3=rq(yall~0+xall, tau = tall[3], method="fn")$coef
b4=rq(yall~0+xall, tau = tall[4], method="fn")$coef
b5=rq(yall~0+xall, tau = tall[5], method="fn")$coef
cheab=mean(abs(yall-xall%*%b3))
chemse=mean((yall-xall%*%b3)^2)
t2=Sys.time()
t_PK=t2-t1

########################################
#######2- ball  conquer#################
########################################
t3=Sys.time()
bh1=conquer(xallg, yall, tau = tall[1], kernel = "Gaussian")$coef
bh2=conquer(xallg, yall, tau = tall[2], kernel = "Gaussian")$coef
bh3=conquer(xallg, yall, tau = tall[3], kernel = "Gaussian")$coef
bh4=conquer(xallg, yall, tau = tall[4], kernel = "Gaussian")$coef
bh5=conquer(xallg, yall, tau = tall[5], kernel = "Gaussian")$coef
conquerab=mean(abs(yall-xall%*%bh3))
conquermse=mean((yall-xall%*%bh3)^2)
t4=Sys.time()
t_H=t4-t3


#####################################
#######3- ball  admm#################
#####################################
t5=Sys.time()
bg=QR.admm(xallg,yall,tall[1])
bp1=c(bg$b,bg$beta)

bg=QR.admm(xallg,yall,tall[2])
bp2=c(bg$b,bg$beta)

bg=QR.admm(xallg,yall,tall[3])
bp3=c(bg$b,bg$beta)

bg=QR.admm(xallg,yall,tall[4])
bp4=c(bg$b,bg$beta)

bg=QR.admm(xallg,yall,tall[5])
bp5=c(bg$b,bg$beta)

admmab=mean(abs(yall-xall%*%bp3))
admmmse=mean((yall-xall%*%bp3)^2)
t6=Sys.time()
t_P=t6-t5


############4-AVQR################
##################################

t7=Sys.time()
y1<-dataallg1[,1]
x2<-dataallg1[,2]
x3<-dataallg1[,3]
x4<-dataallg1[,4]/100
x5<-dataallg1[,5]
x6<-dataallg1[,6]-mean(dataallg1[,6])
x7<-dataallg1[,7]-mean(dataallg1[,7])
x8<-dataallg1[,8]
xallg1<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV11=rq(y1~xallg1, tau = tall[1])$coef
b_AV12=rq(y1~xallg1, tau = tall[2])$coef
b_AV13=rq(y1~xallg1, tau = tall[3])$coef
b_AV14=rq(y1~xallg1, tau = tall[4])$coef
b_AV15=rq(y1~xallg1, tau = tall[5])$coef


y2<-dataallg2[,1]
x2<-dataallg2[,2]
x3<-dataallg2[,3]
x4<-dataallg2[,4]/100
x5<-dataallg2[,5]
x6<-dataallg2[,6]-mean(dataallg2[,6])
x7<-dataallg2[,7]-mean(dataallg2[,7])
x8<-dataallg2[,8]
xallg2<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV21=rq(y2~xallg2, tau = tall[1])$coef
b_AV22=rq(y2~xallg2, tau = tall[2])$coef
b_AV23=rq(y2~xallg2, tau = tall[3])$coef
b_AV24=rq(y2~xallg2, tau = tall[4])$coef
b_AV25=rq(y2~xallg2, tau = tall[5])$coef


y3<-dataallg3[,1]
x2<-dataallg3[,2]
x3<-dataallg3[,3]
x4<-dataallg3[,4]/100
x5<-dataallg3[,5]
x6<-dataallg3[,6]-mean(dataallg3[,6])
x7<-dataallg3[,7]-mean(dataallg3[,7])
x8<-dataallg3[,8]
xallg3<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV31=rq(y3~xallg3, tau = tall[1])$coef
b_AV32=rq(y3~xallg3, tau = tall[2])$coef
b_AV33=rq(y3~xallg3, tau = tall[3])$coef
b_AV34=rq(y3~xallg3, tau = tall[4])$coef
b_AV35=rq(y3~xallg3, tau = tall[5])$coef


y4<-dataallg4[,1]
x2<-dataallg4[,2]
x3<-dataallg4[,3]
x4<-dataallg4[,4]/100
x5<-dataallg4[,5]
x6<-dataallg4[,6]-mean(dataallg4[,6])
x7<-dataallg4[,7]-mean(dataallg4[,7])
x8<-dataallg4[,8]
xallg4<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV41=rq(y4~xallg4, tau = tall[1])$coef
b_AV42=rq(y4~xallg4, tau = tall[2])$coef
b_AV43=rq(y4~xallg4, tau = tall[3])$coef
b_AV44=rq(y4~xallg4, tau = tall[4])$coef
b_AV45=rq(y4~xallg4, tau = tall[5])$coef


y5<-dataallg5[,1]
x2<-dataallg5[,2]
x3<-dataallg5[,3]
x4<-dataallg5[,4]/100
x5<-dataallg5[,5]
x6<-dataallg5[,6]-mean(dataallg5[,6])
x7<-dataallg5[,7]-mean(dataallg5[,7])
x8<-dataallg5[,8]
xallg5<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV51=rq(y5~xallg5, tau = tall[1])$coef
b_AV52=rq(y5~xallg5, tau = tall[2])$coef
b_AV53=rq(y5~xallg5, tau = tall[3])$coef
b_AV54=rq(y5~xallg5, tau = tall[4])$coef
b_AV55=rq(y5~xallg5, tau = tall[5])$coef


y6<-dataallg6[,1]
x2<-dataallg6[,2]
x3<-dataallg6[,3]
x4<-dataallg6[,4]/100
x5<-dataallg6[,5]
x6<-dataallg6[,6]-mean(dataallg6[,6])
x7<-dataallg6[,7]-mean(dataallg6[,7])
x8<-dataallg6[,8]
xallg6<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV61=rq(y6~xallg6, tau = tall[1])$coef
b_AV62=rq(y6~xallg6, tau = tall[2])$coef
b_AV63=rq(y6~xallg6, tau = tall[3])$coef
b_AV64=rq(y6~xallg6, tau = tall[4])$coef
b_AV65=rq(y6~xallg6, tau = tall[5])$coef


y7<-dataallg7[,1]
x2<-dataallg7[,2]
x3<-dataallg7[,3]
x4<-dataallg7[,4]/100
x5<-dataallg7[,5]
x6<-dataallg7[,6]-mean(dataallg7[,6])
x7<-dataallg7[,7]-mean(dataallg7[,7])
x8<-dataallg7[,8]
xallg7<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV71=rq(y7~xallg7, tau = tall[1])$coef
b_AV72=rq(y7~xallg7, tau = tall[2])$coef
b_AV73=rq(y7~xallg7, tau = tall[3])$coef
b_AV74=rq(y7~xallg7, tau = tall[4])$coef
b_AV75=rq(y7~xallg7, tau = tall[5])$coef


y8<-dataallg8[,1]
x2<-dataallg8[,2]
x3<-dataallg8[,3]
x4<-dataallg8[,4]/100
x5<-dataallg8[,5]
x6<-dataallg8[,6]-mean(dataallg8[,6])
x7<-dataallg8[,7]-mean(dataallg8[,7])
x8<-dataallg8[,8]
xallg8<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV81=rq(y8~xallg8, tau = tall[1])$coef
b_AV82=rq(y8~xallg8, tau = tall[2])$coef
b_AV83=rq(y8~xallg8, tau = tall[3])$coef
b_AV84=rq(y8~xallg8, tau = tall[4])$coef
b_AV85=rq(y8~xallg8, tau = tall[5])$coef


y9<-dataallg9[,1]
x2<-dataallg9[,2]
x3<-dataallg9[,3]
x4<-dataallg9[,4]/100
x5<-dataallg9[,5]
x6<-dataallg9[,6]-mean(dataallg9[,6])
x7<-dataallg9[,7]-mean(dataallg9[,7])
x8<-dataallg9[,8]
xallg9<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV91=rq(y9~xallg9, tau = tall[1])$coef
b_AV92=rq(y9~xallg9, tau = tall[2])$coef
b_AV93=rq(y9~xallg9, tau = tall[3])$coef
b_AV94=rq(y9~xallg9, tau = tall[4])$coef
b_AV95=rq(y9~xallg9, tau = tall[5])$coef


y10<-dataallg10[,1]
x2<-dataallg10[,2]
x3<-dataallg10[,3]
x4<-dataallg10[,4]/100
x5<-dataallg10[,5]
x6<-dataallg10[,6]-mean(dataallg10[,6])
x7<-dataallg10[,7]-mean(dataallg10[,7])
x8<-dataallg10[,8]
xallg10<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV101=rq(y10~xallg10, tau = tall[1])$coef
b_AV102=rq(y10~xallg10, tau = tall[2])$coef
b_AV103=rq(y10~xallg10, tau = tall[3])$coef
b_AV104=rq(y10~xallg10, tau = tall[4])$coef
b_AV105=rq(y10~xallg10, tau = tall[5])$coef


y11<-dataallg11[,1]
x2<-dataallg11[,2]
x3<-dataallg11[,3]
x4<-dataallg11[,4]/100
x5<-dataallg11[,5]
x6<-dataallg11[,6]-mean(dataallg11[,6])
x7<-dataallg11[,7]-mean(dataallg11[,7])
x8<-dataallg11[,8]
xallg11<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV111=rq(y11~xallg11, tau = tall[1])$coef
b_AV112=rq(y11~xallg11, tau = tall[2])$coef
b_AV113=rq(y11~xallg11, tau = tall[3])$coef
b_AV114=rq(y11~xallg11, tau = tall[4])$coef
b_AV115=rq(y11~xallg11, tau = tall[5])$coef


y12<-dataallg12[,1]
x2<-dataallg12[,2]
x3<-dataallg12[,3]
x4<-dataallg12[,4]/100
x5<-dataallg12[,5]
x6<-dataallg12[,6]-mean(dataallg12[,6])
x7<-dataallg12[,7]-mean(dataallg12[,7])
x8<-dataallg12[,8]
xallg12<-cbind(x2,x3,x4,x5,x6,x7,x8)
b_AV121=rq(y12~xallg12, tau = tall[1])$coef
b_AV122=rq(y12~xallg12, tau = tall[2])$coef
b_AV123=rq(y12~xallg12, tau = tall[3])$coef
b_AV124=rq(y12~xallg12, tau = tall[4])$coef
b_AV125=rq(y12~xallg12, tau = tall[5])$coef


b_AV1=(n1*b_AV11+n2*b_AV21+n3*b_AV31+n4*b_AV41+n5*b_AV51+n6*b_AV61+n7*b_AV71+n8*b_AV81+n9*b_AV91+n10*b_AV101+n11*b_AV111+n12*b_AV121)/N
b_AV2=(n1*b_AV12+n2*b_AV22+n3*b_AV32+n4*b_AV42+n5*b_AV52+n6*b_AV62+n7*b_AV72+n8*b_AV82+n9*b_AV92+n10*b_AV102+n11*b_AV112+n12*b_AV122)/N
b_AV3=(n1*b_AV13+n2*b_AV23+n3*b_AV33+n4*b_AV43+n5*b_AV53+n6*b_AV63+n7*b_AV73+n8*b_AV83+n9*b_AV93+n10*b_AV103+n11*b_AV113+n12*b_AV123)/N
b_AV4=(n1*b_AV14+n2*b_AV24+n3*b_AV34+n4*b_AV44+n5*b_AV54+n6*b_AV64+n7*b_AV74+n8*b_AV84+n9*b_AV94+n10*b_AV104+n11*b_AV114+n12*b_AV124)/N
b_AV5=(n1*b_AV15+n2*b_AV25+n3*b_AV35+n4*b_AV45+n5*b_AV55+n6*b_AV65+n7*b_AV75+n8*b_AV85+n9*b_AV95+n10*b_AV105+n11*b_AV115+n12*b_AV125)/N
AVab=mean(abs(yall-xall%*%b_AV3))
AVse=mean((yall-xall%*%b_AV3)^2)

t8=Sys.time()
t_AV=t8-t7

############5-AVCQR################
##################################
t9=Sys.time()
taucqr=1:5/6
bcqr1=cqr.admm(xallg1,y1,taucqr)$beta
bcqr2=cqr.admm(xallg2,y2,taucqr)$beta
bcqr3=cqr.admm(xallg3,y3,taucqr)$beta
bcqr4=cqr.admm(xallg4,y4,taucqr)$beta
bcqr5=cqr.admm(xallg5,y5,taucqr)$beta
bcqr6=cqr.admm(xallg6,y6,taucqr)$beta
bcqr7=cqr.admm(xallg7,y7,taucqr)$beta
bcqr8=cqr.admm(xallg8,y8,taucqr)$beta
bcqr9=cqr.admm(xallg9,y9,taucqr)$beta
bcqr10=cqr.admm(xallg10,y10,taucqr)$beta
bcqr11=cqr.admm(xallg11,y11,taucqr)$beta
bcqr12=cqr.admm(xallg12,y12,taucqr)$beta
xx=t(xallg1)%*%xallg1+t(xallg2)%*%xallg2+t(xallg3)%*%xallg3+t(xallg4)%*%xallg4+t(xallg5)%*%xallg5+t(xallg6)%*%xallg6+t(xallg7)%*%xallg7+t(xallg8)%*%xallg8+t(xallg9)%*%xallg9+t(xallg10)%*%xallg10+t(xallg11)%*%xallg11+t(xallg12)%*%xallg12
xxb=t(xallg1)%*%xallg1%*%bcqr1+t(xallg2)%*%xallg2%*%bcqr2+t(xallg3)%*%xallg3%*%bcqr3+t(xallg4)%*%xallg4%*%bcqr4+t(xallg5)%*%xallg5%*%bcqr5+t(xallg6)%*%xallg6%*%bcqr6+t(xallg7)%*%xallg7%*%bcqr7+t(xallg8)%*%xallg8%*%bcqr8+t(xallg9)%*%xallg9%*%bcqr9+t(xallg10)%*%xallg10%*%bcqr10+t(xallg11)%*%xallg11%*%bcqr11+t(xallg12)%*%xallg12%*%bcqr12
b_cqr=ginv(xx)%*%xxb

cqrab=mean(abs(yall-xallg%*%b_cqr))
cqrse=mean((yall-xallg%*%b_cqr)^2)
t10=Sys.time()
t_cqr=t10-t9

############6-Volgushev et al 2019################
##################################################
t11=Sys.time()
Vol=function(tau,b)
{
  p=length(b[,1])
  Bx=function(x)
  {
    tall=seq(0.1,0.9,0.2)
    Bx=matrix(c(1,x,x^2,(x-tall[1])^2,(x-tall[2])^2,(x-tall[3])^2,(x-tall[4])^2,(x-tall[5])^2),nrow=1,ncol=8)
    return(Bx)
  }
  tall=seq(0.1,0.9,0.2)
  
  B1=Bx(tall[1])
  B2=Bx(tall[2])
  B3=Bx(tall[3])
  B4=Bx(tall[4])
  B5=Bx(tall[5])

  
  a1=t(B1)%*%B1
  a2=t(B2)%*%B2
  a3=t(B3)%*%B3
  a4=t(B4)%*%B4
  a5=t(B5)%*%B5

  BB=a1+a2+a3+a4+a5
  B=cbind(t(B1),t(B2),t(B3),t(B4),t(B5))
  b0=matrix(1:5,nrow=5,ncol=1)
  
  aall=matrix(0,nrow=8,ncol=p)
  for (j in 1:p)
  {aall[,j]=ginv(BB)%*%(B%*%b[j,])}
  
  Bxx=as.vector(Bx(tau))
  hatb=t(aall)%*%Bxx
  return(hatb)
}
b=cbind(b_AV1,b_AV2,b_AV3,b_AV4,b_AV5)
b_vol1=Vol(tall[1],b)
b_vol2=Vol(tall[2],b)
b_vol3=Vol(tall[3],b)
b_vol4=Vol(tall[4],b)
b_vol5=Vol(tall[5],b)
volab=mean(abs(yall-xall%*%b_vol3))
volse=mean((yall-xall%*%b_vol3)^2)
t12=Sys.time()
t_vol=t12-t11+t_AV


##################################
#######7- b0  ###################
##################################
b_LE=matrix(0,nrow=8,ncol=5)
t_LE=matrix(0,nrow=5,ncol=1)
b_DS=matrix(0,nrow=8,ncol=5)
t_DS=matrix(0,nrow=5,ncol=1)
for (tii in 1:5)
{
  tau=tall[tii]
  
t13=Sys.time()
y0=yall[1:n1]
x0=xall[1:n1,]  
b0=rq(y0~0+x0, tau = tau)$coef
t14=Sys.time()
t_b0=t14-t13
##################################
#######end 2-b0  ###############
##################################

##################################
#######8- LEQR  ##################
##################################
t15=Sys.time()
n=n1
M=ceiling(2+log(log(sqrt((p+1)/N))/log((p+1)/n))/log(2)) ###迭代数
blqrg=b0
for(mi in 1:M)
{
    hg1=max(sqrt((p+1)/N),((p+1)/n)^(2^(mi-2)))
    sg1=(yall-xall%*%blqrg)/hg1
    S1=(sg1>=1)+(0.5+15/16*(sg1-2/3*sg1^3+0.2*sg1^5))*(abs(sg1)<1)
    Sp1=15/16*(1-2*sg1^2+sg1^4)*(abs(sg1)<1)
    g11=S1+tau-1+yall/hg1*Sp1
    Ul=t(xall)%*%g11
    
    txg11l=matrix(0,nrow=N,ncol=(p+1))
    for (pi in 1:(p+1)) 
    {
      txg11l[,pi]=xall[,pi]*Sp1
    }
    Vl=t(xall)%*%txg11l/hg1
  blqrg=ginv(Vl)%*%Ul
}
b_LE[,tii]=blqrg
t16=Sys.time()
t_LE[tii]=t16-t15+t_b0
##################################
#######END- LEQR  ################
##################################  

##################################
#######9- DSQR  ##################
##################################
t17=Sys.time()
Q=ceiling(1+1.5*log(n1/N)/log((p+1)*log(n1)/n1)) ###迭代数
blqrg=b0
h=((p+1)/N)^{0.4}
h1=((p+1)*log(n)/n)^{1/3}
S2=matrix(0,nrow=(p+1),ncol=(p+1))
txg=matrix(0,nrow=n1,ncol=(p+1))
for(qi in 1:Q)
{
  ag=as.vector(dnorm((y0-x0%*%blqrg)/h1))
  for (pi in 1:(p+1)) 
  {
    txg[,pi]=x0[,pi]*ag
  }
  S2=t(x0)%*%txg/n1/h1
  S1=t(xall)%*%(pnorm(-(yall-xall%*%blqrg)/h)-tau)/N
  blqrg=blqrg-ginv(S2)%*%S1
}
b_DS[,tii]=blqrg
t18=Sys.time()
t_DS[tii]=t18-t17+t_b0
}

LElab=mean(abs(yall-xall%*%b_LE[,3]))
LElse=mean((yall-xall%*%b_LE[,3])^2)
t_LEall=sum(t_LE)

DSlab=mean(abs(yall-xall%*%b_DS[,3]))
DSlse=mean((yall-xall%*%b_DS[,3])^2)
t_DSall=sum(t_DS)
##################################
#######5- DSQR  ##################
##################################
t19=Sys.time()
bls=lm(yall ~ xallg)$coef
LSlab=mean(abs(yall-xall%*%bls))
t20=Sys.time()
t_ls=t20-t19
