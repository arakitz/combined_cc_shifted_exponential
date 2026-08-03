### SE n=1, R:8/8
theta0<-0
lambda0<-1
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)}
a1<-0.0005
theta0-lambda0*log(0.5)
UCLsh<-theta0-lambda0*log(a1)
#UWL<-63.548
for(UWL in seq(83.00,83.01,by=0.001)){
  p00<-F0X(UWL)
  p01<-(F0X(UCLsh)-F0X(UWL))
  Q0<-matrix(0,ncol=8,nrow=8)
  Q0[,1]<-rep(p00,8)
  for(i in 1:7){Q0[i,i+1]<-p01}
  e1<-rep(0,8);e1[1]<-1
  ID<-diag(8)
  l1<-rep(1,8)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1)
  print(c(UCLsh,UWL,ARL0))  
}

########## OOC performance ##################
a1<-0.0005
UCLsh<-theta0-lambda0*log(a1)
UWL<-0.730
for(theta1 in c(0,0.25,0.5)){
  for(lambda1 in c(1,1.25,1.5,1.75,2,3)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
    
#    theta0-lambda0*log(0.5)
    p10<-F1X(UWL)
    p11<-(F1X(UCLsh)-F1X(UWL))
    Q1<-matrix(0,ncol=8,nrow=8)
    Q1[,1]<-rep(p10,8)
    for(i in 1:7){Q1[i,i+1]<-p11}
    e1<-rep(0,8);e1[1]<-1
    ID<-diag(8)
    l1<-rep(1,8)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1)
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1))
  }
}
######################################################
library(cubature)
###############################
theta0<-0 # IC theta
lambda0<-1 # IC lambda
thetamin<-0
thetamax<-1.5
lambdamin<-1
lambdamax<-2
ARLearl<-function(z){
  z1<-z[1]
  z2<-z[2]
  theta1<-z1
  lambda1<-z2
  F0X<-function(x){
    (x>=theta0)*(1-exp(-(x-theta0)/lambda0))
  }
  # OOC case
  F1X<-function(x){
    (x>=theta1)*(1-exp(-(x-theta1)/lambda1))
  }
  
  a1<-0.0005
  MDL<-theta0-lambda0*log(0.5)
  UCLsh<-theta0-lambda0*log(a1)
  UWL<-0.730

  p10<-F1X(UWL)
  p11<-(F1X(UCLsh)-F1X(UWL))
  Q1<-matrix(0,ncol=8,nrow=8)
  Q1[,1]<-rep(p10,8)
  for(i in 1:7){Q1[i,i+1]<-p11}
  e1<-rep(0,8);e1[1]<-1
  Id<-diag(8)
  l1<-rep(1,8)
  
  ARL1<-as.vector(e1%*%solve(Id-Q1)%*%l1) # ARL
  return(ARL1)
}
#######################
VALUE <- hcubature(ARLearl,lowerLimit=c(thetamin,lambdamin),
                   upperLimit=c(thetamax,lambdamax))$integral   
EARL1 <- VALUE/((thetamax-thetamin)*(lambdamax-lambdamin))
EARL1 
##########################################
QARLearl<-function(z){
  z1<-z[1]
  z2<-z[2]
  theta1<-z1
  lambda1<-z2
  F0X<-function(x){
    (x>=theta0)*(1-exp(-(x-theta0)/lambda0))
  }
  # OOC case
  F1X<-function(x){
    (x>=theta1)*(1-exp(-(x-theta1)/lambda1))
  }

  a1<-0.0005
  MDL<-theta0-lambda0*log(0.5)
  UCLsh<-theta0-lambda0*log(a1)
  UWL<-0.730
  
  p10<-F1X(UWL)
  p11<-(F1X(UCLsh)-F1X(UWL))
  Q1<-matrix(0,ncol=8,nrow=8)
  Q1[,1]<-rep(p10,8)
  for(i in 1:7){Q1[i,i+1]<-p11}
  e1<-rep(0,8);e1[1]<-1
  Id<-diag(8)
  l1<-rep(1,8)
  
  ARL1<-as.vector(e1%*%solve(Id-Q1)%*%l1) # ARL
#  return(ARL1)
  QARL1<-ARL1*(1/(1+(theta1-theta0)^2))*(1/(1+(lambda1-lambda0)^2))  
  return(QARL1)

}
#######################
QVALUE <- hcubature(QARLearl,lowerLimit=c(thetamin,lambdamin),
                    upperLimit=c(thetamax,lambdamax))$integral   
EQARL1 <- QVALUE/((thetamax-thetamin)*(lambdamax-lambdamin))
EQARL1 
##########################################
c(EARL1,EQARL1) 