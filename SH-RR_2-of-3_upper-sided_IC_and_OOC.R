### SE n=1, R:2/3
theta0<-0
lambda0<-1
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)}
a1<-0.0005
theta0-lambda0*log(0.5)
UCLsh<-theta0-lambda0*log(a1)
#UWL<-350
for(UWL in seq(350.4,350.41,by=0.001)){
  CL<-theta0-lambda0*log(0.5)
  p00<-F0X(CL)
  p01<-(F0X(UWL)-F0X(CL))
  p02<-(F0X(UCLsh)-F0X(UWL))
  r1<-c(p00+p01,p02,0)
  r2<-c(p00,0,p01)
  r3<-c(p00+p01,0,0)
  Q0<-rbind(r1,r2,r3)
  e1<-c(1,0,0)
  ID<-diag(3)
  l1<-c(1,1,1)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1)
  print(c(UCLsh,UWL,ARL0))  
}

########### OOC #####################
a1<-0.0005
theta0-lambda0*log(0.5)
UCLsh<-theta0-lambda0*log(a1)
UWL<-3.405
for(theta1 in c(0,0.25,0.5)){
  for(lambda1 in c(1,1.25,1.5,1.75,2,3)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
#    a1<-0.0015
#    theta0-lambda0*log(0.5)
#    LCLsh<-theta0-lambda0*log(1-a1)
#    LWL<-12.031
    CL<-theta0-lambda0*log(0.5)
    p10<-F1X(CL)
    p11<-(F1X(UWL)-F1X(CL))
    p12<-(F1X(UCLsh)-F1X(UWL))
    r1<-c(p10+p11,p12,0)
    r2<-c(p10,0,p11)
    r3<-c(p10+p11,0,0)
    Q1<-rbind(r1,r2,r3)
    e1<-c(1,0,0)
    ID<-diag(3)
    l1<-c(1,1,1)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1)
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1))
  }
}
#############################
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
  CL<-theta0-lambda0*log(0.5)
  UCLsh<-theta0-lambda0*log(a1)
  UWL<-3.405
  
  p10<-F1X(CL)
  p11<-(F1X(UWL)-F1X(CL))
  p12<-(F1X(UCLsh)-F1X(UWL))
  r1<-c(p10+p11,p12,0)
  r2<-c(p10,0,p11)
  r3<-c(p10+p11,0,0)
  Q1<-rbind(r1,r2,r3)
  e1<-c(1,0,0)
  ID<-diag(3)
  l1<-c(1,1,1)
  ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # ARL
  
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
  CL<-theta0-lambda0*log(0.5)
  UCLsh<-theta0-lambda0*log(a1)
  UWL<-3.405
  
  p10<-F1X(CL)
  p11<-(F1X(UWL)-F1X(CL))
  p12<-(F1X(UCLsh)-F1X(UWL))
  r1<-c(p10+p11,p12,0)
  r2<-c(p10,0,p11)
  r3<-c(p10+p11,0,0)
  Q1<-rbind(r1,r2,r3)
  e1<-c(1,0,0)
  ID<-diag(3)
  l1<-c(1,1,1)
  ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # ARL
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
