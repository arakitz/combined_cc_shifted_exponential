### SE n=1, R:2/3
theta0<-10
lambda0<-100
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)}
a1<-0.0015
theta0-lambda0*log(0.5)
LCLsh<-theta0-lambda0*log(1-a1)
LWL<-12.031
CL<-theta0-lambda0*log(0.5)
p00<-1-F0X(CL)
p01<-(-F0X(LWL)+F0X(CL))
p02<-(-F0X(LCLsh)+F0X(LWL))
r1<-c(p00+p01,p02,0)
r2<-c(p00,0,p01)
r3<-c(p00+p01,0,0)
Q0<-rbind(r1,r2,r3)
e1<-c(1,0,0)
ID<-diag(3)
l1<-c(1,1,1)
ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1)
ARL0
########### OOC #####################
a1<-0.0005
theta0-lambda0*log(0.5)
LCLsh<-theta0-lambda0*log(1-a1)
LWL<-13.379
for(theta1 in c(10,7,5)){
  for(lambda1 in c(100,90,80,70,60,50)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
#    a1<-0.0015
#    theta0-lambda0*log(0.5)
#    LCLsh<-theta0-lambda0*log(1-a1)
#    LWL<-12.031
    CL<-theta0-lambda0*log(0.5)
    p10<-1-F1X(CL)
    p11<-(-F1X(LWL)+F1X(CL))
    p12<-(-F1X(LCLsh)+F1X(LWL))
    r1<-c(p10+p11,p12,0)
    r2<-c(p10,0,p11)
    r3<-c(p10+p11,0,0)
    Q1<-rbind(r1,r2,r3)
    e1<-c(1,0,0)
    ID<-diag(3)
    l1<-c(1,1,1)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1)
    print(c(a1,LCLsh,LWL,theta1,lambda1,ARL1))
  }
}
#############################
library(cubature)
###############################
theta0<-10 # IC theta
lambda0<-100 # IC lambda
thetamin<-7
thetamax<-10
lambdamin<-70
lambdamax<-100
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
  LCLsh<-theta0-lambda0*log(1-a1)
  LWL<-13.379
  
  p10<-1-F1X(CL)
  p11<-(-F1X(LWL)+F1X(CL))
  p12<-(-F1X(LCLsh)+F1X(LWL))
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
  LCLsh<-theta0-lambda0*log(1-a1)
  LWL<-13.379
  
  p10<-1-F1X(CL)
  p11<-(-F1X(LWL)+F1X(CL))
  p12<-(-F1X(LCLsh)+F1X(LWL))
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