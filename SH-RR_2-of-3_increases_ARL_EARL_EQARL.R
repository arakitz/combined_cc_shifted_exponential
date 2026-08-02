# Below we provide the commands that calculate the
# control and warning limits of the upper one-sided SH-RR:2/3 chart
###############################
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # cdf of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, for the UCLsh
UCLsh<-theta0-lambda0*log(a1) # Shewhart control limit
# numerical search to find the UWL, here is the final search
for(UWL in seq(350.4,350.41,by=0.001)){
  CL<-theta0-lambda0*log(0.5) # median of the SE(theta0,lambda0), as center line
  # transition probabilities
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
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # calculate ARL using Markov chain method
  print(c(UCLsh,UWL,ARL0))  # prints the Shewhart limit, the warning limit and the IC ARL
}

########### OOC #####################
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
UWL<-350.407 # Warning limit, user must give the correct  value, obtained above
listARL<-c() # empty vector to store ARLs
listQARL<-c() # empty vector to store quantities needed to calculate EQARL
# define the interval for the shifts in theta
thetamin<-10
thetamax<-13
# define the interval for the shifts in lambda
lambdamin<-100
lambdamax<-130
for(theta1 in seq(thetamin,thetamax,length=10)){
  for(lambda1 in seq(lambdamin,lambdamax,length=10)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # OOC CDF of SE(theta1,lambda1)
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
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # calculate OOC ARL using the Markov chain method
    QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
    #    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1))
    listARL<-c(listARL,ARL1)
    listQARL<-c(listQARL,QARL1)
  }
}
EARL<-mean(listARL)
EQARL<-mean(listQARL)
print(round(c(a1,UCLsh,UWL,theta0,lambda0,EARL,EQARL),digits=5))
########### OOC ARL ##################
for(theta1 in c(10,13,15)){
  for(lambda1 in c(100,120,140,160,180,200)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
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