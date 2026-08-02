# Below we provide the commands to find the 
# the pair of (UCLsh, UWL) for the SH-RR:4/5 chart.
# Also we calculate the EARL, EQARL and the OOC ARL performance
################################################################
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # CDF of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, for calculating UCLsh
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
# manual numerical search to find the UWL, below it is the final stage of the search
for(UWL in seq(183.16,183.17,by=0.001)){
  CL<-theta0-lambda0*log(0.5) # median of the SE(theta0,lambda0), center line
  # transition probability limits
  p00<-F0X(CL)
  p01<-(F0X(UWL)-F0X(CL))
  p02<-(F0X(UCLsh)-F0X(UWL))
  r1<-c(p00+p01,p02,0,0,0,0,0,0,0,0)
  r2<-c(p00,0,p02,0,p01,0,0,0,0,0)
  r3<-c(p00,0,0,p02,0,0,0,p01,0,0)
  r4<-c(p00,0,0,0,0,0,0,0,0,p01)
  r5<-c(p00+p01,0,0,0,0,p02,0,0,0,0)
  r6<-c(p00,0,0,0,p01,0,p02,0,0,0)
  r7<-c(p00,0,0,0,0,0,0,p01,0,0)
  r8<-c(p00+p01,0,0,0,0,0,0,0,p02,0)
  r9<-c(p00,0,0,0,p01,0,0,0,0,0)
  r10<-c(p00+p01,0,0,0,0,0,0,0,0,0)
  Q0<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
  e1<-rep(0,10);e1[1]<-1
  ID<-diag(10)
  l1<-rep(1,10)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # calculate the IC ARL using the Markov chain method
  print(c(UCLsh,UWL,ARL0)) # print the Shewhart limit, the warning limit and the IC ARL
}

############# EARL & EQARL ###############
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
UWL<-183.17 # Warning limit, user must provide the correct value found above
CL<-theta0-lambda0*log(0.5) # center line of the chart
listARL<-c() # empty vector to store the ARLs
listQARL<-c() # empty vector to store quantities necessary to calculate the EQARL
# define the interval of shifts in theta
thetamin<-10
thetamax<-13
# define the interval of shifts in lambda
lambdamin<-100
lambdamax<-130
for(theta1 in seq(thetamin,thetamax,length=10)){
  for(lambda1 in seq(lambdamin,lambdamax,length=10)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # CDF of SE(theta1,lambda1)
    p10<-F1X(CL)
    p11<-(F1X(UWL)-F1X(CL))
    p12<-(F1X(UCLsh)-F1X(UWL))
    r1<-c(p10+p11,p12,0,0,0,0,0,0,0,0)
    r2<-c(p10,0,p12,0,p11,0,0,0,0,0)
    r3<-c(p10,0,0,p12,0,0,0,p11,0,0)
    r4<-c(p10,0,0,0,0,0,0,0,0,p11)
    r5<-c(p10+p11,0,0,0,0,p12,0,0,0,0)
    r6<-c(p10,0,0,0,p11,0,p12,0,0,0)
    r7<-c(p10,0,0,0,0,0,0,p11,0,0)
    r8<-c(p10+p11,0,0,0,0,0,0,0,p12,0)
    r9<-c(p10,0,0,0,p11,0,0,0,0,0)
    r10<-c(p10+p11,0,0,0,0,0,0,0,0,0)
    Q1<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) # one-step transition probabilities matrix
    e1<-rep(0,10);e1[1]<-1
    ID<-diag(10)
    l1<-rep(1,10)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # ARL calculation uwing Markov chain method
    QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
    listARL<-c(listARL,ARL1)
    listQARL<-c(listQARL,QARL1)
  }
}
EARL<-mean(listARL)
EQARL<-mean(listQARL)
print(round(c(a1,UCLsh,UWL,theta0,lambda0,EARL,EQARL),digits=5))
############ OOC ARL #################################
for(theta1 in c(10,13,15)){
  for(lambda1 in c(100,120,140,160,180,200)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
    p10<-F1X(CL)
    p11<-(F1X(UWL)-F1X(CL))
    p12<-(F1X(UCLsh)-F1X(UWL))
    r1<-c(p10+p11,p12,0,0,0,0,0,0,0,0)
    r2<-c(p10,0,p12,0,p11,0,0,0,0,0)
    r3<-c(p10,0,0,p12,0,0,0,p11,0,0)
    r4<-c(p10,0,0,0,0,0,0,0,0,p11)
    r5<-c(p10+p11,0,0,0,0,p12,0,0,0,0)
    r6<-c(p10,0,0,0,p11,0,p12,0,0,0)
    r7<-c(p10,0,0,0,0,0,0,p11,0,0)
    r8<-c(p10+p11,0,0,0,0,0,0,0,p12,0)
    r9<-c(p10,0,0,0,p11,0,0,0,0,0)
    r10<-c(p10+p11,0,0,0,0,0,0,0,0,0)
    Q1<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
    e1<-rep(0,10);e1[1]<-1
    ID<-diag(10)
    l1<-rep(1,10)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1)
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1))
  }
}