# Below we provide the commands that can be used to determine
# the UCLsh and UWL for the SH-RR:8/8 upper one-sided chart
# Also, we calculate the OOC ARL performance, as well as the EARL, EQARL
#####################################
theta0<-10 # IC origin
lambda0<-100 # OOC origin
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # CDF of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, needed to determine the Shewhart limit
UCLsh<-theta0-lambda0*log(a1)
# manual numerical search to find the UWL, below is the final part of the searching procedure
for(UWL in seq(83.00,83.01,by=0.001)){
  p00<-F0X(UWL)
  p01<-(F0X(UCLsh)-F0X(UWL))
  Q0<-matrix(0,ncol=8,nrow=8)
  Q0[,1]<-rep(p00,8)
  for(i in 1:7){Q0[i,i+1]<-p01}
  e1<-rep(0,8);e1[1]<-1
  ID<-diag(8)
  l1<-rep(1,8)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the ARL is calculated by using the Markov chain method
  print(c(UCLsh,UWL,ARL0))  # print the Shewhart limit, the warning limit and the IC ARL
}

########## OOC performance ##################
UCLsh<-theta0-lambda0*log(a1)
UWL<-83.004 # the UWL must be provided by the user and it is the value obtained previously
listARL<-c() # an empty vector to store the ARLs
listQARL<-c() # an empty vector to store the necessary quantities to calculate the EQARL
# define the interval of the shifts in theta
thetamin<-10
thetamax<-13
# define the interval of the shifts in lambda
lambdamin<-100
lambdamax<-130
for(theta1 in seq(thetamin,thetamax,length=10)){
  for(lambda1 in seq(lambdamin,lambdamax,length=10)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # CDF of the SE(theta1,lambda1)
    p10<-F1X(UWL)
    p11<-(F1X(UCLsh)-F1X(UWL))
    Q1<-matrix(0,ncol=8,nrow=8)
    Q1[,1]<-rep(p10,8)
    for(i in 1:7){Q1[i,i+1]<-p11}
    e1<-rep(0,8);e1[1]<-1
    ID<-diag(8)
    l1<-rep(1,8)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # OOC ARL, calculated by using the Markov chain method
    QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
    listARL<-c(listARL,ARL1)
    listQARL<-c(listQARL,QARL1)
  }
}
EARL<-mean(listARL)
EQARL<-mean(listQARL)
print(round(c(a1,UCLsh,UWL,theta0,lambda0,EARL,EQARL),digits=5))
############### OOC ARL #######################################
for(theta1 in c(10,13,15)){
  for(lambda1 in c(100,120,140,160,180,200)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)}
    p10<-F1X(UWL)
    p11<-(F1X(UCLsh)-F1X(UWL))
    Q1<-matrix(0,ncol=8,nrow=8)
    Q1[,1]<-rep(p10,8)
    for(i in 1:7){Q1[i,i+1]<-p11}
    e1<-rep(0,8);e1[1]<-1
    ID<-diag(8)
    l1<-rep(1,8)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1)
    print(c(UCLsh,UWL,theta0,lambda0,theta1,lambda1,ARL1))
  }
}