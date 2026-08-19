# below we provide the commands needed to calculate
# the IC and the OOC ARL performance of the SH-RR:8/8 upper one-sided chart
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # IC cdf of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, needed to determine the upper Shewhart limit
UCLsh<-theta0-lambda0*log(a1) # the upper Shewhart limit
# below we determine, manually, the UWL so as the IC ARL, given the upper Shewhart limit
# results in IC ARL equal to the nominal value
# here, you see the final stage of the procedure
for(UWL in seq(83.00,83.01,by=0.001)){
  # below we determine the one-step transition probabilities
  # we calculate the ARL by using the Markov chain method
  p00<-F0X(UWL)
  p01<-(F0X(UCLsh)-F0X(UWL))
  Q0<-matrix(0,ncol=8,nrow=8)
  Q0[,1]<-rep(p00,8)
  for(i in 1:7){Q0[i,i+1]<-p01}
  e1<-rep(0,8);e1[1]<-1
  ID<-diag(8)
  l1<-rep(1,8)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the IC ARL
  print(c(UCLsh,UWL,ARL0))  # print the results to verify if the pair (UWL,UCLsh) gives the desired IC performance
}
########## OOC ARL performance ##################
UWL<-83.004 # here the user has to provide the UWL obtained earlier
# below we calculate the OOC ARL for various OOC values theta1 and lambda1
for(theta1 in c(10,13,15)){
  for(lambda1 in c(100,120,140,160,180,200)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # OOC cdf of the SE(theta1,lambda1)
    # below we determine the one-step transition probabilities
    # we calculate the ARL by using the Markov chain method
    p10<-F1X(UWL)
    p11<-(F1X(UCLsh)-F1X(UWL))
    Q1<-matrix(0,ncol=8,nrow=8)
    Q1[,1]<-rep(p10,8)
    for(i in 1:7){Q1[i,i+1]<-p11}
    e1<-rep(0,8);e1[1]<-1
    ID<-diag(8)
    l1<-rep(1,8)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # the OOC ARL
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1)) # print the results
  }
}
#################### END ##################################
