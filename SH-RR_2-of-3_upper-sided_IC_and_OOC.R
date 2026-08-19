# Calculation of IC and OOC ARL of the
# SH-RR:2/3 upper one-sided chart
theta0<-0 # IC origin
lambda0<-1 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # IC cdf of SE(theta0,lambda0)
a1<-0.0005 # alpha1, needed to determine the Shewhart limit
UCLsh<-theta0-lambda0*log(a1) # the Shewhart upper limit
# below we determine manually the upper warning limit UWL
# this is the final stage of determination
for(UWL in seq(350.4,350.41,by=0.001)){
  CL<-theta0-lambda0*log(0.5) # center line of the chart
  # below we construnct the one-step transition probabilities matrix
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
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # we calculate ARL by using the Markov chain method
  print(c(UCLsh,UWL,ARL0)) # print the Shewhart limit, the UWL and the ARL0  
}
########### OOC ARL performance #####################
UWL<-3.405 # give the UWL value obtained in the previous step
# calculate the OOC ARL values for various OOC values theta1 and lambda1
for(theta1 in c(0,0.25,0.5)){
  for(lambda1 in c(1,1.25,1.5,1.75,2,3)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # the OOC cdf of SE(theta1,lambda1)
    CL<-theta0-lambda0*log(0.5) # the center line of the chart
    # below we calculate the one-step transition probabilities matrix
    # the ARL calculation is made via the Markov chain method
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
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # this is the ARL
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1)) # print the results
  }
}
######### END ##########################################
