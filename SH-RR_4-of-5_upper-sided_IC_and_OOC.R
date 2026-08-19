# Calculation of the IC and the OOC ARL performance
# for the SH-RR:4/5 upper one-sided chart
theta0<-0 # IC origin
lambda0<-1 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # IC cdf of the SE(theta0,lambda0) distribution
a1<-0.0005 # alpha1, needed to determine the upper Shewhart limit UCLsh
UCLsh<-theta0-lambda0*log(a1) # the upper Shewhart limit
# below we provide a manual way to determine the upper warning limit
# this is the last stage of the procedure
for(UWL in seq(183.16,183.17,by=0.001)){
  CL<-theta0-lambda0*log(0.5) # center line of the chart
  # below we calculate the transition probabilities needed to construct the one-step transition probabilities matrix
  # we use the Markov chain method to calculate the ARL of the chart
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
  Q0<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) # the IC transition probabilities matrix Q
  e1<-rep(0,10);e1[1]<-1
  ID<-diag(10)
  l1<-rep(1,10)
  ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the IC ARL
  print(c(UCLsh,UWL,ARL0)) # print the IC ARL along with the Shewhart limit and the UWL
}
############# OOC ARL performance ###############
UWL<-1.732 # user has to provide the value obtained in the previous step
CL<-theta0-lambda0*log(0.5) # center line of the chart
# below we calculate the OOC ARL for various OOC values theta1 and lambda1
for(theta1 in c(0,0.25,0.5)){
  for(lambda1 in c(1,1.25,1.5,1.75,2,3)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # the OOC cdf of the SE(theta1,lambda1)
    # below we construct the OOC transition probabilities matrix Q
    # we use the Markov chain method to calculate the OOC ARL of the chart
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
    Q1<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) # the transition probability matrix
    e1<-rep(0,10);e1[1]<-1
    ID<-diag(10)
    l1<-rep(1,10)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # the OOC ARL
    print(c(a1,UCLsh,UWL,theta1,lambda1,ARL1)) # print the results
  }
}
############# END ################################
