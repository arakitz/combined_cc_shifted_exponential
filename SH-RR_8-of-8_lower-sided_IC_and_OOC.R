# below we provide the commands needed to calculate
# the IC and the OOC ARL performance of the SH-RR:8/8 lower one-sided chart
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # IC cdf of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, needed to determine the lower Shewhart limit
LCLsh<-theta0-lambda0*log(1-a1) # the lower Shewhart limit
# below, the use has to try several values of LWL until the IC ARL, given
# the lower Shewhart limit, equals the nominal IC ARL value
# (manual determination of the LWL)
LWL<-75.727
# below we determine the one-step transition probabilities
# we calculate the ARL by using the Markov chain method
p00<-1-F0X(LWL)
p01<-(-F0X(LCLsh)+F0X(LWL))
Q0<-matrix(0,ncol=8,nrow=8)
Q0[,1]<-rep(p00,8)
for(i in 1:7){Q0[i,i+1]<-p01}
e1<-rep(0,8);e1[1]<-1
ID<-diag(8)
l1<-rep(1,8)
ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the IC ARL
print(c(LCLsh,LWL,ARL0)) # print the results to check if the pair (LCLsh,LWL) results in IC ARL close to the nominal value
########## OOC ARL performance ##################
# below we calculate the OOC ARL
# for various OOC values theta1 and lambda1
for(theta1 in c(10,7,5)){
  for(lambda1 in c(100,90,80,70,60,50)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # OOC cdf of the SE(theta1,lambda1)
    # below we determine the one-step transition probabilities
    # we calculate the ARL by using the Markov chain method
    p10<-1-F1X(LWL)
    p11<-(-F1X(LCLsh)+F1X(LWL))
    Q1<-matrix(0,ncol=8,nrow=8)
    Q1[,1]<-rep(p10,8)
    for(i in 1:7){Q1[i,i+1]<-p11}
    e1<-rep(0,8);e1[1]<-1
    ID<-diag(8)
    l1<-rep(1,8)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # the OOC ARL
    print(c(a1,LCLsh,LWL,theta1,lambda1,ARL1)) # print the results
  }
}
################# END #####################################
