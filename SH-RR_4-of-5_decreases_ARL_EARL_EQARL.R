# EARL and EQARL calculation for the SH-RR:4/5 lower one-sided chart
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # CDF of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, for determining the lower Shewhart limit
LCLsh<-theta0-lambda0*log(1-a1) # The lower Shewhart control limit
LWL<-29.478 # the lower warning limit
CL<-theta0-lambda0*log(0.5) # the center line of the chart (median of the SE(theta0,lambda0) distribution)
# create the transition probabilities matrix
p00<-1-F0X(CL)
p01<-(-F0X(LWL)+F0X(CL))
p02<-(-F0X(LCLsh)+F0X(LWL))
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
ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the IC ARL (calculated by the Markov chain method)
############# OOC performance in terms of EARL, EQARL ###############
listARL<-c() # empty vector to store the ARLs
listQARL<-c() # empty vector to store quantities needed to calculate the EQARL
# define the interval for the shifts in origin
thetamin<-7
thetamax<-10
# define the interval for the shifts in scale
lambdamin<-70
lambdamax<-100
for(theta1 in seq(thetamin,thetamax,length=10)){
  for(lambda1 in seq(lambdamin,lambdamax,length=10)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # CDF of the SE(theta1,lambda1)
    # create the transition probabilities matrix
    p10<-1-F1X(CL)
    p11<-(-F1X(LWL)+F1X(CL))
    p12<-(-F1X(LCLsh)+F1X(LWL))
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
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # OOC ARL
    QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
    listARL<-c(listARL,ARL1)
    listQARL<-c(listQARL,QARL1)
  }
}
EARL<-mean(listARL)
EQARL<-mean(listQARL)
print(round(c(a1,LCLsh,LWL,theta0,lambda0,EARL,EQARL),digits=5))
print('######### END #############')