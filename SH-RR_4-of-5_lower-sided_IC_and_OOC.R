# below we provide the commands that calculate the IC
# and the OOC ARL performance of the lower one-sided SH-RR:4/5 chart
theta0<-10 # IC origin
lambda0<-100 # IC scale
F0X<-function(x){(1-exp(-(x-theta0)/lambda0))*(x>=theta0)} # IC cdf of the SE(theta0,lambda0)
a1<-0.0005 # alpha1, needed to determine the lower Shewhart limit
LCLsh<-theta0-lambda0*log(1-a1) # the lower Shewhart limit
# now the user has to try different LWL values so as, for the given lower Shewhart limit
# the IC ARL is equal to the nominal value
# (manual determination of the UWL)
LWL<-29.478
CL<-theta0-lambda0*log(0.5) # center line of the chart
# below we construct the IC one-step transition probabilities matrix
# the ARL is calculate via the Markov chain method
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
Q0<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) # the IC transition probabilities matrix
e1<-rep(0,10);e1[1]<-1
ID<-diag(10)
l1<-rep(1,10)
ARL0<-as.vector(e1%*%solve(ID-Q0)%*%l1) # the IC ARL
print(c(LCLsh,LWL,ARL0)) # print the results to check that the pair (LCLsh,LWL) gives IC ARL close to the nominal value
############# OOC ARL performance ###############
# now we can calculate the OOC ARL performance of the chart
# for various OOC values theta1 and lambda1
for(theta1 in c(10,7,5)){
  for(lambda1 in c(100,90,80,70,60,50)){
    F1X<-function(x){(1-exp(-(x-theta1)/lambda1))*(x>=theta1)} # OOC cdf of the SE(theta1,lambda1)
    # below we construct the transition probabilities matrix Q
    # needed to calculate the ARL
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
    Q1<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) # the transition probabilities matrix
    e1<-rep(0,10);e1[1]<-1
    ID<-diag(10)
    l1<-rep(1,10)
    ARL1<-as.vector(e1%*%solve(ID-Q1)%*%l1) # the OOC ARL
    print(c(a1,LCLsh,LWL,theta1,lambda1,ARL1)) # print the results 
  }
}
################### END ##########################
QVALUE <- hcubature(QARLearl,lowerLimit=c(thetamin,lambdamin),
                    upperLimit=c(thetamax,lambdamax))$integral   
EQARL1 <- QVALUE/((thetamax-thetamin)*(lambdamax-lambdamin))
EQARL1 
##########################################
c(EARL1,EQARL1) 
