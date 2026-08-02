# EWMA increases with SH limit -- OOC -- EARL calculation
# below we provide the code that is used
# for the calculation of the EARL, EQARL 
#########################################################################
sims<-5000
theta0<-10 # IC origin
lambda0<-100 # IC scale
a1<-0.0005 # alpha1 for the UCLsh
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
h2<-153.085 # control limit of the EWMA
lam<-0.05 # smoothing parameter
mu0X<-theta0+lambda0 # IC mean of the process
sigma0X<-lambda0 # IC standard deviation of the process
listARL<-c() # an empty vector to store the ARLs
listQARL<-c() # an empty vector to store the quantities for calculating EQARL
# define the interval for the shifts in origin
thetamin<-10
thetamax<-13
# define the interval for the shifts in scale
lambdamin<-100
lambdamax<-130
for(theta1 in seq(thetamin,thetamax,length=10)){
for(lambda1 in seq(lambdamin,lambdamax,length=10)){
listrl1<-c() # empty vector to store the RLs
  # EWMA
  for(j0 in 1:sims){
    c0<-mu0X # starting value
    j<-1 # counter
    x<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
    c1<-lam*x+(1-lam)*c0 # EWMA
    # use a while and count points until the 1st OOC signal
    while(c1<h2&x<UCLsh){
      j<-j+1
      c0<-c1
      x<-rexp(1,rate=(1/lambda1))+theta1
      c1<-lam*x+(1-lam)*c0 # EWMA
    }
    listrl1[j0]<-j
  }
ARL1<-mean(listrl1)
QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
listARL<-c(listARL,ARL1)
listQARL<-c(listQARL,QARL1)
}
}
EARL<-mean(listARL)
EQARL<-mean(listQARL)
print(round(c(a1,UCLsh,lam,h2,theta0,lambda0,EARL,EQARL),digits=5))
######################################
print('######### END #############')