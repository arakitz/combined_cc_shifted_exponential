# EWMA decreases with SH limit -- OOC -- EARL and EQARL calculation
theta0=10 # IC origin
lambda0=100 # IC scale
a1=0.0005 # alpha1, for determing the LCLsh
lam=0.10 # smoothing parameter
h2=63.01566 # LCLe 
sims=5000 # simulation runs
LCLsh<-theta0-lambda0*log(1-a1) # Shewhart lower control limit
mu0X<-theta0+lambda0 # IC process mean
sigma0X<-lambda0 # IC process standard deviation
listARL<-c() # empty vector to store the ARLs
listQARL<-c() # empty vector to store quantities needed to calculate EQARL
# define the interval for the shifts in origin
thetamin<-7
thetamax<-10
# define the interval for the shifts in scale
lambdamin<-70
lambdamax<-100
################################################
for(theta1 in seq(thetamin,thetamax,length=10)){
for(lambda1 in seq(lambdamin,lambdamax,length=10)){
listrl1<-c() # empty vector to store the RLs
  # EWMA
  for(j0 in 1:sims){
    c0<-mu0X # starting value
    j<-1 # counter
    x<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
#print(x)
    c1<-lam*x+(1-lam)*c0 # EWMA
    # use a while and count points until the 1st OOC signal
    while(c1>h2&x>LCLsh){
      j<-j+1
      c0<-c1
      x<-rexp(1,rate=(1/lambda1))+theta1
#print(x)
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
print(round(c(a1,LCLsh,lam,h2,theta0,lambda0,EARL,EQARL),digits=5))
print('######### END #############')