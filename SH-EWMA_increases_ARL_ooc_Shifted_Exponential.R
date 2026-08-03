# EWMA increases with SH limit -- OOC
theta0<-10 # IC theta
lambda0<-100 # IC lambda
a1<-0.0005 # alpha1, for determing the Shewhart limit
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
lam<-0.30 # reference value
h2<- 283.02 # EWMA limit
sims<-30000 # simulation runs
mu0X<-theta0+lambda0 # IC process mean
sigma0X<-lambda0 # IC process standard deviation
# below we calculate the OOC ARLS
# change theta1, lambda1 according to the IC settings and needs
for(theta1 in c(10,13,15)){
for(lambda1 in c(100,120,140,160,180,200)){
listrl<-c() # empty vector to store the RLs
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
    listrl[j0]<-j
  }
  ARL0<-mean(listrl)
  cat("SH-UCL:",UCLsh,"EWMA-UCL:",h2,"omega:",lam,"theta1:",theta1,"lambda1:",lambda1,"ARL:",ARL0,"\n")
}
}
