# EWMA increases with SH limit, determine the values of the design parameters
theta0<-10 # IC theta
lambda0<-100 # IC lambda
a1<-0.0005 # alpha1, for determing Shewhart limit
UCLsh<-theta0-lambda0*log(a1) # Shewhart limit
lam<-0.30 # reference value
mu0X<-theta0+lambda0 # IC process mean
sigma0X<-lambda0 # IC process standard deviation
# manual search for the UCLe
# below is the final stage of searching
for(h2 in seq(282.95,283.05,by=0.01)){
  sims<-25000 # simulation runs
  listrl<-c() # empty vector to store the RLs
  # EWMA
  for(j0 in 1:sims){
    c0<-mu0X # starting value
    j<-1 # counter
    x<-rexp(1,rate=(1/lambda0))+theta0 # an obs from SE -- IC
    # z<-(x-theta0)/lambda0
    c1<-lam*x+(1-lam)*c0 # EWMA
    # use a while and count points until the 1st OOC signal
    while(c1<h2&x<UCLsh){
      j<-j+1
      c0<-c1
      x<-rexp(1,rate=(1/lambda0))+theta0
      # z<-(x-theta0)/lambda0
      c1<-lam*x+(1-lam)*c0 # EWMA
    }
    listrl[j0]<-j
  }
  ARL0<-mean(listrl)
  #  ARL0
  #  SDRL0<-sd(listrl)
  #  SDRL0
  #  SDRL0/sqrt(sims)
  print(c(UCLsh,h2,lam,ARL0))
}
