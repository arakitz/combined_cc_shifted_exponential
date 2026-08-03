alpha1<-0.0005 # lower than 0.002
theta0<-10 # IC theta
lambda0<-100 # IC lambda
mu0X<-theta0+lambda0
sigma0X<-lambda0
UCLsh<-theta0-lambda0*log(alpha1)
k<-1.25 # reference value
h1<-4.39 # decision interval
sims<-30000 # simulation runs
for(theta1 in c(10,13,15)){
for(lambda1 in c(100,120,140,160,180,200)){
  listrl1<-c() # empty vector to store the RLs
  # CUSUM
  for(j0 in 1:sims){
    c0<-0# starting value
    j<-1 # counter
    z<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
    x<-(z-theta0)/lambda0
    c1<-max(0,c0+x-1-k) # CUSUM
    # use a while and count points until the 1st OOC signal
    while(c1<h1&z<UCLsh){
      j<-j+1
      c0<-c1
      z<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
      x<-(z-theta0)/lambda0
      c1<-max(0,c0+x-1-k) # CUSUM
    }
    listrl1[j0]<-j
  }
  ARL1<-mean(listrl1)
  cat("SH-UCL:",UCLsh,"k:",k,"h:",h1,"theta1:",theta1,"lambda1:",lambda1,"ARL1:",ARL1,"\n")
}
}
