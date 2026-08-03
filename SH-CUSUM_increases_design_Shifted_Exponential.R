alpha1<-0.0015 # lower than 0.002
theta0<-10 # IC theta
lambda0<-100 # IC lambda
UCLsh<-theta0-lambda0*log(alpha1)
for(h1 in seq(5.35,5.35,by=0.01)){
  #h<-4.5
  sims<-25000 # simulation runs
  k<-1.0 # reference value
  mu0X<-theta0+lambda0
  sigma0X<-lambda0
  listrl1<-c() # empty vector to store the RLs
  # CUSUM
  for(j0 in 1:sims){
    c0<-0# starting value
    j<-1 # counter
    z<-rexp(1,rate=(1/lambda0))+theta0 # an obs from SE -- IC
    x<-(z-theta0)/lambda0
    c1<-max(0,c0+x-1-k) # CUSUM
    # use a while and count points until the 1st OOC signal
    while(c1<h1&z<UCLsh){
      j<-j+1
      c0<-c1
      z<-rexp(1,rate=(1/lambda0))+theta0 # an obs from SE -- IC
      x<-(z-theta0)/lambda0
      c1<-max(0,c0+x-1-k) # CUSUM
    }
    listrl1[j0]<-j
  }
  ARL0L<-mean(listrl1)
  print(c(UCLsh,h1,k,ARL0L))
}