# EARL & EQARL calculation
# for the upper one-sided SH-CUSUM chart
#####################################
alpha1=0.0005 # alpha1, for determining the UCLsh
theta0=10 # IC origin
lambda0=100 # IC scale
k=0 # k+ reference value
h1=22.88 # h+ decision interval
sims=25000 # number of simulation runs
mu0X<-theta0+lambda0 # IC process mean
sigma0X<-lambda0 # IC process standard deviation
UCLsh<-theta0-lambda0*log(alpha1) # Shewhart limit
listARL<-c() # empty vector to store ARLs
listQARL<-c() # empty vector to store quantities needed to calculate EQARL
# define the interval for the shifts in theta
thetamin<-10
thetamax<-13
# define the interval for the shifts in lambda
lambdamin<-100
lambdamax<-130
  for(theta1 in seq(thetamin,thetamax,length=10)){
    for(lambda1 in seq(lambdamin,lambdamax,length=10)){
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
      QARL1<-ARL1*(1/(1+(theta0-theta1)^2))*(1/(1+(lambda0-lambda1)^2))
      listARL<-c(listARL,ARL1)
      listQARL<-c(listQARL,QARL1)
    }
  }
  EARL<-mean(listARL)
  EQARL<-mean(listQARL)
  print(round(c(alpha1,theta0,lambda0,k,h1,EARL,EQARL),digits=5))
#######################
print('######### END ########')