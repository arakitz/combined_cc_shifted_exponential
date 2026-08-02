# EARL & EQARL calculation
# for the SH-CUSUM lower-sided (for detecting decreases)
alpha1=0.0005 # alpha1 for determining Shewhart limit
theta0=10 # IC origin
lambda0=100 # IC scale
k=0.0 # k-, reference value
h1=20.934 # h-, decision interval
sims=5000 # number of simulation runs
mu0X<-theta0+lambda0 # IC process mean
sigma0X<-lambda0 # IC process standard deviation
LCLsh<-(alpha1>=0.0005)*theta0-lambda0*log(1-alpha1) # Shewhart limit, lower limit 
listARL<-c() # empty vector to store the ARLs
listQARL<-c() # empty vector to store quantities needed for the calculation of EQARL
# define the interval of shifts in theta
thetamin<-7
thetamax<-10
# define the interval of shifts in lambda
lambdamin<-70
lambdamax<-100
  for(theta1 in seq(thetamin,thetamax,length=10)){
    for(lambda1 in seq(lambdamin,lambdamax,length=10)){
      listrl1<-c() # empty vector to store the RLs
      ############ CUSUM
      for(j0 in 1:sims){
        c0<-0 # starting value
        j<-1 # counter
        z<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
        x<-(z-theta0)/lambda0
        c1<-max(0,c0-x+1-k) # CUSUM
        ############ use a while and count points until the 1st OOC signal
        while(c1<=h1&z>LCLsh){
          j<-j+1
          c0<-c1
          z<-rexp(1,rate=(1/lambda1))+theta1 # an obs from SE -- IC
          x<-(z-theta0)/lambda0
          c1<-max(0,c0-x+1-k) # CUSUM
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
  print(round(c(alpha1,LCLsh,theta0,lambda0,k,h1,EARL,EQARL),digits=5))

#########################
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.0,h1=20.934,sims=3000)
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.15,h1=8.212,sims=3000)
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.25,h1=5.245,sims=3000)
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.50,h1=1.905,sims=3000)
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.75,h1=0.544,sims=3000)
#fun1(alpha1=10^(-12),theta0=10,lambda0=100,k=0.85,h1=0.2512,sims=3000)
#########################
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.0,h1=35.01,sims=3000)
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.15,h1=11.32,sims=3000)
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.25,h1=7.02,sims=3000)
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.50,h1=2.44,sims=3000)
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.75,h1=0.69,sims=3000)
fun1(alpha1=0.0015,theta0=10,lambda0=100,k=0.85,h1=0.307,sims=3000)
###############################
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.0,h1=27.65,sims=10^5)
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.15,h1=9.79,sims=10^5)
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.25,h1=6.12,sims=10^5)
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.50,h1=2.16,sims=10^5)
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.75,h1=0.616,sims=10^5)
fun1(alpha1=0.001,theta0=10,lambda0=100,k=0.85,h1=0.2753,sims=10^5)
###############################
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.0,h1=23.61,sims=10^5)
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.15,h1=8.85,sims=10^5)
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.25,h1=5.58,sims=10^5)
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.50,h1=2.011,sims=10^5)
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.75,h1=0.575,sims=10^5)
fun1(alpha1=0.0005,theta0=10,lambda0=100,k=0.85,h1=0.262,sims=10^5)
################## END ####################################