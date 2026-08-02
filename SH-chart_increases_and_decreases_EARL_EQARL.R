# The SH-chart for shifted exponential distribution
# we start with the upper one-sided chart
theta0<-10 # IC theta
lambda0<-100 # IC lambda
# define the interval of shifts in theta
AR1<-10 # thetamin and
AR2<-13 # thetamax
# define the interval of shifts in lambda
BR1<-100 # lambdamin
BR2<-130 #lambdamax
a<-0.002 # nominal False Alarm Rate (FAR)
UCL<-theta0-lambda0*log(a) # Upper control limit
################ EARL calculation via simulation ###############################
sims=5000 # simulation runs
dim1=10 # how many points to take from each interval of theta and lambda OOC values
y1<-seq(AR1,AR2,length=dim1)
y2<-seq(BR1,BR2,length=dim1)
listARLs<-c() # an empty vector to store the ARLs
listQARLs<-c() # an empty vector to store quantities needed to calculate EQARL
for(i1 in 1:dim1){
  for(i2 in 1:dim1){
    listrl1<-c()# empty vector to put the RLs ? OOC.
    thetaOUT<-y1[i1]
    lambdaOUT<-y2[i2]
    for(j0 in 1:sims){
      j<-1 # counting the points until first OOC signal
      x<-rexp(1,rate=(1/lambdaOUT))+thetaOUT # a single obs from SE ? OOC
      while(x<UCL){
        j<-j+1
        x<-rexp(1,rate=(1/lambdaOUT))+thetaOUT
      }
      listrl1[j0]<-j # here the condition on while is violated
    }
    ARL1<-mean(listrl1) # ARL as the sample mean
    listARLs<-c(listARLs,ARL1)
    QARL1<-ARL1*(1/(1+(theta0-thetaOUT)^2))*(1/(1+(lambda0-lambdaOUT)^2))
    listQARLs<-c(listQARLs,QARL1)
  }
}
earl<-mean(listARLs) # EARL
eqarl<-mean(listQARLs) # EQARL
cat('th0:',theta0,'lam0:',lambda0,
    'UCL:',UCL,'thmin:',AR1,'thmax:',AR2,'lammin:',BR1,'lammax:',BR2,'EARL:',earl,'EQARL:',eqarl,'\n')
############# The lower one-sided chart #####################
theta0<-10 # IC theta
lambda0<-100 # IC lambda
# define the interval for shifts in theta
AR1<-7 # thetamin
AR2<-10 # thetamax
# define the interval for shifts in lambda
BR1<-70 # lambdamin
BR2<-100 # lambdamax
######################
a<-0.002 # nominal False Alarm Rate (FAR)
LCL<-theta0-lambda0*log(1-a) # lower control limit
###### EARL calculation ######################
dim1=10 # how many points to consider from each interval
y1<-seq(AR1,AR2,length=dim1)
y2<-seq(BR1,BR2,length=dim1)
listARLs<-c() # empty vector to store the ARLs
listQARLs<-c() # an empty vector to store quantities needed to calculate EQARL
for(i1 in 1:dim1){
  for(i2 in 1:dim1){
    listrl1<-c()# empty vector to put the RLs ? OOC.
    thetaOUT<-y1[i1]
    lambdaOUT<-y2[i2]
    for(j0 in 1:sims){
      j<-1 # counting the points until first OOC signal
      x<-rexp(1,rate=(1/lambdaOUT))+thetaOUT # a single obs from SE ? OOC
      while(x>LCL){
        j<-j+1
        x<-rexp(1,rate=(1/lambdaOUT))+thetaOUT
      }
      listrl1[j0]<-j # here the condition on while is violated
    }
    ARL1<-mean(listrl1) # ARL as the sample mean of the 5000
    listARLs<-c(listARLs,ARL1)
    QARL1<-ARL1*(1/(1+(theta0-thetaOUT)^2))*(1/(1+(lambda0-lambdaOUT)^2))
    listQARLs<-c(listQARLs,QARL1)
  }
}
earl<-mean(listARLs) # EARL
eqarl<-mean(listQARLs) # EQARL
cat('th0:',theta0,'lam0:',lambda0,
    'LCL:',LCL,'thmin:',AR1,'thmax:',AR2,'lammin:',BR1,'lammax:',BR2,'EARL:',earl,'EQARL:',eqarl,'\n')
#######################