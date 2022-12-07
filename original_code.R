############################################################
#
# Codes sent Dec 06 2022 written by Ruth Joy for Gabriel Guaquiere
# 
# 
############################################################
install.packages("mvtnorm");library(mvtnorm)


####################
# offline, compute system noise (sig.e) and obs error variance (sig.o.bar)
# note we set obs error variance as a constant (sig.o.bar)
rm(list=ls())
####################

# setwd("/Users/rutherford/Desktop/SMRU/Jenny Smith")
JS.4121=read.csv("4121_filtervelocity_PostThesis.csv")
idx=362:239315 # index only for non NA values at begin/end of time series

# linear infill for NA's in Calib_Depth
S=approx(idx,JS.4121$Calib_Depth[idx],xout=(idx),method="linear",rule=2)

# create data series vertical "speed"
depth=S$y
speed=diff(depth,lag=1)

####################
# Compute ACF
####################

# set up stuff
ns=length(speed)
ts=1:ns
width=360; #% window with in index units (=5 sec) for 30 minutes. 
overlap=180; #% overlap of adjacent windows. 
nw=floor((length(ts)/(width-overlap) )-1); #% number of windows to use
maxlag=150; #% maximum lag to use for auto-correlation function (ACF)

ACFacc=matrix(NA,maxlag+1,nw); #% accumulator for ACF, 
ACVacc=matrix(NA,maxlag+1,nw); #% accumulator for ACV 
speedacc=matrix(NA,nw,1) # accumulator for mean of absolute value of speed[window]
varacc=matrix(NA,nw,1) # accumulator for variance of speed[]
intcpt=matrix(NA,nw,1) # accumulator for estimate of observation error for each window
ar1acc=matrix(NA,nw,1) # accumulator for ar1 coeff from arima function
## moving window shifts by 400 units each step through k 
## if k=1 then window=(401:1200)-400 = 1:800
k=1 
for (k in k:nw){
	window=((k)*(width-overlap)+1):((k+1)*(width-overlap)+overlap); 
	window=window-overlap #if k=1, then 1:800	
    spd=speed[window]; #% extract the speed in the window
    ACVacc[,k]=acf(spd, lag.max = maxlag,type = c("covariance"),na.action=na.pass, plot=FALSE)$acf; 
    ACFacc[,k]=acf(spd, lag.max = maxlag,type = c("correlation"),na.action=na.pass, plot=FALSE)$acf; 
# insert a few 'ifelse' statements to 
# 1) account for it entire record is 0's
# 2) convergence issues then use index k-1 value
    ar1acc[k]=ifelse(is.numeric(try(arima(spd, order=c(1,0,0))$coef[1],TRUE)),arima(spd, order=c(1,0,0))$coef[1],ifelse(sum(spd)==0,0,ar1acc[k-1]))
	varacc[k]=var(spd)
	speedacc[k]=mean(abs(spd))
	y=ACVacc[(1:5)+1,k]
	y.lm=lm(y~c(1:5)+c(1:5)^2)
	intcpt[k]=y.lm$coef[1]
	if (k%%25==0) print(k) # print every 25th 'iw'
}

##############
# Compute the mean observation error variance
##############
sig.o=varacc-intcpt
sig.o.bar=mean(sig.o)
sig.x=varacc-sig.o.bar

##############
# make model for process variance vs. speed
##############
y=sig.x[sig.x>0]
x=speedacc[sig.x>0] # mean(abs(speed)), but only use those average speeds where process error is > 0
y.lm=lm(y~0+x+I(x^2))
	beta1=y.lm$coefficients[1]
	beta2=y.lm$coefficients[2]
sig.x.hat=beta1*speedacc+beta2*(speedacc^2)

##############
# make model for system noise variance vs. speed
# derived from process variance above
##############
sig.e = sig.x.hat*(1-ar1acc^2)
# set min threshold for sig.e to prevent sig.e from being ~ 0 
sig.e[sig.e<0.035]=0.035  

##############
# Plot to see that sys noise variance scales with speed. 
##############
dev.off()
quartz(height=5,width=10);par(oma=c(3,5,3,5),mar=rep(0,4))
plot(speed,axes=F,col="gray75",type="l");axis(2,las=2)
	mtext(side=2,"speed",line=3)
par(new=T)
plot(sig.e,type="o",axes=F,cex=.6,pch=16,ylab="")
	axis(4,las=2)
	mtext(side=4,"sig.e",line=3)


############################################################
# normal mixture code
normix=function(n,p,sig1,sig2)
{
# distribution 1
# distribution 2
idx=rbinom(n,1,p)
e = c(rnorm(n,mean=0,sd=sig1)[idx==0], rnorm(n,mean=0,sd=sig2)[idx==1])
return(e)
}


############################################################
# mixture random walk codes
rw=function(xold,sigp)
{
#% single stage transition of the (deterministic) model
#% inputs: 
#% xold - state at time t-1
#% theta - parameter vector
#% outputs:
#% xnew - state at time t+1
  
a1=xold[3]; xnewp1=xold[3]+sigp*rnorm(1) # update a1
a2=xold[4]; xnewp2=xold[4]+sigp*rnorm(1) # update a2
# sig.n=exp(xold[5]); xnewp3=xold[5]+sigp*xold[5]*rnorm(1) # update system noise as a function of size at ast time interval
			# could also use sigp to scale the forcing.

#% model: bivariate AR(1)
# update state using "old" values
	D=matrix(c(a1, a2, 1, 0),ncol=2,byrow=T);
	n=matrix(c(normix(n=1,p=0.025,sig1=sig.n,sig2=10*sig.n),0));
xnews=D%*%matrix(xold[1:2],2,1)+n;

xnew=matrix(c(xnews, xnewp1, xnewp2),4,1)
return(xnew)
}


############################################################
# particle filter codes 
pfilter=function(param,psig,y)
{
# accumulators for mean and variance of filtered state
xbar=matrix(NA,nt,nstate); # filter mean

# setup initial particle positions
xas=sig.i*matrix(rnorm(np*ns),np,ns); # original state np x 2
xap1=matrix(param[1],np,1) # parameter 1
xap2=matrix(param[2],np,1) # parameter 2

xa=cbind(xas,xap1,xap2)

#% compute mean and variance for initial conditions
xbar[1,]=colMeans(xa); 


# ==== FILTER/ SMOOTHER IMPLEMENTATION
for (k in 1:nt)  # start time loop
{
	# === prediction step - yields forecast ensemble
xf=t(apply(xa,1,rw,psig)) # predict state forward one time step

	#% === measurement step - yields nowcast ensemble
	# assign weights to each particle (prop to MVN likelihood)
w=apply(matrix(xf[,1],np,1),1,dmvnorm,x=y[k],sigma=R)

	# carry out weighted resample on forecast ensemble
m=sample(1:np,size=np,prob=w,replace=T);

	# resample forecast to yield analysis
xa=xf[m,] 			# here: "a" means analysis, "f" means forecast)
xbar[k,]=apply(xa,2,mean); 	# compute mean

} # end k, for

return(xbar)
}


############################################
# 
# Running my adaptation of Ionides MIF
# You'll have to bring in 3 things for this to work (maybe more?)
# (1) your vertical speed data, 
# (2) sig.o.bar (obs error), 
# (3) sig.e (system error)
#
############################################

numw=floor((length(speed)/(width-overlap) )-1); 
paramMIF=matrix(NA,numw,2) 
paramMIF[1,]=c(0,0)
np=2000; # number of particles, or sample size (put as small number (np=10) until you're sure code is working)
nt=width # no. of timesteps 

ns=2 # number of original state variables
nparam=2 # number of dynamic parameters to be made part of state
nstate=ns+nparam # total size of augmented state

niter=12
iter=1

sig.o=sqrt(sig.o.bar) 	 # observation error std dev
#Id=matrix(0,ncol(y),ncol(y)); diag(Id)=1; R=(sig.o^2)*Id; 
R=matrix(sig.o^2)

sig.n=sqrt(sig.e) 	 # system noise variance

alpha=0.4 # rampdown effect
sigpi=0.1 # initial value for forcing for augmented state parameters
sigp=sapply(1:(niter-1),function(i) (sigpi)*alpha^i)

sig.i=1 				 # variance for the initial conditions
theta=c(sig.i)		 # vector of static parameters (passed to model subroutine)

############################################
# LOOPING THROUGH "iw" windows
############################################
for (iw in 1:numw){
# for (iw in 2:25){

window=((iw-1)*(width-overlap)+1):(iw*(width-overlap)+overlap); 
	y=matrix(speed[window])
	y=sweep(y,1,mean(y));

# ==== RE-SETUP INSIDE LOOP FOR EACH STEP OF "iw"

paramvec=matrix(NA,niter,2) # niter x 2 matrix for a1, a2 for niter iter's
paramvec[1,]=c(1,0)

# MIF unweighted : niter ITERATIONS TO GET TO MLE of STATIC a1, a2 parameters
for (iter in 1:(niter-1))
{
# updates system error, and sigp
X=pfilter(param=paramvec[iter,],psig=sigp[iter],y)
paramvec[iter+1,]=colMeans(X[,3:4])	
# print(iter)
}

paramMIF[iw,]=paramvec[niter,]
print(iw)
} # end iw

