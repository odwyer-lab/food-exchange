# as requested by journal staff, I am combining into a single file the four scripts used for Figs 2 and 3 in the main ms, and SI 1 and SI 2 in the Supplementary information. Each script will be separated by % ======= and titled appropriately

#============

#CODE TO GENERATE FIG 2 plots

library(deSolve)
#set number of locations
m<-10
#set parameter values
# what is equilibrium standing resource pool (equal in our model in all locations)
r<-1
# consumption parameter c (defining the per capita, per unit resource consumption rate of an individual)
c<-1
# equilibrium res pool and c combine to determine the per capita mortality rate in absence of food
mu<-r*c
#build matrix P via whatever rules you wish, but here I will pick the raw entries of P from the absolute values of a normal distribution. I'll use these in multiple P matrices, each time rescaling these raw values to show the effect of large/smaller overall exchange but with the same network structure and relative values
rawPentries<-abs(rnorm(m^2,mean=0,sd=1))
#set the number of environmental phases considered. At the beginning of each phase, rho_i will be redrawn in each location, representing the effect of a fluctuating environment.
nphase=25
# define the interval that each phase will last for.  I chose this to be set by the timescale defined by mortality (so that this interval was longer than this timescale). If shorter, then the population will bnot be able to `catch up'  with the environmental changes.  That alternative regime is certainly of interest too but not what I am considering here.
interval<-10/mu
# determine three (or more) different "scales" to rescale the strength of exchange by.
exchrange<-c(0,0.3,1)
# draw the mean value of each rho in each of the N different locations.  Parameter values (and in fact this choice of distribution) arbitrary, just chosen for convenience.
meanrho<-rnbinom(m,size=10,mu=100)
#define the ODEs we'll be using with reference to the main text. Note that this requires package "deSolve".
CR<-function(time,abundances,parameters){
  R<-abundances[1:m]
  N<-abundances[(m+1):(2*m)]
  P<-matrix(parameters$P,nrow=m,ncol=m)
  mu<-parameters$mu
  rho<-parameters$rho
  dRdt<-rho-R*N+P%*%N-N*colSums(P)
  dNdt<-N*R-mu*N
  return(list(c(dRdt,dNdt)))
}
# six- plot grid layout to show three different outcomes for both resources and consumers
layout(matrix(c(1:6),nrow=3,ncol=2,byrow=TRUE))
for(exch in exchrange){
  #compute P matrix by rescaling raw entries appropriately in each case
  P<-matrix(exch*rawPentries,nrow=m,ncol=m)
  diag(P)<-0
  # compute Q and M matrices from P, as described in the text of the paper
  Q<-(1/mu)*(P-diag(colSums(P)))
  M<-diag(1,m)-Q
  #solve for the equilibrium population sizes, N_i, which I'll put in a vector called n. 
  n<-(1/mu)*solve(M)%*%meanrho
  # set initial conditions for R to have some random perturbations.  Or could set them equal to the equilibrium values if desired.
  R<-matrix(rep(r,m)+0.1*runif(m),nrow=1,ncol=m,byrow=T)
  # set initial conditions for N equal to equilibrium values.
  N<-matrix(n,nrow=1,ncol=m,byrow=T)
  times<-0
  for(phase in 1:nphase){
    #now cycling through the dynamics of each phase.  For each phase we'll solve the approriate ODEs using the final conditions from the previous phase as the new initial conditions. Most of the focus here is by design on the new equilibrium values that population sizes will reach
    # at beginning of each phase, calculate rho in each location by drawing a value from a negative binomial with mean = meanrho, and dispersion parameter "size" chosen to be proportional to mean rho, such that variance = 11*meanrho here, thus providing an example distribution showing Taylor's temporal power law with z=1, and beta= 11.
    rho<-sapply(meanrho,function(x){
      rnbinom(n=1,size=x/10,mu=x) 
    })
    #set time range to be the interval length, and adust the timestep as desired
    time_range<-seq(0,interval,by=interval/100)
    #set new initial abundances
    initial_abundances<-c(R[dim(R)[1],],N[dim(N)[1],])
    # solve ODES from these initial conditions using previously-defined function CR
    output_LV<-ode(initial_abundances,time_range,func=CR,parms=list(P=c(P),rho=c(rho),mu=c(mu)))
    #extract the relevant quantities---time points (only the new ones), resource values, population values.
    times<-c(times[-length(times)],times[length(times)]+output_LV[,1])
    R<-rbind(R[-dim(R)[1],],output_LV[,2:(m+1)])
    N<-rbind(N[-dim(N)[1],],output_LV[,(m+2):(2*m+1)])
  }
  
  matplot(times,N[,1:m],pch=1,lty=1,type="l",main="Population Sizes through time",xlab="Time",ylab="Abundances")
  matplot(times,R[,1:m],pch=1,lty=1,type="l",ylim=c(0.5*r,1.5*r),main="Resource Pools through time",xlab="Time",ylab="Abundances")
}
#==========

#CODE TO GENERATE FIG 3 Plots

library(deSolve)
#set number of locations
m<-10
#set parameter values
# what is equilibrium standing resource pool (equal in our model in all locations)
r<-1
# consumption parameter c (defining the per capita, per unit resource consumption rate of an individual)
c<-1
# equilibrium res pool and c combine to determine the per capita mortality rate in absence of food
mu<-r*c
#define the ODEs we'll be using with reference to the main text. Note that this requires package "deSolve".
CR<-function(time,abundances,parameters){
  R<-abundances[1:m]
  N<-abundances[(m+1):(2*m)]
  P<-matrix(parameters$P,nrow=m,ncol=m)
  mu<-parameters$mu
  rho<-parameters$rho
  dRdt<-rho-R*N+P%*%N-N*colSums(P)
  dNdt<-N*R-mu*N
  return(list(c(dRdt,dNdt)))
}
#set the number of phases _within_ a time interval with a given set of rho_i values. Each of these phases will experience a perturbation to resource pools and/or population, and then allow the population return to equilibrium (or not, if it is unstable). Note that this is distinct from the uncorrelated fluctuations to rho_i used in Figure 2 of the manuscript.
nphase=50
# interval will be the time period before the next perturbation to population sizes/resource pools.  Designed to allow the populations to approximately return to equilibrium (if it is stable).
interval<-10*mu
# In this script we will investigate either a ring-like exchange process, which we expect to be unstable below a critical resource inflow rate, and a symmetric version of the same ring, whocih we expect to be stable. Choose as desired.
sym<-FALSE
# define the matrix P
P<-4*cbind(rbind(rep(0,m-1),diag(m-1)),c(1,rep(0/m,m-1)))
if(sym){P<-0.5*(P+t(P))}
#define the matrices Q and M as per main text
Q<-(1/mu)*(P-diag(colSums(P)))
M<-diag(1,m)-Q
#choose the population sizes. This is somewhat arbitrary, but bear in mind that the transition from instability to instability will only appear if you rescale resource inflow over the right range. I.e. your range of resource inflow rates has to cover either side of this transition, and whether it does or not will in part be defined by how big the resource flows are to begin with, combined with how much they are rescaled by at each phase, combined with how many phases you look at.
n<-10
# ninit will be the initial values of population abundances. This is not dynamically important but I use it later for plot limits
ninit<-n
# compute the rho values corresponding to this set of n. Note that as discussed in the main paper, check whether this set of n is possible with a given P matrix and positive values of rho_i. I.e. make sure that your rho computed below has positive values, if not then there is no set of positive rho_i which (when combined with some alternative choice of P to my ring structure above) will lead to your choice of n.
rho<-mu*M%*%rep(n,m)
R<-matrix(rep(r,m)*(0.95+0.1*runif(m)),nrow=1,ncol=m,byrow=T)
N<-matrix(n,nrow=1,ncol=m,byrow=T)
times<-0
for(phase in 1:nphase){
  time_range<-seq(0,interval,by=interval/100)
  # at beginning of each interval perturb both R and N (or choose what you wish)
  initial_abundances<-c(R[dim(R)[1],]*(0.95+0.1*runif(m)),N[dim(N)[1],]*(0.95+0.1*runif(m)))
  
  output_LV<-ode(initial_abundances,time_range,func=CR,parms=list(P=c(P),rho=c(rho),mu=c(mu)))
  
  #extract the relevant quantities (time points and res/pop values)
  times<-c(times[-length(times)],times[length(times)]+output_LV[,1])
  R<-rbind(R[-dim(R)[1],],output_LV[,2:(m+1)])
  N<-rbind(N[-dim(N)[1],],output_LV[,(m+2):(2*m+1)])
  #add perturbation every phase, but only every 10 phases reduce sypply system wide by a chosen factor (I choise 0.7 foe the figure, but this is arbitrary)
  if(phase%%10==0){
    #compute rescaled rho_i vector
    rho<-0.7*rho
    #compute new equilibrium N_i
    n<-(1/mu)*solve(M)%*%rho
  }
}
layout(matrix(c(1:2),nrow=1,ncol=2,byrow=TRUE))
matplot(times,N[,1:m],pch=1,lty=1,type="l",ylim=c(0.1*ninit,1.1*ninit),main="Population Sizes through time",xlab="Time",ylab="Abundances",log="y")
matplot(times,R[,1:m],pch=1,lty=1,type="l",ylim=c(0.5*r,1.5*r),main="Resource Pools through time",xlab="Time",ylab="Abundances")

#==============

#CODE TO GENERATE SI FIG 1 Plots


library(deSolve)
#set number of locations
m<-10
#set parameter values
# what is equilibrium standing resource pool (equal in our model in all locations)
r<-1
# consumption parameter c (defining the per capita, per unit resource consumption rate of an individual)
c<-1
# equilibrium res pool and c combine to determine the per capita mortality rate in absence of food
mu<-r*c
#choose the exponent y, which will determine the dependence of trade rates on the standing resource pool.
y<-0.5
#define the ODEs we'll be using with reference to the SI text. Note that this requires package "deSolve".
CR<-function(time,abundances,parameters){
  R<-abundances[1:m]
  N<-abundances[(m+1):(2*m)]
  P<-matrix(parameters$P,nrow=m,ncol=m)
  mu<-parameters$mu
  rho<-parameters$rho
  dRdt<-rho-R*N+P%*%((c*R/mu)^y*N)-((c*R/mu)^y*N)*colSums(P)
  dNdt<-N*R-mu*N
  return(list(c(dRdt,dNdt)))
}

#set the number of phases _within_ a time interval with a given set of rho_i values. Each of these phases will experience a perturbation to resource pools and/or population, and then allow the population return to equilibrium (or not, if it is unstable). Note that this is distinct from the uncorrelated fluctuations to rho_i used in Figure 2 of the manuscript.
nphase=50
# interval will be the time period before the next perturbation to population sizes/resource pools.  Designed to allow the populations to approximately return to equilibrium (if it is stable).
interval<-10*mu
# In this script we will investigate either a ring-like exchange process, which we expect to be unstable below a critical resource inflow rate, and a symmetric version of the same ring, whocih we expect to be stable. Choose as desired.
sym<-FALSE
# define the matrix P
P<-4*cbind(rbind(rep(0,m-1),diag(m-1)),c(1,rep(0/m,m-1)))
if(sym){P<-0.5*(P+t(P))}
#define the matrices Q and M as per main text
Q<-(1/mu)*(P-diag(colSums(P)))
M<-diag(1,m)-Q
#choose the population sizes. This is somewhat arbitrary, but bear in mind that the transition from instability to instability will only appear if you rescale resource inflow over the right range. I.e. your range of resource inflow rates has to cover either side of this transition, and whether it does or not will in part be defined by how big the resource flows are to begin with, combined with how much they are rescaled by at each phase, combined with how many phases you look at.
#
# note that with the new more general dependence of exchange on R, the location of the transition from stability to instability seems to be systematically lower with larger y.  I kept all other params the same as in Fig 3, so this means that in general I had to start with smaller population sizes as y got bigger. You could rescale parameters to infalte pop sizes to more realistic values if desired.
n<-2.5
# ninit will be the initial values of population abundances. This is not dynamically important but I use it later for plot limits
ninit<-n
# compute the rho values corresponding to this set of n. Note that as discussed in the main paper, check whether this set of n is possible with a given P matrix and positive values of rho_i. I.e. make sure that your rho computed below has positive values, if not then there is no set of positive rho_i which (when combined with some alternative choice of P to my ring structure above) will lead to your choice of n.
rho<-mu*M%*%rep(n,m)
R<-matrix(rep(r,m)*(0.95+0.1*runif(m)),nrow=1,ncol=m,byrow=T)
N<-matrix(n,nrow=1,ncol=m,byrow=T)
times<-0
for(phase in 1:nphase){
  time_range<-seq(0,interval,by=interval/100)
  # at beginning of each interval perturb both R and N (or choose what you wish)
  initial_abundances<-c(R[dim(R)[1],]*(0.95+0.1*runif(m)),N[dim(N)[1],]*(0.95+0.1*runif(m)))
  
  output_LV<-ode(initial_abundances,time_range,func=CR,parms=list(P=c(P),rho=c(rho),mu=c(mu)))
  
  #extract the relevant quantities (time points and res/pop values)
  times<-c(times[-length(times)],times[length(times)]+output_LV[,1])
  R<-rbind(R[-dim(R)[1],],output_LV[,2:(m+1)])
  N<-rbind(N[-dim(N)[1],],output_LV[,(m+2):(2*m+1)])
  #add perturbation every phase, but only every 10 phases reduce sypply system wide by a chosen factor (I choise 0.7 foe the figure, but this is arbitrary)
  if(phase%%10==0){
    #compute rescaled rho_i vector
    rho<-0.7*rho
    #compute new equilibrium N_i
    n<-(1/mu)*solve(M)%*%rho
  }
}
layout(matrix(c(1:2),nrow=1,ncol=2,byrow=TRUE))
matplot(times,N[,1:m],pch=1,lty=1,type="l",ylim=c(0.1*ninit,1.1*ninit),main=paste("Population Sizes through time\n","y = ",y,sep=""),xlab="Time",ylab="Abundances",log="y")
matplot(times,R[,1:m],pch=1,lty=1,type="l",ylim=c(0.5*r,1.5*r),main=paste("Resource Pools through time\n","y = ",y,sep=""),xlab="Time",ylab="Abundances")

#========

#CODE TO GENERATE SI Fig 2 plots


library(deSolve)
#set number of locations
m<-10
#set parameter values
# what is equilibrium standing resource pool (equal in our model in all locations)
r<-1
# consumption parameter c (defining the per capita, per unit resource consumption rate of an individual)
c<-1
# equilibrium res pool and c combine to determine the per capita mortality rate in absence of food
mu<-r*c
#choose the exponent y, which will determine the dependence of trade rates on the standing resource pool.
y1<-0.25
#choose the exponent x, which will determine the dependence of trade rates on the current pop size.
y2<-0.75
#define the ODEs we'll be using with reference to the SI text. Note that this requires package "deSolve".
CR<-function(time,abundances,parameters){
  R<-abundances[1:m]
  N<-abundances[(m+1):(2*m)]
  P<-matrix(parameters$P,nrow=m,ncol=m)
  mu<-parameters$mu
  rho<-parameters$rho
  dRdt<-rho-R*N+P%*%(((c*R/mu)^y1)*N^y2)-(((c*R/mu)^y1)*N^y2)*colSums(P)
  dNdt<-N*R-mu*N
  return(list(c(dRdt,dNdt)))
}

#set the number of phases _within_ a time interval with a given set of rho_i values. Each of these phases will experience a perturbation to resource pools and/or population, and then allow the population return to equilibrium (or not, if it is unstable). Note that this is distinct from the uncorrelated fluctuations to rho_i used in Figure 2 of the manuscript.
nphase=50
# interval will be the time period before the next perturbation to population sizes/resource pools.  Designed to allow the populations to approximately return to equilibrium (if it is stable).
interval<-10*mu
# In this script we will investigate either a ring-like exchange process, which we expect to be unstable below a critical resource inflow rate, and a symmetric version of the same ring, whocih we expect to be stable. Choose as desired.
layout(matrix(c(1:4),nrow=2,ncol=2,byrow=TRUE))
for(sym in c(FALSE,TRUE)){
  if(sym){type<-"Symmetrized Kula"
  }else{
    type<-"Kula"}
# define the matrix P
P<-4*cbind(rbind(rep(0,m-1),diag(m-1)),c(1,rep(0/m,m-1)))
if(sym){P<-0.5*(P+t(P))}
#define the matrices Q and M as per main text
Q<-(1/mu)*(P-diag(colSums(P)))
M<-diag(1,m)-Q
#choose the population sizes. This is somewhat arbitrary, but bear in mind that the transition from instability to instability will only appear if you rescale resource inflow over the right range. I.e. your range of resource inflow rates has to cover either side of this transition, and whether it does or not will in part be defined by how big the resource flows are to begin with, combined with how much they are rescaled by at each phase, combined with how many phases you look at.
#
# note that with the new more general dependence of exchange on R, the location of the transition from stability to instability seems to be systematically lower with larger y.  I kept all other params the same as in Fig 3, so this means that in general I had to start with smaller population sizes as y got bigger. You could rescale parameters to infalte pop sizes to more realistic values if desired.
n<-3
# ninit will be the initial values of population abundances. This is not dynamically important but I use it later for plot limits
ninit<-n
# compute an initial set of rho values 
rho<-mu*M%*%rep(n,m)
  R<-matrix(rep(r,m)*(0.95+0.1*runif(m)),nrow=1,ncol=m,byrow=T)
  N<-matrix(n,nrow=1,ncol=m,byrow=T)
  times<-0
for(phase in 1:nphase){
time_range<-seq(0,interval,by=interval/100)
# at beginning of each interval perturb both R and N (or choose what you wish)
initial_abundances<-c(R[dim(R)[1],]*(0.95+0.1*runif(m)),N[dim(N)[1],]*(0.95+0.1*runif(m)))

output_LV<-ode(initial_abundances,time_range,func=CR,parms=list(P=c(P),rho=c(rho),mu=c(mu)))

#extract the relevant quantities (time points and res/pop values)
times<-c(times[-length(times)],times[length(times)]+output_LV[,1])
R<-rbind(R[-dim(R)[1],],output_LV[,2:(m+1)])
N<-rbind(N[-dim(N)[1],],output_LV[,(m+2):(2*m+1)])
#add perturbation every phase, but only every 10 phases reduce sypply system wide by a chosen factor (I choise 0.7 foe the figure, but this is arbitrary)
if(phase%%10==0){
  #compute a rescaled rho_i vector
  rho<-0.7*rho
#compute new equilibrium N_i
#n<-(1/mu)*solve(M)%*%rho
}
}
matplot(times,N[,1:m],pch=1,lty=1,type="l",ylim=c(0.01*ninit,1.1*ninit),main=paste("Population Sizes through time\n",type,", y1 = ",y1,", y2 = ",y2,sep=""),xlab="Time",ylab="Abundances",log="y")
matplot(times,R[,1:m],pch=1,lty=1,type="l",ylim=c(0.5*r,1.5*r),main=paste("Resource Pools through time\n",type,", y1 = ",y1,", y2 = ",y2,sep=""),xlab="Time",ylab="Abundances")
}