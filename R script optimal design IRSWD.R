
#############################################################################################################
### Optimal allocation to treatment sequences
### Individually randomized stepped wedge design
### Constant attrition over time
###
### author: Mirjam Moerbeek, Utrecht University, The Netherlands
### lat update: June 1, 2022
#############################################################################################################

### load alabama library
library(alabama)

### Inequality constraint function
### Upper and lower bounds for proportions 
upper=1 # user-specified value of upper boundary (should be <1 and >lower bound)
lower=0 # user-specified value of lower boundary (should be >0 and <upper bound)
hin <- function(x)
{
  nr.seq <- length(x)
  h <- rep(NA, 2*nr.seq)
  h[1:nr.seq] <- upper-x            
  h[(nr.seq+1):(2*nr.seq)] <- x-lower 
  h
}

### Equality constraint function
### Proportions should sum to 1
heq=function(x)
  sum(x)-1

### Function to calculate variance of treatment effect
f.covmat=function(x)
{
  nr.seq=length(x)
  nr.per=nr.seq+1
  exponent <- abs(matrix(1:nr.per - 1, nrow = nr.per, ncol = nr.per, byrow = TRUE) - (1:nr.per - 1))
  VV <- rho^exponent # correlation matrix
  XX=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],rep(1,nr.per))
  
  prop.vec=(1-attrition)^seq(0,nr.seq)-(1-attrition)^seq(1,nr.seq+1)
  prop.vec[nr.seq+1]=1-sum(prop.vec[1:(nr.seq)])
  
  inf.mat=matrix(0,nr.per+1,nr.per+1)
  for(ii in 1:nr.seq)
    {
    XX[ii,(nr.per+1)] <- 0
    for(jj in nr.per:1)
      {
      XXa=XX[1:jj,]
      VVa=VV[1:jj,1:jj]
      if(jj==1)
        {
        XXa=matrix(XXa,nrow=1)
        VVa=matrix(VVa,nrow=1)
      }
      VVainv=solve(VVa)
      inf.mat=inf.mat+x[ii]*prop.vec[jj]*t(XXa)%*%VVainv%*%XXa
     }
  }
  cov.mat=solve(inf.mat)
  cov.mat[(nr.per+1),(nr.per+1)]
}


### User-specified attrition rate between two adjacent time points
attrition=0.0 
### User-specified intraclass correlation
rho=0.4
### User-specified number of sequences
nr.seq=4
### Vector with values to start the optimization procedure
### Note: here the uniform design is used, but any other set of start values can be chosen by the user
start=rep(1/nr.seq,nr.seq)
### Find optimal allocation to sequences
solution=constrOptim.nl(par=start,fn=f.covmat,heq=heq,hin=hin)
### Print optimal allocation to sequences
solution$par

### Calculate efficiency of uniform allocation relative to optimal allocation
solution$value/f.covmat(rep(1/nr.seq,nr.seq))   


