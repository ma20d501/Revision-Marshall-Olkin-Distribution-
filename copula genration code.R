install.packages("rootSolve")
install.packages("fitdistrplus")
install.packages("EnvStats")
install.packages("readxl")
install.packages("copula")
install.packages("fpc")
install.packages("ggplot2")
install.packages("invgamma")
install.packages("scatterplot3d")
install.packages("numDeriv")
install.packages("graDiEnt")
install.packages("MASS")
install.packages("pracma")


library(graDiEnt)
library(numDeriv)
library(rootSolve)
library(fitdistrplus)
library(EnvStats)
library(readxl)
library(copula)
library(fpc)
library(ggplot2)
library(invgamma)
library(copula)
library(scatterplot3d)
library(MASS)
library(pracma)

# Random sample from bivariate weibull distribution 
alpha1=2
beta1=3
alpha2=3
beta2=4
rho=0.6
parv=c(alpha1, beta1, alpha2, beta2, rho)
#parv
n=50
xv=matrix(0,nrow=n,ncol=2);
for(i in 1:n)
{
  u=runif(1,0,1);
  t=runif(1,0,1);
  #d=1-q1*exp(-alpha1*d1^beta1)
  nc= normalCopula(rho, dim = 2, dispstr = "ex") #for normal copula
  # nc=tCopula(r, dim = 4, dispstr = "ex", df = 4, df.fixed = TRUE, df.min = 0.01) #FOR t Copula
  # nc=archmCopula("clayton",r,dim=4);
  # nc=archmCopula("frank",r, dim=4);
  uv <- rCopula(1, nc);
  x=(-alpha1*log(1-uv[1]))^(1/beta1)
  y=(-alpha2*log(1-uv[2]))^(1/beta2)
  xv[i,]=c(x,y);
  if(x*y==0)
  {
    i=i-1
  }
}
mydata<-data.frame(xv)
x<-mydata[,1]
y<-mydata[,2]
############Test statistics####################
MOC.Test<-function(x,y){ # here X is an array
  n = length(x)
  Z<-rep()
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[i]>x[j]+x[k] & y[i]>y[j]+x[k])
      }
    }
  }
  delta1<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[i]>x[k]+x[j] & y[i]>y[k]+x[j])
      }
    }
  }
  delta2<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[j]>x[i]+x[k] & y[j]>y[i]+x[k])
      }
    }
  }
  delta3<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[j]>x[k]+x[i] & y[j]>y[k]+x[i])
      }
    }
  }
  delta4<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[k]>x[i]+x[j] & y[k]>y[i]+x[j])
      }
    }
  }
  delta5<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[k]>x[j]+x[i] & y[k]>y[j]+x[i])
      }
    }
  }
  delta6<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[i]>x[j]+y[k] & y[i]>y[j]+y[k])
      }
    }
  }
  delta7<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[j]>x[i]+y[k] & y[j]>y[i]+y[k])
      }
    }
  }
  delta8<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[k]>x[i]+y[j] & y[k]>y[i]+y[j])
      }
    }
  }
  delta9<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[i]>x[k]+y[j] & y[i]>y[k]+y[j])
      }
    }
  }
  delta10<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[j]>x[k]+y[i] & y[j]>y[k]+y[i])
      }
    }
  }
  delta11<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-2)){
    for(j in (i+1):n-1){
      for(k in (j+1):n){
        s1<-s1+as.integer(x[k]>x[j]+y[i] & y[k]>y[j]+y[i])
      }
    }
  }
  delta12<-s1/choose(n,4)
  s1<-0
  for(i in 1:(n-3)){
    for(j in (i+1):n-2){
      for(k in (j+1):n-1){
        for(l in (k+1):n){
          s1<-s1+as.integer(prod(x[l]>x[k] & y[l]>y[k] & x[j]>y[i] & y[j]>y[i] ) )
        }
      }
    }
  }
  delta001<-s1/choose(n,4)
  delta=((4*delta1+4*delta2+4*delta3+4*delta4+4*delta5+4*delta6+4*delta7+4*delta8+4*delta9+4*delta10+4*delta11+4*delta12
          -12*delta001)/24)
  return(delta)
}

# Sample size for Monte Carlo simulations
n < 50
# Number of Monte Carlo simulations
n_simulations<- 100
# Significance level (alpha) for the test
alpha <- 0.01
# Create empty vector to store test statistic results
test_statistics <- numeric(n_simulations)
# Perform Monte Carlo simulations
deltas= rep()
for (i in 1:n_simulations) {
  # Generate random variables x and y
  theta1 <- 2.0
  theta2 <- 0.5
  theta12 <- 0.4
  s1=runif(n)
  s2=runif(n)
  x=-1*log(s1)/(theta1+theta12)
  y=-1*log(s2)/(theta2+theta12)
  deltas[i]= MOC.Test(x,y)
}
deltas=sort(deltas)
c1= quantile(deltas, 0.025) # lower quantile
c2= quantile(deltas, 0.975) # upper quantile
## power of the test.
deltap= rep()
for (i in 1:n_simulations) {
  # Generate random variables x and y
  # Random sample from bivariate weibull distribution 
  alpha1=2
  beta1=3
  alpha2=3
  beta2=4
  rho=0.6
  parv=c(alpha1, beta1, alpha2, beta2, rho)
  #parv
  n=50
  xv=matrix(0,nrow=n,ncol=2);
  for(i in 1:n)
  {
    u=runif(1,0,1);
    t=runif(1,0,1);
    #d=1-q1*exp(-alpha1*d1^beta1)
    nc= normalCopula(rho, dim = 2, dispstr = "ex") #for normal copula
    # nc=tCopula(r, dim = 4, dispstr = "ex", df = 4, df.fixed = TRUE, df.min = 0.01) #FOR t Copula
    # nc=archmCopula("clayton",r,dim=4);
    # nc=archmCopula("frank",r, dim=4);
    uv <- rCopula(1, nc);
    x=(-alpha1*log(1-uv[1]))^(1/beta1)
    y=(-alpha2*log(1-uv[2]))^(1/beta2)
    xv[i,]=c(x,y);
    if(x*y==0)
    {
      i=i-1
    }
  }
  mydata<-data.frame(xv)
  x<-mydata[,1]
  y<-mydata[,2]
  deltap[i]= MOC.Test(x,y)
}
deltap
power= mean((deltap[i]< c1))+mean((deltap[i]> c2))

power= mean(as.integer(deltap[i]< c1))+mean(as.integer(deltap[i]> c2))

require(mvtnorm)
S <- matrix(c(1,.8,.8,1),2,2) #Correlation matrix
AB <- rmvnorm(mean=c(0,0),sig=S,n=1000) #Our gaussian variables
U <- pnorm(AB) #Now U is uniform - check using hist(U[,1]) or hist(U[,2])
x <- qgamma(U[,1],2) #x is gamma distributed
y <- qbeta(U[,2],1,2) #y is beta distributed
plot(x,y) #They correlate!



# Load the copula package
library(copula)
# Set the parameters for the t-copula
df <- 3  # Degrees of freedom
rho <- 0.5  # Correlation parameter
# Create a t-copula object
tc <- tCopula(param = rho, df = df)
# Generate random samples from the t-copula
n_samples <- 1000  # Number of samples
samples <- rCopula(n_samples, tc)
x<- samples[,1]
y<- samples[,2]
y



#' Random Number Generation for Copula Functions.
#' @description Random Number Generation for Bivariate Copula Functions. Returns a number of pairs of random values 
#' from a bivariate copula with marginal distributions X and Y. 
#' @param typeCopula Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"}, \code{"AMH"},
#' \code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.
#' @param theta A numeric value for the space parameter.
#' @param typeX Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and
#' \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1_X A numeric value for the first parameter of the first marginal distribution.
#' @param num2_X A numeric value for the second parameter of the first marginal distribution.
#' Only required for two parameter distributions.
#' @param typeY Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and \code{"Gamma"}.
#' Defaults to \code{"Exp"}.
#' @param num1_Y A numeric value for the first parameter of the second marginal distribution.
#' @param num2_Y A numeric value for the second parameter of the second marginal distribution. 
#' Only required for two parameter distributions. 
#' @param nsim Number of observations to be generated.
#' 
#' @return
#' 2-dimensional random vector with the results of the simulation.
#' @examples 
#' 
#' res<-rcopula(typeCopula = 'clayton', theta = 2, typeX='Exp', num1_X=0.9, 
#'              typeY='Exp', num1_Y=0.3, nsim=1000)
#' 
#' res
#' 
#' res2<-rcopula(typeCopula = 'AMH', theta = 2, typeX='Norm', num1_X=0.9, num2_X=0.3, 
#'               typeY='Gamma', num1_Y=3, num2_Y=2, nsim=1000)
#'               
#' res2[,2]
#' 
#' @author Gustavo Soutinho, Luis Meira-Machado

rcopula<- function(typeCopula = 'clayton', theta = 1, typeX='Exp', num1_X=1, num2_X=NULL, 
                   typeY='Exp', num1_Y=1, num2_Y=NULL, nsim=500){
  
  TAB<-NULL
  
  
  for(i in 1:nsim){
    
    #i<-1
    
    v1<-runif(1,0,1)
    
    v2<-runif(1,0,1)
    
    res<-copula(v1, v2, theta=theta,type=typeCopula,typeX=typeX, num1_X=num1_X, num2_X=num2_X,
                typeY=typeY, num1_Y=num1_Y, num2_Y=num2_Y)
    
    x<-res[1]
    
    y<-res[2]
    
    TAB<-rbind(TAB,cbind(i, x, y))
    
  }
  
  TAB<-as.data.frame(TAB)
  
  colnames(TAB)<-c('ID','X','Y')
  
  #res<-list(tab=TAB, typeCopula=typeCopula,teta=teta,typeX=typeX, num1_X=num1_X, num2_X=num2_X, 
  #          typeY=typeY, num1_Y=num1_Y, num2_Y=num2_Y,
  #          nsim=nsim)
  
  
  res<-TAB[,2:3]
  
  return(res)
  
}
res<-rcopula(typeCopula = 'clayton', theta = 2, typeX='Exp', num1_X=0.9, 
             typeY='Exp', num1_Y=0.3, nsim=1000)

res

res2<-rcopula(typeCopula = 'AMH', theta = 2, typeX='Norm', num1_X=0.9, num2_X=0.3, 
              typeY='Gamma', num1_Y=3, num2_Y=2, nsim=1000)

res2[,2]
