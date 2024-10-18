#Required Package
library(mvtnorm)
# Number of bivariate samples
n <- 25 # 25,50,75,100 #sample size
# Significance level (alpha) for the test
alpha <- 0.1  # 0.02,0.2,0.10 (level of significance)
# Number of Monte Carlo simulations
n_simulations<- 1000

# Define the parameters for the bivariate normal distribution
meanlog <- c(0, 0)     # Mean vector on the log scale for both variables
sdlog <- c(1, 1)       # Standard deviation vector on the log scale for both variables
# Create a covariance matrix (for demonstration, it's diagonal for independence)
cov_matrix <- diag(2)
# Generate bivariate random samples from the normal distribution
samples_normal <- rmvnorm(n, mean = meanlog, sigma = cov_matrix)
# Transform the normal samples to lognormal
samples_lognormal <- exp(samples_normal)
# Extract x and y values
x <- samples_lognormal[, 1]
y <- samples_lognormal[, 2]
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

# Create empty vector to store test statistic results
test_statistics <- numeric(n_simulations)

# Perform Monte Carlo simulations
deltas= rep()
for (i in 1:n_simulations) {
  # Generate random variables x and y
  lambda1 <- 2.0
  lambda2 <- 0.5
  lambda12 <- 0.4
  # Generate independent exponential random variables
  x1 <- rexp(n, rate = lambda1)
  x2 <- rexp(n, rate = lambda2)
  x3 <- rexp(n, rate = lambda3)
  # Generate the bivariate sample
  x <- pmin(x1, x3)
  y <- pmin(x2, x3)
  # Calculate sums
  sum_x <- sum(x)
  sum_y <- sum(y)
  # Calculate n3 (number of times x_i == y_i)
  n3 <- sum(x == y)
  # Method of Moments Estimators
  theta1 <- (n / sum_x) - (n3 / sum_y) / (1 + n3 / n)
  theta2 <- (n / sum_y) - (n3 / sum_x) / (1 + n3 / n)
  theta12 <- n3 * ((1 / sum_x) + (1 / sum_y)) / (1 + n3 / n)
  # Output the estimates
  theta1
  theta2
  theta12
  
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
  x <- samples_lognormal[, 1]
  y <- samples_lognormal[, 2]
  
  deltap[i]= MOC.Test(x,y)
}
power= mean(as.integer(deltap< c1))+mean(as.integer(deltap> c2))




