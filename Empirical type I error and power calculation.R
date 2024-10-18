# ~~~~~~~~~~~~~~~~~~~~ R program for Marshall Olkin distribution ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~ testing Pareto Type I distribution with censored observations ~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~ Date: 15-09-2024~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set the parameters (replace with your desired values)
lambda1 <-  5
lambda2 <- 1.7
lambda3 <- 1 
# Number of samples to generate
n <- 50  #(50, 100,150)   # sample size

# Number of Monte Carlo simulations
N <- 1000 

# Significance level (alpha) for the test
alpha <- 0.05    # 0.01, 0.02, 0.05

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
test_statistics <- numeric(N)

# Perform Monte Carlo simulations
deltas= rep()
for (i in 1:N) {
  # Generate random variables x and y
  s1=runif(n)
  s2=runif(n)
  x=-1*log(s1)/(theta1+theta12)
  y=-1*log(s2)/(theta2+theta12)
  
  deltas[i]= MOC.Test(x,y)
}
deltas=sort(deltas)
c1= quantile(deltas, 0.025)
c2= quantile(deltas, 0.975)
## power
deltap= rep()
for (i in 1:N) {
  # Generate random variables x and y
  s1=runif(n)
  s2=runif(n)
  x=-1*log(s1)/(theta1+theta12)
  y=-1*log(s2)/(theta2+theta12)
  
  deltap[i]= MOC.Test(x,y)
}

power= mean(as.integer(deltap< c1))+mean(as.integer(deltap> c2))
power




























 