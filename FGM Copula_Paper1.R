############## FGM Copula  #############
Cpla5  #FGM Copula
n <- 100
rNC5 <- r.cpla(Cop=Cpla5,par=c(0.5,30),n)
x<-rNC5[,1]
y<-rNC5[,2]  
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
  n < 100
  
  # Number of Monte Carlo simulations
  n_simulations<- 1000
  
  # Significance level (alpha) for the test
  alpha <- 0.02
  
  # Create empty vector to store test statistic results
  test_statistics <- numeric(n_simulations)
  
  # Perform Monte Carlo simulations
  deltas= rep()
  for (i in 1:n_simulations) {
    # Generate random variables x and y
    theta1 <- 2.0
    theta2 <- 0.5
    theta12 <- 0.4
    gamma1 <- theta12 / (theta1 + theta12)
    gamma2 <- theta12 / (theta2 + theta12)
    
    # Generate marginal random variables (exponential)
    x1 <- rexp(n, rate = theta1)
    y1 <- rexp(n, rate = theta2)
    
    ### Distribution function 
    # Parameters of the exponential distribution
    # Values of x for which you want to calculate the CDF
    x1_values <- seq(0, 10, by = 0.001)
    y1_values <- seq(0, 10, by = 0.001)
    # Calculate the CDF using pexp()
    Fx <- 1 - exp(-theta1 * x1_values)
    Fy <- 1 - exp(-theta2 * y1_values)
    #########Surivaival functions##
    s1=runif(n)
    s2=runif(n)
    
    ########## inversion ##
    x=-1*log(s1)/(theta1+theta12)
    y=-1*log(s2)/(theta2+theta12)
    
    deltas[i]= MOC.Test(x,y)
  }
  deltas=sort(deltas)
  c1= quantile(deltas, 0.005) # lower quantile
  c2= quantile(deltas, 0.995) # upper quantile
  ## power of the test.
  deltap= rep()
  for (i in 1:n_simulations) {
    # Generate random variables x and y
    x<-rNC5[,1]
    y<-rNC5[,2] 
    deltap[i]= MOC.Test(x,y)
  }
  power= mean(as.integer(deltap< c1))+mean(as.integer(deltap> c2))
  
  
  
  
