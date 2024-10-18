## MArshal olkin copula ##
library(copula)
MC<- 1000 # no. of repition
power<- rep()
alpha <- 0.10 # 0.02,0.2,0.10 (level of significance)
sam= c(25,50,75,100)
for(d in 1:4){
  n<-sam[d]
    # Function to generate Marshall-Olkin copula samples
    moCopula <- function(n, dim) {
      # Generate uniform samples for each dimension
      u <- matrix(runif(n * dim), nrow = n, ncol = dim)
      # Generate Marshall-Olkin copula samples
      y <- matrix(0, nrow = n, ncol = dim)
      for (i in 1:n) {
        y[i, 1] <- u[i, 1]
        for (j in 2:dim) {
          y[i, j] <- min(u[i, j], prod(1 - y[i, 1:(j - 1)]))
        }
      }
      return(y)
    }# Example: Generate Marshall-Olkin copula samples with dimensionality 2 and 100 samples
    samples <- moCopula(n, 2)
    x<- samples[,1]
    y<- samples[,2]
    ## test statistics 
    MOC.Test<-function(x,y){ # here X is an array
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
    }
    ###### 
    deltas= rep() 
    for (r in 1:MC) {
      # Generate random variables x and y
      samples <- moCopula(n, 2)
      x<- samples[,1]
      y<- samples[,2]
      deltas[r]= MOC.Test(x,y)
    }
    deltas=sort(deltas)
    c1= quantile(deltas, 0.025)
    c2= quantile(deltas, 0.975)
    ## power
    deltap= rep()
    for (r in 1:MC) {
      # Generate random variables x and y
      samples <- moCopula(n, 2)
      x<- samples[,1]
      y<- samples[,2]
      deltap[r]= MOC.Test(x,y)
    }
    power[d]= mean(as.integer(deltap[deltap!='NaN']< c1))+mean(as.integer(deltap[deltap!='NaN']> c2))
  }
power

