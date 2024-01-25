syspss_r <- function(x,n){ 
  N <- length(x)
  U <- sample(N,N)
  xx <- x[U]
  z <- n*cumsum(xx)/sum(xx)
  r <- runif(1)
  s <- numeric(n)
  k <- 1
  for(i in 1:N){
    if(z[i]>=r){
      s[k] <- U[i]
      r <- r+1
      k <- k+1
    }
  }
  return(s)
}

syspss_r_pij<-function(x,s,nsims=10){
  N<-length(x)
  n<-length(s)
  p<-matrix(0,n,n)
  for(k in 1:nsims){
    ss<-syspss_cpp(x,n)
    for(i in 1:(n-1)){
      for(j in (i+1):n){ 
        if(min(abs(ss-s[i]))+min(abs(ss-s[j]))==0) 
          p[i,j]<-p[i,j]+1 
        }
    }
    }
  #p<-(p+t(p))/nsims
  return(p)
}


# 
# set.seed(123)
# x <- rnorm(2000)
# z <- x + max(x)
# d <- syspss_cpp(z, 10)
# mm <- syspss_pij_cpp(z, d, 100000)

