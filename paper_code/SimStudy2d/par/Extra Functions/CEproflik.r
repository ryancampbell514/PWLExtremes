#=================================
# For conditional extremes fit

HT.nll.PL<-function(theta,dat,u,which)
{
  theta[1]->a
  theta[2]->b
  
  # Assume exponential margins and positive dependence
   if(a< 0||a>1||b>=1){return(10e10)}
  # let x always be the large variable
  x<-dat[,which]
  y<-dat[,c(1,2)[-which]]
  
  exc<-x>u
  x<-x[exc]
  y<-y[exc]
  
  mu<-mean((y-a*x)/x^b)
  sig<-sqrt(mean(((y-a*x-mu*x^b)/x^b)^2))
  
  nll<--sum(dnorm(y,mean=a*x+mu*x^b, sd=sig*x^b,log=T))
  return(nll)
}
