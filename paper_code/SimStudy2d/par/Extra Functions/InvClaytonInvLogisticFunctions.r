condsim.ievl<-function(x,a,U)
{
  dummy<-function(y)
  {
    1-x^(1/a-1)*(x^(1/a)+y^(1/a))^(a-1) *exp(x-(x^(1/a)+y^(1/a))^(a)) - U
  }
  ur<-uniroot(dummy,interval = c(0,20))
  return(ur$root)
}

condsim.ievl<-Vectorize(condsim.ievl,vectorize.args = c("x","U"))


PX1.givenx2<-function(x1,x2,beta)
{
  1-exp((beta+1)*x2)*(exp(beta*x1)+exp(beta*x2)-1)^(-1/beta-1)
}

PX3.givenx2<-function(x3,x2,alpha)
{
  a<-alpha
  1-x2^(1/a-1)*(x2^(1/a)+x3^(1/a))^(a-1) *exp(x2-(x2^(1/a)+x3^(1/a))^(a))
}

Joint.cdf<-function(x123,alpha,beta)
{
  x123[1]->x1
  x123[2]->x2
  x123[3]->x3
  
  to.int<-function(v)
  {
    PX1.givenx2(x1=x1,x2=v,beta=beta)*PX3.givenx2(x3=x3,x2=v,alpha=alpha)*exp(-v)
  }
  int<-integrate(to.int,0,x2)
  return(int$value)
}