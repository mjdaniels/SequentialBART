
serBart=function(x,y,nd=1000,burn=500,m=200,sigest=NA, sigdf=3, sigquant=.90,kfac=2.0,
                 fname="",nmissing=0,xmiss=NULL,vartype,z,bistart,binum,type=0,beta=NULL,V=NULL)
{
  cat("***** Running serBart\n")

  nu = sigdf
  sigq = qchisq(1.0-sigquant,nu)
  lambda = (sigest*sigest*sigq)/nu

  # New Code
  if(type==0) mifValues <- cpp_bart(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
  else if(type==1) mifValues <- cpp_bart_y(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
  else if(type==2) mifValues <- cpp_bart_y1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),beta,t(V),lambda, type)

  retlist <- list(mif=mifValues)
  return(retlist)
}
