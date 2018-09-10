# Sequential BART
#
# @param x Dataset of covariate matrix with missing values (NAs).
# @param y Response (fully observed).
# @param nd number of iterations skipped
# @param burn number of iterations for burn-in

# @param m  the number of trees, default = 200
# @param sigdf sig df value, default = 3
# @param sigquant sign quant values, default = 0.90
# @param kfac kd fac value, default = 2.0


# @param sigest asdfsdf
# @param nmissing defualt=0
# @param xmiss default = Null
# @param vartype vartpe
# @param z z
# @param bistart k
# @param binum k
# @param type k
# @param beta k
# @param V k
#
# @return Numeric double values

# @export


serBart=function(x,y,burn=500,nd=1000, nmissing=0,xmiss=NULL,sigest=NA,vartype,z,
                 bistart,binum,type=1,beta=NULL,V=NULL, seed)
{
  cat("***** Running serBart\n")

  m=200
  kfac=2.0


  sigquant=.90
  sigdf=3

  nu = sigdf
  sigq = qchisq(1.0-sigquant,nu)
  lambda = (sigest*sigest*sigq)/nu

  # New Code
  if(type==0|type==3) mifValues <- cpp_bart(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, seed)
  else if(type==4) mifValues <- cpp_bart1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, seed)
  else if(type==1) mifValues <- cpp_bart_y(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, seed)
  else if(type==2) mifValues <- cpp_bart_y1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),beta,t(V),lambda, seed)

  retlist <- list(mif=mifValues)

  return(retlist)
}


