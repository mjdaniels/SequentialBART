# Sequential BART
#
# @param x Dataset of covariate matrix with missing values (NAs).
# @param y Response (fully observed).
# @param nd number of iterations skipped
# @param burn number of iterations for burn-in
# @param m m value
# @param sigest asdfsdf
# @param sigdf sig df value
# @param sigquant sign quant values
# @param kfac kd fac value
# @param fname fanme value
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

serBart=function(x,y,nd=1000,burn=500,m=200,sigest=NA, sigdf=3, sigquant=.90,kfac=2.0,
                 fname="",nmissing=0,xmiss=NULL,vartype,z,bistart,binum,type=0,beta=NULL,V=NULL)
{
  cat("***** Running serBart\n")

  # print("####################")
  # print(Sys.getenv())
  # print("####################")

  nu = sigdf
  sigq = qchisq(1.0-sigquant,nu)
  lambda = (sigest*sigest*sigq)/nu

  # New Code
  if(type==0) mifValues <- cpp_bart(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
  else if(type==1) mifValues <- cpp_bart_y(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
  else if(type==2) mifValues <- cpp_bart_y1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),beta,t(V),lambda, type)

  #k<<- NULL
  retlist <- list(mif=mifValues)
  #k <<- list(mif=mifValues)
  return(retlist)
}


# #' @importFrom stats glm lm binomial na.omit qnorm vcov qchisq
# #' @importFrom LaplacesDemon rbern
# #' @importFrom msm rtnorm
# #' @importFrom Rcpp sourceCpp
#
# #' @useDynLib bartpkg1
