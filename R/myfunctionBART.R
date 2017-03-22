#' Sequential Bayesian Additive Regression Trees Model
#'
#' A flexible Bayesian nonparametric model that is used as imputation tool for missing covariates.

#' @param xx Dataset of covariate matrix with missing values (NAs).
#' @param yy Response (fully observed).
#' @param datatype a vector indicating the type of covariates (0=continuous, 1=binary).
#' @param type 0=no reponse, 1=continuous response (linear regression used for imputation) and 2=binary response (logistic regression used for imputation)
#' @param numskip number of iterations skipped
#' @param burn number of iterations for burn-in
#' @param m m value
#' @param sigdf sig df value
#' @param sigquant sign quant values
#' @param kfac kd fac value
#'
#' @return Imputed Dataset Values
#' @export
#' @importFrom stats glm lm binomial na.omit qnorm vcov qchisq


#' @importFrom LaplacesDemon rbern
#' @importFrom msm rtnorm

#' @useDynLib sbart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
seqBART<- function( xx, yy, datatype, type=1, numskip=199,burn=1000, m=200,sigdf=3, sigquant=.90, kfac=2.0)
{
  set.seed(12345)

  pp<-ncol(xx)
  n<-nrow(xx)
  summis<-rep(NA,pp)
  for(i in 1:pp){summis[i]<-sum(is.na(xx[,i]))}
  summis1<-summis
  summis1[datatype==1]<-summis[datatype==1]+n
  summis1[datatype==0]<-summis[datatype==0]+2*n
  summis1[summis==0]<-0

  nvar<-pp-length(summis[summis==0])
  xx_reorder<-xx[,order(summis1)]
  vartype<-datatype[order(summis1)]


  if (type==1|type==2) {
    xx_reorder<-cbind(xx_reorder,yy)
    vartype<-c(vartype,0)
    nvar<-nvar+1
    pp<-pp+1}

  ## missing data indicator
  xmiss=matrix(rep(0,(nvar)*n),ncol=(nvar))
  for (j in 1:(nvar)){
    xmiss[is.na(xx_reorder[,pp+1-j]),j]=1
  }


  sighat<-rep(NA,nvar)
  for(j in (pp):(pp-nvar+1)){
    data1<-data.frame(xx_reorder[,1:(j-1)])
    lmf = lm(xx_reorder[,j]~.,na.action=na.omit,data=data1)
    sighat[pp+1-j] = summary(lmf)$sigma
  }

  xx_reorder1<-xx_reorder
  binum<-length(summis[summis>0&datatype==1])
  bistart<-0
  z<-c(0,0)
  if (binum>0) {
    bistart<-length(summis[summis==0])
    z<-matrix(NA,nrow=n,ncol=binum)
    for(i in 1:binum){
      bip<-mean(xx_reorder1[,(bistart+i)],na.rm=TRUE)
      xx_reorder1[is.na(xx_reorder1[,(bistart+i)]),(bistart+i)]<-rbern(sum(is.na(xx_reorder1[,(bistart+i)])),bip)
      z[xx_reorder1[,bistart+i]==0,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==0),mean=qnorm(bip),1,lower=-Inf,upper=0)
      z[xx_reorder1[,bistart+i]==1,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==1),mean=qnorm(bip),1,lower=0,upper=Inf)
    }
  }

  X<-matrix(NA,n,pp-1)
  x0=xx_reorder1[,1:pp-1]#in-sample x
  y0=xx_reorder1[,pp]

  m<-apply(xx_reorder1,2,mean,na.rm=TRUE)
  m[1:(pp-nvar)]<-0
  m[vartype==1]<-0
  for (j in 1:(pp-1)){
    X[,j]<-x0[,j]-m[j]
  }
  for (j in 1:(pp-1)) {
    X[is.na(X[,j]),j]=0
  }

  if (type==0)
    {
    y<-y0-m[pp]
    y[is.na(y)]=0
    }
  else {
    y<-y0
  }

 beta<-NULL
  V<-NULL
  if (type==2) {
    data1<-data.frame(y,X)
    beta<-glm(y~.,data=data1, family=binomial)$coefficients
    V <- vcov(glm(y~.,data=data1, family=binomial))
  }

  numimpute<-5
  numskip<-numskip+1
  nd<-numimpute*numskip

  ##### run BART
  parfit = serBart(x=X,y,burn=burn, nd=nd,nmissing=nvar-1,xmiss=xmiss, sigest=sighat,vartype=vartype,z=z,
                  bistart=bistart,binum=binum,type=type,beta=beta,V=V)

  miff<-parfit$mif
  mif_matrix<-t(matrix(miff,ncol=(burn+nd)))
  mif_m1<-mif_matrix[(burn+1):(burn+nd),]
  xtrans<-y
  for(j in 1:(nvar-1)){
    xtrans<-cbind(xtrans,X[,pp-j])
  }
  xtrans<-c(t(xtrans))
  miss_ind<-c(t(xmiss))

  #####imputed dataset
  impu<-function(j){
    xtrans[miss_ind==1]<-mif_m1[j,]
    xtrans1<-t(matrix(xtrans,ncol=n))
    xtrans2<-xtrans1[,(nvar):1]
    ximpute<-xx_reorder[,1:(pp-nvar)]
    for(j in 2:(nvar+1)){
      ximpute<-cbind(ximpute,xtrans2[,j-1]+m[pp-nvar-1+j])}
    if (type==1|type==2) {
      ximpute<-ximpute[,-(pp)]}
    ximpute<-ximpute[,order(order(summis1))]
    return(ximpute)
  }

  for (i in 1:numimpute) {
    assign(paste("imputed",i,sep=""),impu(i*numskip))
  }
   #imputedValues<-NULL
   imputedValues<-list(imputed1 <- imputed1,imputed2 <- imputed2,imputed3 <-imputed3,
                imputed4<- imputed4,imputed5<- imputed5)

   return(imputedValues)
}

# serBart=function(cstem,x,y,nd=1000,burn=500,m=200,sigest=NA, sigdf=3, sigquant=.90,kfac=2.0,
#                  fname="",nmissing=0,xmiss=NULL,vartype,z,bistart,binum,type=0,beta=NULL,V=NULL)
# {
#     cat("***** Running serBart\n\n")
#
#     nu = sigdf
#     sigq = qchisq(1.0-sigquant,nu)
#     lambda = (sigest*sigest*sigq)/nu
#
#     # New Code
#     if(type==0) mifValues <- cpp_bart(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
#     else if(type==1) mifValues <- cpp_bart_y(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),lambda, type)
#     else if(type==2) mifValues <- cpp_bart_y1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),beta,t(V),lambda, type)
#
#     retlist <- list(mif=mifValues)
#     return(retlist)
# }
