#' Sequential Bayesian Additive Regression Trees Model
#'
#' A flexible Bayesian nonparametric model that is used as imputation tool for missing covariates.

#' @param x Dataset of covariate matrix with missing values (NAs).
#' @param y Response (fully observed).
#' @param x.type A vector indicating the type of covariates (0=continuous, 1=binary).
#' @param y.type 0=no reponse, 1=continuous response (linear regression used for imputation) and 2=binary response (logistic regression used for imputation), default value = 1
#' @param numimpute Number of Imputed Datasets, default = 5
#' @param numskip Number of iterations skipped, default = 199
#' @param burn Number of iterations for burn-in, default value = 1000
#' @param sigest For continuous variable, an estimate of error variance, sigma^2, used to calibrate an inverse-chi-squared prior used on that parameter. If not supplied, the least-squares estimate is derived instead. See sigquant for more information. Not applicable when variable is binary. Default value is NA.
#' @param seed_dist default value is 12345
#' @param seed_draws defualt value is 99
#' @importFrom stats glm lm binomial na.omit qnorm vcov qchisq
#' @importFrom LaplacesDemon rbern
#' @importFrom msm rtnorm
#' @importFrom Rcpp sourceCpp
#'
#' @useDynLib sbart, .registration = TRUE
#'
#' @return Imputed Dataset Values named as 'ximpute'.
#' @export

seqBART<- function(x, y, x.type, y.type=1, numimpute=5, numskip=199,burn=1000, sigest=NA, seed_dist=12345, seed_draws=99)
{
  set.seed(seed_dist)
  xx<-x
  yy<-y
  pp<-ncol(xx)
  n<-nrow(xx)
  summis<-rep(NA,pp)
  for(i in 1:pp){summis[i]<-sum(is.na(xx[,i]))}
  summis1<-summis
  summis1[x.type==1]<-summis[x.type==1]+n
  summis1[x.type==0]<-summis[x.type==0]+2*n
  summis1[summis==0]<-0

  nvar<-pp-length(summis[summis==0])
  xx_reorder<-xx[,order(summis1)]
  vartype<-x.type[order(summis1)]

##########
  if (y.type==1|y.type==2|y.type==3) {
    xx_reorder<-cbind(xx_reorder,yy)
    vartype<-c(vartype,0)
    nvar<-nvar+1
    pp<-pp+1}

  if (y.type==4) {
    xx_reorder<-cbind(xx_reorder,yy)
    vartype<-c(vartype,1)
    nvar<-nvar+1
    pp<-pp+1}

############
  # if (y.type==1|y.type==2) {
  #   xx_reorder<-cbind(xx_reorder,yy)
  #   vartype<-c(vartype,0)
  #   nvar<-nvar+1
  #   pp<-pp+1}

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
  binum<-length(summis[summis>0&x.type==1])
  if (y.type==4) {binum<-length(summis[summis>0&x.type==1])+1}

  bistart<-0
  z<-c(0,0)

  # if (binum>0) {
  #   bistart<-length(summis[summis==0])
  #   z<-matrix(NA,nrow=n,ncol=binum)
  #   for(i in 1:binum){
  #     bip<-mean(xx_reorder1[,(bistart+i)],na.rm=TRUE)
  #     xx_reorder1[is.na(xx_reorder1[,(bistart+i)]),(bistart+i)]<-rbern(sum(is.na(xx_reorder1[,(bistart+i)])),bip)
  #     z[xx_reorder1[,bistart+i]==0,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==0),mean=qnorm(bip),1,lower=-Inf,upper=0)
  #     z[xx_reorder1[,bistart+i]==1,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==1),mean=qnorm(bip),1,lower=0,upper=Inf)
  #   }
  # }
###########

  if (y.type==4) {
    if (binum>1) {
      bistart<-length(summis[summis==0])
      z<-matrix(NA,nrow=n,ncol=binum)
      for(i in 1:(binum-1)){
        bip<-mean(xx_reorder1[,(bistart+i)],na.rm=TRUE)
        xx_reorder1[is.na(xx_reorder1[,(bistart+i)]),(bistart+i)]<-rbern(sum(is.na(xx_reorder1[,(bistart+i)])),bip)
        z[xx_reorder1[,bistart+i]==0,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==0),mean=qnorm(bip),1,lower=-Inf,upper=0)
        z[xx_reorder1[,bistart+i]==1,i]<-rtnorm(sum(xx_reorder1[,bistart+i]==1),mean=qnorm(bip),1,lower=0,upper=Inf)
      }
      for(i in binum:(binum)){
        bip<-mean(xx_reorder1[,(pp)],na.rm=TRUE)
        xx_reorder1[is.na(xx_reorder1[,(pp)]),(pp)]<-rbern(sum(is.na(xx_reorder1[,(pp)])),bip)
        z[xx_reorder1[,pp]==0,i]<-rtnorm(sum(xx_reorder1[,pp]==0),mean=qnorm(bip),1,lower=-Inf,upper=0)
        z[xx_reorder1[,pp]==1,i]<-rtnorm(sum(xx_reorder1[,pp]==1),mean=qnorm(bip),1,lower=0,upper=Inf)
      }
    } else {
      bistart=pp-1
      z<-matrix(NA,nrow=n,ncol=1)
      bip<-mean(xx_reorder1[,(pp)],na.rm=TRUE)
      z[xx_reorder1[,pp]==0,1]<-rtnorm(sum(xx_reorder1[,pp]==0),mean=qnorm(bip),1,lower=-Inf,upper=0)
      z[xx_reorder1[,pp]==1,1]<-rtnorm(sum(xx_reorder1[,pp]==1),mean=qnorm(bip),1,lower=0,upper=Inf)
    }
  } else {
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
  }




############




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





  if (y.type==0|y.type==3)
    {
    y<-y0-m[pp]
    y[is.na(y)]=0
    }
  else {
    y<-y0
  }

 beta<-NULL
  V<-NULL
  if (y.type==2) {
    data1<-data.frame(y,X)
    beta<-glm(y~.,data=data1, family=binomial)$coefficients
    V <- vcov(glm(y~.,data=data1, family=binomial))
  }

  #numimpute<-5
  numskip<-numskip+1
  nd<-numimpute*numskip

  ##### run BART
  parfit = serBart(x=X,y,burn=burn, nd=nd,nmissing=nvar-1,xmiss=xmiss, sigest=sighat,vartype=vartype,z=z,
                  bistart=bistart,binum=binum,type=y.type,beta=beta,V=V, seed = seed_draws)

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
    if (y.type==1|y.type==2|y.type==3|y.type==4) {
      ximpute<-ximpute[,-(pp)]}
    ximpute<-ximpute[,order(order(summis1))]
    return(ximpute)
  }

  # for (i in 1:numimpute) {
  #   assign(paste("imputed",i,sep=""),impu(i*numskip))
  #
  # }
  #
  #  imputedValues<-list(imputed1 <- imputed1,imputed2 <- imputed2,imputed3 <-imputed3,
  #               imputed4<- imputed4,imputed5<- imputed5)
  #
  #  return(imputedValues)
#############
  retlist<-NULL
  for (i in 1:numimpute) {
    assign(paste("imputed",i,sep=""),impu(i*numskip))
    retlist<-c(retlist,list(impu(i*numskip)))
  }
  names(retlist) <- paste("imputed", 1:numimpute, sep = "")

  return(retlist)
  #############

}

