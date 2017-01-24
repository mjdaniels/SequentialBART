#' Seq BART
#'
#' A flexible Bayesian nonparametric model that is used as imputation tool for missing covariates.
#' @param cstem location of the directory
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
#'
#' @importFrom LaplacesDemon rbern
#' @importFrom msm rtnorm

#' @useDynLib bartpkg1
#' @importFrom Rcpp sourceCpp
#'
serBARTfunc3Jan<- function(cstem,xx,yy,datatype,type,
                            numskip=199,burn=1000,
                            m=200,sigdf=3, sigquant=.90,
                            kfac=2.0)
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
  nvaris <<- nvar
  xx_reorder<-xx[,order(summis1)]
  xxreorder_1<<- xx_reorder
  vartype<-datatype[order(summis1)]


  if (type==1|type==2) {
    xx_reorder<-cbind(xx_reorder,yy)
    xxreorder_2<<- xx_reorder
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
  sighatis <<- sighat
  xx_reorder1<-xx_reorder
  xxreorder_3<<- xx_reorder1
  binum<-length(summis[summis>0&datatype==1])
  binumis <<-binum
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

  xxreorder_4<<- xx_reorder1
  X<-matrix(NA,n,pp-1)
  x0=xx_reorder1[,1:pp-1]#in-sample x
  y0=xx_reorder1[,pp]
  ppis <<- pp
  y0is <<- y0
  m<-apply(xx_reorder1,2,mean,na.rm=TRUE)
  m[1:(pp-nvar)]<-0
  m[vartype==1]<-0
  for (j in 1:(pp-1)){
    X[,j]<-x0[,j]-m[j]
  }
  for (j in 1:(pp-1)) {
    X[is.na(X[,j]),j]=0
  }

###################################


  #pppp <<- ncol(X)
  #Xis <<- X
  #write(X,file="/Users/as82986/Desktop/X.txt")
  #transX <- t(X)
  #write(transX,file="/Users/as82986/Desktop/TransX.txt",ncolumns=pppp)
  #newX <<- read.csv("/Users/as82986/Desktop/TransX.txt")

####################################
thisis_y0<<-y0
  if (type==0)
    {
    y<-y0-m[pp]
    y[is.na(y)]=0
    }
  else {
    y<-y0
  }
 thisis_y<<-y


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
  #cat("BURN" ,burn)
  parfit =serBart(cstem,x=X,y,burn=burn,
                  nd=nd,nmissing=nvar-1,xmiss=xmiss,
                  sigest=sighat,vartype=vartype,z=z,
                  bistart=bistart,binum=binum,type=type,beta=beta,V=V)




  #print(class(parfit))
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
  myanswer1<- NULL
  myanswer1<<-list(imputed1 <- imputed1,imputed2 <- imputed2,imputed3 <-imputed3,
                imputed4<- imputed4,imputed5<- imputed5)
  ImputedValuesSet1 <- NULL

  ImputedValuesSet1 <<- imputed1
  #write(imputed1, file = "imputed1file.txt", ncolumns = 4)
  #print("You have reached the end of the function BART")
  # myanswer2 <-list(imputed1 <- imputed1, imputed2 <- imputed2, imputed3 <-imputed3,
  #                  imputed4<- imputed4, imputed5<- imputed5)
  #return(ImputedValuesSet1)

}

serBart=function(cstem,x,y,nd=1000,burn=500,m=200,sigest=NA, sigdf=3, sigquant=.90,kfac=2.0,
                 fname="",nmissing=0,xmiss=NULL,vartype,z,bistart,binum,type=0,beta=NULL,V=NULL)
{
    cat("***** Running serBart\n")
    # print(vartype)

    p=ncol(x)
    #cat( "Value of p is : ", p)
    x_has_values <<- x
    xroot=paste(cstem,fname,"x1.txt",sep="")
    # write(x,file="/Users/as82986/Desktop/x_nocol1.txt")
    # write(t(x),file="/Users/as82986/Desktop/tx_nocol1.txt")
    # write(x,file="/Users/as82986/Desktop/x_col1.txt", ncolumns=4)
    write(t(x),file=xroot,ncolumns=p)

    yroot=paste(cstem,fname,"y1.txt",sep="")
   # print( "y is sent without beinng transposed")
    write(y,file=yroot,ncolumns=1)


    xmissroot=paste(cstem,fname,"xmiss1.txt",sep="")
   # print( "xmiss is sent as transposed")
    write(t(xmiss),file=xmissroot,ncolumns=nmissing+1)

    vartyperoot=paste(cstem,fname,"vartype1.txt",sep="")
  #  print( "vartype is sent as transposed")
    write(t(vartype),file=vartyperoot,ncolumns=p+1)


    #vroot=paste(cstem,fname,"z.txt",sep="") # fix this####################
    #write(t(V),file=vroot,ncolumns=p+1) #

    zroot=paste(cstem,fname,"z1.txt",sep="")
   #  print( "z is sent as transposed")
    if(binum==0){
      write(t(z),file=zroot,ncolumns=1)
    }
    else{
      write(t(z),file=zroot,ncolumns=binum)
      }
    #--------------------------------------------------

     ffname=paste(cstem,fname,"mif.txt",sep="")# fix this####################

    nu=sigdf
    sigq=qchisq(1.0-sigquant,nu)
    sigqis<<-sigq
    sigestis<<- sigest
    lambda = (sigest*sigest*sigq)/nu #lambda parameter for sigma prior #
    Lambda_thiscode_senttocpp <<- lambda
    #typeof(lambda)
    #[1] "double"
    # LAmbda here is a vecotr that has double precision values
    # is.vector = True
    # is.double= True
    # and u can call each value like arrays

    # paste creates a vecotr that has character values  > typeof(mine2)
    #  [1] "character"
    # > is.vector(mine2)
    #  [1] TRUE
    # > is.character(mine2)
    #  [1] TRUE
    # BUT it gives a length = 1
    #> length(mine2)
    # [1] 1
    # soooooo,basially... u get only one entity
    # Run bart
    #--------------------------------------------------
    # Previous Code
    # if (type==0) {
    #   cmd = paste(cstem,'mi_bart',sep='') ##############'mi_bart' is an executable here. You create it at loc/executable and use it later by feeding values to it
                                            ### Instead you can call the cpp func directly here by feedind values to that func

    #   cmd = paste(cmd,xroot,yroot,"0",nd,"0.0",burn,m,nu,kfac,nmissing,xmissroot,bistart,vartyperoot,zroot,ffname)
    # }
    # else if (type==1)
    #{
    #   cmd = paste(cstem,'mi_bart_y',sep='')
    #   cmd = paste(cmd,xroot,yroot,"0",nd,"0.0",burn,m,nu,kfac,nmissing,xmissroot,bistart,vartyperoot,zroot,ffname)
    # }
    # else if (type==2)
    #{
    #   betaroot=paste(cstem,fname,"beta.txt",sep="")
    #   vroot=paste(cstem,fname,"z.txt",sep="")

    #   write(t(V),file=vroot,ncolumns=p+1)
    #   write(beta,file=betaroot,ncolumns=1)
    #
    #   cmd = paste(cstem,'mi_bart_y1',sep='')
    #   cmd = paste(cmd,xroot,yroot,"0",nd,"0.0",burn,m,nu,kfac,nmissing,xmissroot,bistart,vartyperoot,zroot,betaroot,vroot,ffname)
    # }
    # for (i in 1:(nmissing+1)){
    #   cmd<-paste(cmd,lambda[i])
    # }
    # cat(cmd,'\n')
    # print("here")
    # ### THE most important command...
    # system(cmd) ## gonna create the excutable files from cpp that will eventually create mif.txt
    #             ## 'cmd' here is the command that is passed to a shell and it can be anything the shell regards as executable
    # #--------------------------------------------------
    # #read in results
##################################### costants were removed from the arguments in the .cppbarty call####
    # lambda_string<-" "
    # for (i in 1:(nmissing+1)){
    #   lambda_string<-paste(lambda_string,lambda[i])
    # }
    #--------------------------------------------------
    # New Code
    if(type==0)
    {
      ans <- cpp_bart(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),ffname,lambda, type)
     #ans <- cpp_bart(x,y,nd,burn,m,nu,kfac,nmissing,xmiss,bistart,vartype,z,ffname,lambda)
     print( " I never came here") }
    else if(type==1)
    { ans <- cpp_bart_y(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),ffname,lambda, type)
      #ans <- cpp_bart_y(x,y,nd,burn,m,nu,kfac,nmissing,xmiss,bistart,vartype,z,ffname,lambda)
     #print ("I am herer")
     }
    else if(type==2)
      {
      print( " I never came here2")
        #betaroot=paste(cstem,fname,"beta.txt",sep="") # fix this####################
        #write(beta,file=betaroot,ncolumns=1) # fix this####################
        #vroot=paste(cstem,fname,"z.txt",sep="") # fix this####################
        #write(t(V),file=vroot,ncolumns=p+1) # fix this####################
        ans <- cpp_bart_y1(t(x),y,nd,burn,m,nu,kfac,nmissing,t(xmiss),bistart,t(vartype),t(z),beta,t(V),ffname,lambda)
      }
    # for (i in 1:(nmissing+1)){
    #   cmd<-paste(cmd,lambda[i])
    # }
    # ans <- cpp_bart_y(x,y,nd,burn,m,nu,kfac,nmissing,xmiss,bistart,vartype,z,ffname,lambda)
    #print("mif.txt was created by cpp")
    mifthis = scan(paste(cstem,fname,"mif.txt",sep=""))
    retlist <- list(mif=mifthis)
    return(retlist)
    #KK<- ans
    #mm <- list(mif=KK)
    #return(mm)
}
