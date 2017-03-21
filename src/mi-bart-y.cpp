#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Cholesky"

#include "rng.h"
#include "info.h"


#include "funs.h" // for mcmc + initinput AND for the fit in the main function
#include "func_mcmc.h"
#include "func_initialInput.h"

using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using namespace Eigen;



//--------------------------------
// main program
//--------------------------------
// [[Rcpp::export]]
NumericVector cpp_bart_y (NumericVector new_xroot, NumericVector new_yroot, int new_nd, int new_burn,
                          int new_m, int new_nu, int new_kfac, int  new_nmissing,
                          IntegerVector new_xmissroot, int new_bistart, NumericVector new_vartyperoot,
                          NumericVector new_zroot, NumericVector new_lambda, int new_type)
{
  size_t m;
  double kfac;
  std::vector<double> x,y;
  std::vector<double> z;
  size_t n,p;
  double nu;
  double lambda;
  int bistart = 0;
  int binum = 0;

  NumericVector mymif2;

  uint seed=99; //random number generation
  RNG gen(seed); //this one random number generator is used in all draws

  //read in data
  y = Rcpp::as< std::vector<double> >(new_yroot);
  n = y.size();
  //n = 0;
  if(n<1)
  {
    //cout << "error n<1\n";
    return mymif2;
  }

  //read x
  x = Rcpp::as< std::vector<double> >(new_xroot); // ***aparently x is not a 2d vector.
  p = x.size()/n;
  if(x.size() != n*p)
  {
    //cout << "error: input x file has wrong number of values\n";
    return mymif2;
  }


  m = new_m;
  nu = new_nu;
  kfac = new_kfac;
  size_t burn  = new_burn;
  size_t nd = new_nd;
  size_t nvar = new_nmissing;

  std::vector<int> ind_missing = Rcpp::as< std::vector<int> >(new_xmissroot);

  bistart = new_bistart;

  int* vartype=new int[p+1];

  std::vector<int> vart2 = Rcpp::as< std::vector<int> >(new_vartyperoot);

  int jj=0;
  for (std::vector<int>::const_iterator i = vart2.begin(); i != vart2.end(); ++i)
  {
    vartype[jj]=*i;
    if(jj>=bistart)
    {
      if(*i==1)
      {
        binum++;
      }
    }
    jj++;
  }

  z = Rcpp::as< std::vector<double> >(new_zroot);

  init *all_initial = new init[nvar+1]; //right

  for(size_t i=0;i<nvar+1;i++)
  {
    all_initial[i].ftemp=new double [n];
    all_initial[i].r=new double [n];
    //cout<<all_initial[i].r<<endl;
  }

  for(size_t i=0;i<nvar+1;i++){
    init_input(i, all_initial[i], vartype,  m, kfac, bistart, binum, n, p, z, y, x);
  }

  // mcmc
  double pro_prop=0.0; //proposal probability of MH min(1,alpha)
  double pro_propnum=1.0;//numerator of proposal prob
  double pro_propdenom=1.0;//denominator of proposal prob
  double prop=0.0;// proposal for mh
  dinfo di_temp; //x used to get fitted value
  double* ppredmean=new double[p];
  double* fpredtemp=new double[1];

  // new change1  4/21/2015
  MatrixXd xx(n,p+1);
  VectorXd x1=VectorXd::Constant(n,1);
  VectorXd beta(p+1);
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > yy(y.data(),n,1);
  VectorXd ee(n);
  double* sigsq=new double[1];
  Map<Eigen::Matrix<double,1,1> > sigsq1(sigsq);
  MatrixXd cho(p+1,p+1);
  MatrixXd xxinv(n,n);
  std::vector<double> rnor;
  for(size_t i=0;i<=p;i++) rnor.push_back(0.0);
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > rbeta(rnor.data(),p+1,1);
  double* meantemp1=new double[1];
  double* meantemp2=new double[1];
  Map<Eigen::Matrix<double,1,1> > meantemp11(meantemp1);
  Map<Eigen::Matrix<double,1,1> > meantemp21(meantemp2);
  VectorXd xxtemp(p+1);

  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > xx1(all_initial[0].di.x,n,p);
  xx<<x1,xx1;
  beta=(xx.transpose()*xx).inverse()*xx.transpose()*yy;
  di_temp.n=1;
  di_temp.y=0;
  std::vector<double> tempx;
  for(size_t i=0;i<=p;i++) tempx.push_back(0.0);
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > xxtemp1(tempx.data(),p,1);
  xxtemp<<1,xxtemp1;
  size_t l=0;

  for(size_t i=0;i<(nd+burn);i++){
    // mcmc used to get updated parameters
    for(size_t j=1;j<nvar+1;j++){
      lambda= (new_lambda[j]);
      mcmc(all_initial[j], gen,vartype,  m,  n,  lambda, nu);
    }

    //new change2 -- posterior of beta and sigma
    xx<<x1,xx1;
    ee=yy-xx*beta;

    sigsq1=(ee.transpose()*ee);
    sigsq[0]=sigsq[0]/gen.chi_square(double(n-p-1));
    xxinv=(xx.transpose()*xx).inverse()*sigsq[0];
    cho=(xxinv.llt().matrixL());
    gen.normal(rnor,0,1);
    beta=xxinv/sigsq[0]*xx.transpose()*yy+cho*rbeta;
    jj=0;
    //-------------------------------
    // impute missing values
    for(size_t j=0;j<n;j++){
      //the row of x used to get prediction
      for(size_t k=1;k<=nvar;k++) //change3
      {
        if(ind_missing[j*(nvar+1)+k]==1)
        {
          for(size_t j1=0;j1<p;j1++){
            tempx[j1]=all_initial[0].di.x[j*p+j1];
          }
          if(vartype[p-k]==0){
            prop=gen.normal(all_initial[k].allfit[j],all_initial[k].pi.sigma);
            pro_propdenom=1.0;
            pro_propnum=1.0;
          }
          else{
            prop=double(abs(all_initial[k].z[j]-1));
            pro_propnum=fabs(prop-phi(-all_initial[k].allfit[j]));
            pro_propdenom=fabs(double(all_initial[k].z[j])-phi(-all_initial[k].allfit[j]));
          }
          // change 4
          xxtemp<<1,xxtemp1;
          meantemp11=(xxtemp.transpose()*beta);
          pro_propdenom=pro_propdenom*pn(all_initial[0].y[j],meantemp1[0],sigsq[0]);
          tempx[p-k]=prop;//change the element to proposal
          xxtemp<<1,xxtemp1;
          meantemp21=(xxtemp.transpose()*beta);
          pro_propnum=pro_propnum*pn(all_initial[0].y[j],meantemp2[0],sigsq[0]);

          l=1; //change 5
          while(l<k){
            if(vartype[p-l]==0){
              pro_propdenom=pro_propdenom*pn(all_initial[l].y[j],all_initial[l].allfit[j],pow(all_initial[l].pi.sigma,2));
            }
            else{
              pro_propdenom=pro_propdenom*fabs(double(all_initial[l].z[j])-phi(-all_initial[l].allfit[j]));
            }

            di_temp.p=p-l;
            di_temp.x=&tempx[0];
            di_temp.vartype=vartype;
            ppredmean[l]=0.0;
            for(size_t q=0;q<m;q++)
            {
              fit(all_initial[l].t[q],all_initial[l].xi,di_temp,fpredtemp);
              ppredmean[l] += fpredtemp[0];
            }
            if(vartype[p-l]==0){
              pro_propnum=pro_propnum*pn(all_initial[l].y[j],ppredmean[l],pow(all_initial[l].pi.sigma,2));
            }
            else{
              pro_propnum=pro_propnum*fabs(double(all_initial[l].z[j])-phi(-ppredmean[l]));
            }
            l++;
          }
          if(vartype[p-k]==0){
            pro_prop=std::min(1.0,pro_propnum/pro_propdenom);
          }
          else{
            pro_prop=pro_propnum/(pro_propnum+pro_propdenom);
          }
          //if proposal accepted, update variables
          if(gen.uniform()<pro_prop){
            if(vartype[p-k]==0){
              all_initial[k].y[j]=prop;
            }
            else{
              all_initial[k].z[j]=int(prop);
              if(all_initial[k].z[j]==0){
                all_initial[k].y[j]=0.0;
                while(all_initial[k].y[j]>=0.0){
                  all_initial[k].y[j]=gen.normal(all_initial[k].allfit[j],1);
                }
              }
              else{
                all_initial[k].y[j]=0.0;
                while(all_initial[k].y[j]<=0.0){
                  all_initial[k].y[j]=gen.normal(all_initial[k].allfit[j],1);
                }
              }
            }
            l=0;
            while(l<k){
              all_initial[l].di.x[j*(p-l)-k+p]=prop;
              all_initial[l].allfit[j]=ppredmean[l];
              l++;
            }
            mymif2.push_back(prop);
          }
          else{
            if(vartype[p-k]==0){
              mymif2.push_back(all_initial[k].y[j]);
            }
            else{
              mymif2.push_back(all_initial[k].z[j]);
            }
          }
        }
      }
    }
  }
  // cout << "28Feb"<<endl;
  return mymif2;
}
