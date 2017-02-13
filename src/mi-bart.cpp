#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using std::cout;
using std::endl;

void printthisfunc(std::vector<double> v)
{
  std::cout << "THIS IS inside printfunction   \n" << endl;
  for (std::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
  {std::cout << *i << endl;} cout<<"outtheloop for" << endl;
}
int flag =0;
size_t m;
double kfac=2.0;
std::vector<double> x;
std::vector<double> y;
std::vector<double> z;
size_t n,p;
double nu;
double lambda;
int bistart=0;
int binum=0;

class init {
public:
  xinfo xi;
  pinfo pi;
  dinfo di;
  std::vector<tree> t;
  std::vector<double> y;
  std::vector<int> z;
  //  std::vector<double>  yy;
  double* allfit; //sum of fit of all trees
  double* r; //y-(allfit-ftemp) = y-allfit+ftemp
  double* ftemp; //fit of current tree
};

//void  init_input(size_t index, init *ip_initial, int* vartype, size_t m, double kfac, int bistart, int binum, size_t n, size_t p, std::vector<double> z,std::vector<double> y, std::vector<double> x) //, double *exp_pt
void  init_input(size_t index, init& ip_initial, int* vartype)
  {
  std::vector<tree> t(m);
  ip_initial.t=t;
  ip_initial.di.x = new double [n*(p-index)];
  ip_initial.allfit = new double [n];

  ip_initial.ftemp=new double [n];
  ip_initial.r=new double [n];

  double miny = INFINITY; //use range of y to calibrate prior for bottom node mu's
  double maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
  double ytemp;

  // cout<< "address of di.x = " << ip_initial.di.x<<endl;
  //
  // cout<< "address of allfit" << ip_initial.allfit<<endl;
  // cout<< "address of first value in allfit" << ip_initial.allfit[0]<<endl;
  //
  //

  if(index==0)
  {
    if(vartype[p-index]==0)
    {
      cout<<"continuous"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.y.push_back(y[i]);
      }
    }
    else
    {
      cout<<"binary"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.y.push_back(z[i*binum+p-bistart]);
        ip_initial.z.push_back(int(y[i]));
      }
    }
    // cout<<"ip_initial.y[0] ="<<ip_initial.y[0]<<endl;
    //
    // cout<<"address for x[0] is =" << &x[0]<<endl;
    //ip_initial.di.x=&x[0];

    for( int s = 0 ; s < n*(p-index); s++)
      ip_initial.di.x[s]=x[s];


    // cout<<"ip_initial.di.x is = (the x[0] address or &x[0] )= "<<ip_initial.di.x<<endl;
    // cout<< "ip_initial.di.x[0] = "<<ip_initial.di.x[0]<<endl;
  }

  else
  {
    if(vartype[p-index]==0)
    {
      cout<<"continuous"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.y.push_back(x[(i+1)*p-index]);
        for(size_t j=0;j<p-index;j++)
        {
          ip_initial.di.x[i*(p-index)+j]=x[i*p+j];
        }
      }
    }
    else
    {
      cout<<"binary"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.z.push_back(int(x[(i+1)*p-index]));
        ip_initial.y.push_back(z[i*binum+p-index-bistart]);
        for(size_t j=0;j<p-index;j++)
        {
          ip_initial.di.x[i*(p-index)+j]=x[i*p+j];
        }
      }
    }
  }
  cout<<"\ndata import ok"<<endl;
  for(size_t i=0;i<n;i++) {
    ytemp=ip_initial.y[i];
    if(ytemp<miny) miny=ytemp;
    if(ytemp>maxy) maxy=ytemp;
    allys.sy += ytemp; // sum of y
    allys.sy2 += ytemp*ytemp; // sum of y^2
  }
  allys.n = n;

  cout <<  "allys.sy : sum of y "<< allys.sy <<  endl;
  cout <<  "allys.sy2 : sum of y^2 "<< allys.sy2 <<  endl;

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
  cout << "ybar,shat: " << ybar << ", " << shat <<  endl;

  cout << "\ny read in:\n";
  cout << "n: " << n << endl;
  cout << "y first and last:\n";
  cout << ip_initial.y[0] << ", " << ip_initial.y[n-1] << endl;
  cout <<  "\n";

  cout << "\nx read in:\n";
  cout << "p: " << p << endl;
  cout << "first row: " <<  ip_initial.di.x[0] << " ...  " << ip_initial.di.x[p-1-index] << endl;
  cout << "last row: " << ip_initial.di.x[(n-1)*(p-index)] << " ...  " << ip_initial.di.x[n*(p-index)-1] << endl;
  cout <<  "\n";

  //x cutpoints
  size_t nc=100; //100 equally spaced cutpoints from min to max.
  makexinfo(p-index,n,&(ip_initial.di.x[0]),ip_initial.xi,nc,vartype);
  cout<<ip_initial.di.x[0]<<endl;

  //dinfo
  for(size_t i=0;i<n;i++)
  {
    ip_initial.allfit[i]=ybar;
  }

  for(unsigned int i=0;i<m;i++)
  {
    ip_initial.t[i].setm(ybar/m); //if you sum the fit over the trees you get the fit.
  }

  //prior and mcmc
  ip_initial.pi.pbd=1.0; //prob of birth/death move
  ip_initial.pi.pb=.5; //prob of birth given  birth/death

  ip_initial.pi.alpha=.95; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  ip_initial.pi.beta=2.0; //2 for bart means it is harder to build big trees.

  if(vartype[p-index]==0)
  {
    ip_initial.pi.tau=(maxy-miny)/(2*kfac*sqrt((double)m));
    ip_initial.pi.sigma=shat;
  }
  else
  {
    ip_initial.pi.tau=3.0/(kfac*sqrt((double)m));
    ip_initial.pi.sigma=1.0;
  }

  cout << "\nalpha, beta: " << ip_initial.pi.alpha << ", " << ip_initial.pi.beta << endl;
  cout << "sigma, tau: " << ip_initial.pi.sigma << ", " << ip_initial.pi.tau << endl;

  cout<<"allfit[n-1] = ";
  cout<<ip_initial.allfit[n-1]<<endl;

  ip_initial.di.n=n;
  cout<<"ip_initial.di.n = ";
  cout<<ip_initial.di.n<<endl;

  ip_initial.di.p=p-index;
  cout<<"ip_initial.di.p = ";
  cout<<ip_initial.di.p<<endl;

  ip_initial.di.y=ip_initial.r; //the y for each draw will be the residual
  //cout << ip_initial.di.y.size() << endl;
  cout<<"ip_initial.r = ";
  cout<<ip_initial.r<<endl;

  ip_initial.di.vartype=vartype;
  cout<<"\n"<<endl;

  // cout << " THE VALUE IS HERE" << endl;
  // cout<<ip_initial.di.x[0]<<endl;
  // for( int s = 0 ; s < n*(p-index); s++)
  //   exp_pt[s] = ip_initial.di.x[s];
  // // cout << " inside init_input experiment array["<<s<<"]"<< exp_pt[s]<<endl;

  // for ( int q=0;q<n;q++){
  //   cout<<"allfit["<<q<<"] = "<<ip_initial.allfit[q]<<endl;
  // }



  // for( int s = 0 ; s < n*(p-index); s++)
  //   cout << " inside init_input experiment array["<<s<<"]"<< &ip_initial.di.x[s]<<endl;
  //


  return;
  }


//---------------------------------
// function to get each mcmc run
//--------------------------------

//void mcmc(init *ip_initial, RNG gen,int* vartype, size_t m, double nu, double lambda, size_t n)
void mcmc(init &ip_initial, RNG gen,int* vartype)
{
  for(size_t j=0;j<m;j++){

    // if (flag == 4|| flag == 1){
    //   cout<<" VALUEs OF ALLFIT inside mcmc function for the j = "<<j<<"th turn are"<<endl;
    //   for(size_t k=0;k<n;k++)
    //     {cout<<ip_initial.allfit[k]<<endl;}
    //
    //     //cout<< "y = "<< ip_initial.y[k] << " allfit = "<< ip_initial.allfit[k]<< " ftemp = "<<ip_initial.ftemp[k]<<endl;
    //
    // }

    //cout << "1 ftemp ~ " << ip_initial.ftemp[0] << endl;
    fit(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.ftemp);


    // if (flag == 4|| flag == 1 ){
    //   {for(size_t k=0;k<n;k++)
    //     cout<< "y = "<< ip_initial.y[k] << " allfit = "<< ip_initial.allfit[k]<< " ftemp = "<<ip_initial.ftemp[k]<<endl;
    //   }
    // }
    //cout << "*************************************"<<endl;

    //cout << "2 ftemp ~ " << ip_initial.ftemp[0] << endl;

    for(size_t k=0;k<n;k++)
    {
      ip_initial.allfit[k] = ip_initial.allfit[k]-ip_initial.ftemp[k];
      ip_initial.r[k] = ip_initial.y[k]-ip_initial.allfit[k];
    }

    // if (flag == 4)
    //   {for(size_t k=0;k<n;k++)
    //     {cout << "R & di.y are: "<< ip_initial.r[k]<<" = "<<ip_initial.di.y[k] <<endl;}}


    bd(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.pi,gen);
    drmu(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.pi,gen);

    //cout << "2 ftemp ~ " << ip_initial.ftemp[0] << endl;
    fit(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.ftemp);
    //cout << "3 ftemp ~ " << ip_initial.ftemp[0] << endl;

    for(size_t k=0;k<n;k++)
    { ip_initial.allfit[k] += ip_initial.ftemp[k];
    }

  }
  flag++;

// cout << " AFTER MCMC, allfit[0] is = "<<ip_initial.allfit[0] << endl;
// for(size_t k=0;k<n;k++)
// {
//   cout << "AFTER MCMC\nallfit["<<k<<"] = "<<ip_initial.allfit[k] <<" and ftemp["<<k<<"] = " << ip_initial.ftemp[k] << endl;
// }


  if(vartype[ip_initial.di.p]==0)
    {
      //draw sigma
      double rss;  //residual sum of squares
      double restemp; //a residual
      rss=0.0;
      for(size_t k=0;k<n;k++)
    	{
    	  restemp=ip_initial.y[k]-ip_initial.allfit[k];
    	  rss += restemp*restemp;
    	}
      //cout <<" restamp n rss = " << restemp << " , " << rss << endl;
      ip_initial.pi.sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    }
  else
    {
      for(size_t k=0;k<n;k++){
        if(ip_initial.z[k]==0){
          ip_initial.y[k]=0.0;
	        while(ip_initial.y[k]>=0.0){
	          ip_initial.y[k]=gen.normal(ip_initial.allfit[k],1);
		        }
	        }
        else{
	        ip_initial.y[k]=0.0;
	        while(ip_initial.y[k]<=0.0){
		        ip_initial.y[k]=gen.normal(ip_initial.allfit[k],1);
		        }
	        }
	      }
    }
  // cout<<"pi.sigma = "<<ip_initial.pi.sigma<<
 // cout <<" *********************** \n"<<endl;

 return;
}

//--------------------------------
// main program
//--------------------------------

//' Multiply a number by one
//' @param x A single integer.
//' @export
//'
// [[Rcpp::export]]
int cpp_bart (NumericVector new_xroot, NumericVector new_yroot, int new_nd, int new_burn,
                int new_m, int new_nu, int new_kfac, int  new_nmissing,
                IntegerVector new_xmissroot, int new_bistart, NumericVector new_vartyperoot,
                NumericVector new_zroot, CharacterVector new_ffname, NumericVector new_lambda, int new_type)
{
  uint seed=99; //random number generation
  RNG gen(seed); //this one random number generator is used in all draws

  y = Rcpp::as< std::vector<double> >(new_yroot);
  cout << "\nprinting y here" << endl;
  printthisfunc(y);


  n = y.size();
  if(n<1)
    {
      cout << "error n<1\n";
      return 1;
    }
  cout << "\ny read in:\n";
  cout << "n: " << n << endl;
  cout << "y first and last:\n";
  cout << y[0] << ", " << y[n-1] << endl;

  x = Rcpp::as< std::vector<double> >(new_xroot); //file to read x from
  cout << "\nprinting xhere" << endl;
  printthisfunc(x);
  cout << "x's size = " << x.size() << endl;

   p = x.size()/n;
  if(x.size() != n*p)
    {
      cout << "error: input x file has wrong number of values\n";
      return 1;
      //cout << "Error2 ";
    }
  cout << "\nx read in:\n";
  cout << "p: " << p << endl;
  cout << "first row: " <<  x[0] << " ...  " << x[p-1] << endl;
  cout << "last row: " << x[(n-1)*p] << " ...  " << x[n*p-1] << endl;

  // size_t m;
  // double kfac=2.0;
  // std::vector<double> x;
  // std::vector<double> y;
  // std::vector<double> z;
  // size_t n,p;
  // double nu;
  // double lambda;
  // int bistart=0;
  // int binum=0;

 size_t nvar = new_nmissing; // # of covariates with missing values
 size_t nd = new_nd; // int atoi (const char * str);
 size_t burn  = new_burn;

  m = new_m;
  kfac = new_kfac;
  nu = new_nu;

  binum = 0;

  std::vector<int> ind_missing = Rcpp::as< std::vector<int> >(new_xmissroot);  //file to read missing indicator from
  cout << "printing ind_missing here" << endl;
   for (std::vector<int>::const_iterator i = ind_missing.begin(); i !=ind_missing.end(); ++i)
   {std::cout << *i << " ";}

  cout << "\nind_missing's size = " << ind_missing.size() << endl;
  cout <<"\nburn, nd, number of trees: " << burn << ", " << nd << ", " << m << endl;
  cout <<"\nlambda =  " << lambda << "nu is = " << nu << " kfac is = " << kfac << endl;
  cout <<"\nnvar: " << nvar<<endl;

  bistart = new_bistart;

  std::vector<int> vart2 = Rcpp::as< std::vector<int> >(new_vartyperoot);  //file to read variable type from
 // cout << "printing vart2 here" << endl;
  //for (std::vector<int>::const_iterator i = vart2.begin(); i !=vart2.end(); ++i)
//  {std::cout << *i << endl;}

  int* vartypeArray=new int[p+1];


   int jj=0;
   //while loop from the original code was repalced by for loop here
   for (std::vector<int>::const_iterator i = vart2.begin(); i != vart2.end(); ++i)
   { vartypeArray[jj]=*i;
       if(jj>=bistart){
       if(*i==1){binum++;}
     }
     jj++;
   }

   // for (int s = 0; s < p+1; s++){
   //  cout << "the content of the array = " <<vartypeArray[s] << endl;
   // }

   z = Rcpp::as< std::vector<double> >(new_zroot);
//   cout << "printing z here" << endl;
//   printthisfunc(z);

   cout << "z's size is = " << z.size() << endl;

    std::string new_ffname_conv = Rcpp::as<std::string>(new_ffname);
    std::string filename="";
    filename.append(new_ffname_conv);

  init *all_initial = new init[nvar+1]; //right

  for(size_t i=0;i<nvar+1;i++)
    {
    cout << "calling init input for "<<i<<"th time" << endl;

    //double *expermient_pointer = new double[n*(p-i)];

    //for( int s = 0 ; s < n*(p-i); s++) cout << "1 experiment array["<<s<<"]"<< expermient_pointer[s]<<endl;

    init_input(i,all_initial[i],vartypeArray);//,  m, kfac, bistart, binum, n, p, z, y, x);


    //for (std::vector<int>::const_iterator i = vart2.begin(); i !=vart2.end(); ++i)
    //  {std::cout << *i << endl;}

   // double b = all_initial[i].di.x.size;
    // int b = n*(p-i);
    // cout << " for reference n = " << n << " p is = " << p << " i is = "<< i <<endl;
    // for ( size_t j = 0; j < b; j++)
    //   {//cout << " for reference n = " << n << " p is = " << p << " i is = "<< i <<endl;
    //   cout << " all_initial[i].di.x["<<j<<"] is "<<all_initial[i].di.x[j] << endl;
    //  // cout << " \n2 experiment array["<<j<<"]"<< expermient_pointer[j]<<endl;
    //   }

    //cout << "\nfirst and last row of x"<<endl;
    //cout << "first row: " <<  all_initial[i].di.x[0] << " ...  " << all_initial[i].di.x[p-1-i] << endl;
    //cout << "last row: " << all_initial[i].di.x[(n-1)*(p-i)] << " ...  " << all_initial[i].di.x[n*(p-i)-1] << endl;


    }

  std::ofstream mif(filename.c_str());

  // //mcmc
  double pro_prop=0.0; //proposal probability of MH min(1,alpha)
  double pro_propnum=1.0;//numerator of proposal prob
  double pro_propdenom=1.0;//denominator of proposal prob
  double prop=0.0;// proposal for mh
  dinfo di_temp; //x used to get fitted value
  double* ppredmean=new double[p];
  double* fpredtemp=new double[1];
  di_temp.n=1; di_temp.y=0;
  std::vector<double> tempx;
  for(size_t i=0;i<=p;i++) tempx.push_back(0.0);
  size_t l=0;
  cout << "\nMCMC:\n";
  clock_t tp;
  tp = clock();

  int type = new_type;
  cout << "you have type =" << type<< endl;
  cout << " \n ********************* THE BIGGEST LOOP IS HERE**************************" <<endl;
  for(size_t i=0;i<(nd+burn);i++)
  {
    cout << "i: " << i << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++ i starts +++++++++++++++++++++++++++++++++++++++++++++\n" << endl;

    cout <<"first internal loop starts"<<endl;
    for(size_t j=0;j<nvar+1;j++){
      lambda = (new_lambda[j]);
      cout << " for object number ="<<j<<endl;
      cout << "lambda ["<< j <<"] = " <<lambda << endl;
      //cout << "BEFORE MCMC = "<< all_initial[j].allfit[0] <<endl;

      mcmc(all_initial[j], gen,vartypeArray);//, m, nu, lambda, n);
      cout << " Done first loop, ie mcmc done, for object number ="<<j<<endl;
    }

    // cout <<"\nsecond internal loop starts"<<endl;
     cout << "BEFORE the i= "<<i<<"th loop ="<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<endl;


    for(size_t j=0;j<n;j++){
      cout << " VALUE OF j is = " << j <<endl;
      for(size_t k=0;k<=nvar;k++){
        cout << " JUST TO CHECK THE K " << k <<endl;

        // for(size_t s = 0; s < n; s++)
        //   cout<< "inside last while loop MINuS all_initial is "<<all_initial[k].allfit[s] << endl;

        if(ind_missing[j*(nvar+1)+k]==1){
  	      for(size_t j1=0;j1<p;j1++){
  	        tempx[j1]=all_initial[0].di.x[j*p+j1];
  	        }

            if(vartypeArray[p-k]==0){
  		      prop=gen.normal(all_initial[k].allfit[j],all_initial[k].pi.sigma);
              cout << "here" <<endl;
              pro_propdenom=1.0;
		        pro_propnum=1.0;
		        }
            else{
            prop=double(abs(all_initial[k].z[j]-1));
              cout << "there" <<endl;
              pro_propnum=fabs(prop-phi(-all_initial[k].allfit[j]));
            pro_propdenom=fabs(double(all_initial[k].z[j])-phi(-all_initial[k].allfit[j]));
            }
            cout << "here is prop = " << prop << endl;
            cout << "before while::  all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<endl;
            tempx[p-k]=prop;//change the element to proposal

            // cout << "NOOoooooooooooooooooooooooooooooooooooooo"<<endl;
            // cout << "K IS EQUAl TO = " << k << endl;
            // cout << "herererererererererererererererererererer" <<endl;
            //
		        l=0;
  		      while(l<k){
  		        cout<<"In first while loop"<<endl;
  		         //cout << "the value of l is =" << l <<endl;
  		        if(vartypeArray[p-l]==0){
  		          pro_propdenom=pro_propdenom*pn(all_initial[l].y[j],all_initial[l].allfit[j],pow(all_initial[l].pi.sigma,2));
                }
  		        else{
			          pro_propdenom=pro_propdenom*fabs(double(all_initial[l].z[j])-phi(-all_initial[l].allfit[j]));
			          }
              di_temp.p=p-l; di_temp.x=&tempx[0];
              di_temp.vartype=vartypeArray;
  		        ppredmean[l]=0.0;
  		        for(size_t q=0;q<m;q++){
  		          //if ( i ==0 )
                 // cout << " i is =" << i<< "di_temp is = "<< di_temp.x<<endl;
  		            //cout << " i is =" << i << "\n"<<" all_initial[l].t["<<q<<"] "<< all_initial[l].t[q] <<endl;
  		         // cout << "the value of q is =" << q <<endl;
  			        fit(all_initial[l].t[q],all_initial[l].xi,di_temp,fpredtemp);
        			  ppredmean[l] += fpredtemp[0];
        			  // cout <<" ppredmean is = " <<ppredmean[l]<< endl;
  			        }

  		          //cout <<"FIRST ppredmean for k = "<<k<<" is  = " << ppredmean[l]<< endl;

              if(vartypeArray[p-l]==0){
  		          pro_propnum=pro_propnum*pn(all_initial[l].y[j],ppredmean[l],pow(all_initial[l].pi.sigma,2));
                }
              else{
		            pro_propnum=pro_propnum*fabs(double(all_initial[l].z[j])-phi(-ppredmean[l]));
                }
              //cout <<"SECOND ppredmean for k = "<<k<<" is  = " << ppredmean[l]<< endl;
		            l++;
		        }
  		      cout << "after first while" << endl;
  		      cout << "after while::  all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<endl;

            if(vartypeArray[p-k]==0){
  		          pro_prop=std::min(1.0,pro_propnum/pro_propdenom);
		            }
            else{
                pro_prop=pro_propnum/(pro_propnum+pro_propdenom);
		            }
            double gu = gen.uniform();
  //           cout << "gen.uniform() = " << gu <<endl;
  //           cout << "pro_prop for K = ["<<k<<"] is "<<pro_prop << endl;
   		      if(gu<pro_prop){
   		        cout << "\n\ncame in the last if loop"<<endl;
  		          if(vartypeArray[p-k]==0){
  		            cout << "YOOP"<<endl;
  		            all_initial[k].y[j]=prop;
		              }
  		          else{cout << "NOOOP"<<endl;
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
  		          // for(size_t s = 0; s < n; s++)
  		          //   cout<< "inside last while loop ZERO all_initial is "<<all_initial[k].allfit[s] << endl;

  		          cout << "place2 "<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<"\n\n\n"<<endl;

                l=0;
		            while(l<k){
		              cout << "Final while loop"<<endl;
		              cout << "position for x = " <<j*(p-l)-k+p <<endl;
		              cout << "address of all_initial di.x[j*(p-l)-k+p = "<<j*(p-l)-k+p<<"] "<<&all_initial[k].di.x[j*(p-l)-k+p]<<endl;


		              all_initial[l].di.x[j*(p-l)-k+p]=prop;
		              cout << "place3 "<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<"\n"<<endl;

		              // for(size_t s = 0; s < n; s++)
		              //   cout<< "inside last while loop first all_initial is "<<all_initial[k].allfit[s] << endl;
		              cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! infinity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!= "<<j << " "<< ppredmean[l]<<endl;
		              all_initial[l].allfit[j]=ppredmean[l];
		              //for(size_t s = 0; s < n; s++)
		              //cout<< "inside last while loop second all_initial is "<<all_initial[k].allfit[s] << endl;
		              cout << "place4 "<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<"\n"<<endl;

                  l++;
                }
		            mif<<prop<<endl; //but prop never really changes so why so much work above.....???
            }





  		      else{
  		        //cout << "no"<<endl;
  		          if(vartypeArray[p-k]==0) {
			            mif<<all_initial[k].y[j] <<endl;
  		          }
  		          else{
  		            mif<<all_initial[k].z[j] <<endl;
  		          }
  		      }

  		      cout << "in the IF loop "<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[63]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<"\n"<<endl;
        }
        // cout << "harish"<<endl;
        // for(size_t kk=0;kk<n;kk++)
        //   cout<<all_initial[k].allfit[kk] <<endl;
      }
  	}


    for(size_t kk=0;kk<n;kk++)
      cout<<all_initial[1].allfit[kk] <<endl;

    //cout << "AFTER the i= "<<i<<"th loop"<<" all_initial + 0 = "<<all_initial[0].allfit[0] <<" all_initial + 1 = "<<all_initial[1].allfit[0]<<" all_initial + 2 = "<<all_initial[2].allfit[0]<<"\n"<<endl;
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ i ends @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << endl;

  }
  tp=clock()-tp;
  double thetime = (double)(tp)/(double)(CLOCKS_PER_SEC);
  cout << "time for loop: " << thetime << endl;
  cout << "you had type =" << type<< endl;
  cout << "aarti5feb";
  return 0;
}
