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
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using namespace Eigen;


// static void printthisfunc(std::vector<double> v)
// {
//     std::cout << "You are inside the Printfunction   \n" << endl;
//     for (std::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
//     {std::cout << *i << endl;}
//     cout<<"out of the print function";
// }

static int flag;



static size_t m;
static double kfac;
static std::vector<double> x,y;
static std::vector<double> z;
static size_t n,p;
static double nu;
static double lambda;
static int bistart;
static int binum;


class init {
public:

  // ~init(){
  //   t.clear();
  //   y.clear();
  //   z.clear();
  //   delete[] allfit;
  //   delete[] r;
  //   delete[] ftemp;
  // }
  xinfo xi;
  pinfo pi;
  dinfo di;
  std::vector<tree> t;
  std::vector<double> y;
  std::vector<int> z;
  double* allfit; //sum of fit of all trees
  double* r; //y-(allfit-ftemp)  = y-allfit+ftemp
  double* ftemp; //fit of current tree
};

//-----------------------------------------
// function to get initial data
//-----------------------------------------

static void init_input(size_t index, init *ip_initial, int* vartype)
  {
    std::vector<tree> t(m);
    ip_initial->t=t;
    ip_initial->di.x = new double [n*(p-index)];
    ip_initial->allfit = new double [n];
    ip_initial->ftemp=new double [n];
    ip_initial->r=new double [n];

   //read in data
   //read y NOTE: we assume y is already centered at 0.
   double miny = INFINITY; //use range of y to calibrate prior for bottom node mu's
   double maxy = -INFINITY;
   sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
   double ytemp;

   // cout << "//***********************************************ONE***************************************//" <<endl;
   // cout << "index is equal to = " << index << endl;
   if(index==0)
     {
       // cout << "index was 0" << endl;
       // cout << "***********"<< endl;
       // cout << "[p-index] =" <<p-index <<endl;
       // cout << "vartype[p-index] =" <<vartype[p-index] <<endl;
       // cout << "***********"<<endl;
      if(vartype[p-index]==0)
	     {
    	  // cout<<"continuous"<<endl;
    	   for(size_t i=0;i<n;i++)
    	   {
	        ip_initial->y.push_back(y[i]);
	       }
	     }
      else
	     {
    // 	    cout<<"binary"<<endl;
    //       cout << "binum1 is" << binum << endl;
    //       cout << "bistart1 is" << bistart << endl;
          for(size_t i=0;i<n;i++)
    	     {


     	       ip_initial->y.push_back(z[i*binum+p-bistart]);
    	       ip_initial->z.push_back(int(y[i]));
    	     }
        }
      // cout<<"ip_initial->y[0] = ";
      // cout<<ip_initial->y[0]<<endl;

     // cout<<"ip_initial->di.x = ";  //WHAT IS THIS??
    	ip_initial->di.x=&x[0];

    // 	cout<<"ip_initial->di.x = ";
    //   cout<<ip_initial->di.x[0]<<endl;
   }
   else
   {
       // cout << "index was not = 0" <<endl; cout << "***********" << endl;
       // cout << "[p-index]" <<p-index <<endl;
       // cout << "vartype[p-index]" <<vartype[p-index] <<endl;  cout << "***********"<<endl;
       if(vartype[p-index]==0)
    	 {
    	   //cout<<"continuous"<<endl;
    	   for(size_t i=0;i<n;i++)
  	     {
  	       ip_initial->y.push_back(x[(i+1)*p-index]);
  	       for(size_t j=0;j<p-index;j++)
      		 {
      		   ip_initial->di.x[i*(p-index)+j]=x[i*p+j];
      		 }
         }
    	 }
       else
    	 {
    // 	   cout<<"binary"<<endl;
    //      cout << "binum2 is" << binum << endl;
    //      cout << "bistart2 is" << bistart << endl;

    	   for(size_t i=0;i<n;i++)
  	     {
    	     //cout << " i = " << i << " , binum = " << binum << " , p = " << p << " , index = "<< index << " , bistart = " << bistart<< endl;
  	       ip_initial->z.push_back(int(x[(i+1)*p-index]));
    	     //cout << "i*binum+p-index-bistart = " << i*binum+p-index-bistart << endl;
    	     //cout << "internal Z is here :: " << z[i*binum+p-index-bistart] << "  ";
  	       ip_initial->y.push_back(z[i*binum+p-index-bistart]); //  y is getting wrong z values.. cause z's index is skipping one value each time
  	       for(size_t j=0;j<p-index;j++)
      		 {
      		   ip_initial->di.x[i*(p-index)+j]=x[i*p+j];
      		 }
  	     }
    	 }
     }
// cout << "//**********************************************TWO*********************************************//" <<endl;

//cout<<"\ndata import ok"<<endl;
//cout <<"Size of y = " <<ip_initial->y.size() << " and Size of z = "; cout << ip_initial->z.size();


 if (index==3){

// cout<<"Printing Y's here" <<endl;
// printthisfunc(ip_initial->y);
// cout<<"Printing z's here" <<endl;
//printthisfunc(ip_initial->z);  // wont work cause ip_initial->z is int and printthisfunc is for double
// for (std::vector<int>::iterator i = ip_initial->z.begin(); i !=ip_initial->z.end(); ++i)
// {std::cout << *i << endl;}

 }

//cout << "\n//***********************************************THREE**********************************************//" <<endl;

   for(size_t i=0;i<n;i++){
     ytemp=ip_initial->y[i];
    // cout <<" ytemp = " << ytemp;

   if(ytemp<miny) miny=ytemp;
      if(ytemp>maxy) maxy=ytemp;
      allys.sy += ytemp; // sum of y
     // cout <<" .... The Sum of all ys is here = " << allys.sy  << endl;
     allys.sy2 += ytemp*ytemp; // sum of y^2
   }


   allys.n = n;
   // cout << "\ny read in:\n";
   // cout << "n: " << n << endl;
   // cout << "y first and last:\n";
   // cout << ip_initial->y[0] << ", " << ip_initial->y[n-1] << endl;
   //
   // cout <<  "\nallys.sy : sum of y "<< allys.sy <<  endl;
   // cout <<  "allys.sy2 : sum of y^2 "<< allys.sy2 <<  endl;

   double ybar = allys.sy/n; //sample mean
   double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
//   cout << "ybar,shat: " << ybar << ", " << shat <<  endl;

   //read x
   //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
   // cout << "\nx read in:\n";
   // cout << "p: " << p << endl;
   // cout << "first row: " <<  ip_initial->di.x[0] << " ...  " << ip_initial->di.x[p-1-index] << endl;
   // cout << "last row: " << ip_initial->di.x[(n-1)*(p-index)] << " ...  " << ip_initial->di.x[n*(p-index)-1] << endl;

   //x cutpoints
   size_t nc=100; //100 equally spaced cutpoints from min to max.
   makexinfo(p-index,n,&(ip_initial->di.x[0]),ip_initial->xi,nc,vartype);
//   cout<<ip_initial->di.x[0]<<endl;
   //dinfo
   for(size_t i=0;i<n;i++)
     {
       ip_initial->allfit[i]=ybar;
     // cout<<"allfit[i] has the following values = ";
     // cout<<ip_initial->allfit[i]<<endl;
     }



   for(unsigned int i=0;i<m;i++)
   {
     ip_initial->t[i].setm(ybar/m);
   }

   // for(unsigned int i=0;i<m;i++)
   // {
   //   cout << "here is the setm = "<< ip_initial->t[i].getm() << endl;
   // }

  //prior and mcmc
   ip_initial->pi.pbd=1.0; //prob of birth/death move
   ip_initial->pi.pb=.5; //prob of birth given  birth/death

   ip_initial->pi.alpha=.95; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   ip_initial->pi.beta=2.0; //2 for bart means it is harder to build big trees.

   if(vartype[p-index]==0)
     {
       ip_initial->pi.tau=(maxy-miny)/(2*kfac*sqrt((double)m));
       ip_initial->pi.sigma=shat;
}
   else
     {
       ip_initial->pi.tau=3.0/(kfac*sqrt((double)m));
       ip_initial->pi.sigma=1.0;
     }


   // cout << "\nalpha, beta: " << ip_initial->pi.alpha << ", " << ip_initial->pi.beta << endl;
   // cout << "sigma, tau: " << ip_initial->pi.sigma << ", " << ip_initial->pi.tau << endl;
   //
   // cout<<"allfit[n-1] = ";
   // cout<<ip_initial->allfit[n-1]<<endl;

   ip_initial->di.n=n;
   // cout<<"di.n = ";
   // cout<<ip_initial->di.n<<endl;

   ip_initial->di.p=p-index;
   // cout<<"di.p = ";
   // cout<<ip_initial->di.p<<endl;

   ip_initial->di.y=ip_initial->r; //the y for each draw will be the residual
   // cout<<"r = ";
   // cout<<ip_initial->r<<endl;
   //
   ip_initial->di.vartype=vartype;
   //cout<<"\n"<<endl;

   // for(unsigned int i=0;i<m;i++)
   // {
   //   cout << "here is the setm 2 = "<< ip_initial->t[i].getm() << endl;
   // }


  return;
}


//---------------------------------
// function to get each mcmc run
//--------------------------------
//mcmc(all_initial[j], gen,vartype);
static void mcmc(init *ip_initial, RNG  gen,int* vartype){

for(size_t j=0;j<m;j++)
{
 // if (j==0) cout << "1ftemp = :"<<ip_initial->ftemp[0]<<endl;

  fit(ip_initial->t[j],ip_initial->xi,ip_initial->di,ip_initial->ftemp);
  //cout << "here is the setm 3 = "<< ip_initial->t[j].getm() << endl;

 //if (j==0)cout << "2ftemp = :"<<ip_initial->ftemp[0]<<endl;




  for(size_t k=0;k<n;k++)
  {
    //cout << "here is  ip_initial" << ip_initial->allfit[k] << endl;
    ip_initial->allfit[k] = ip_initial->allfit[k]-ip_initial->ftemp[k];
    ip_initial->r[k] = ip_initial->y[k]-ip_initial->allfit[k];
  }

  bd(ip_initial->t[j],ip_initial->xi,ip_initial->di,ip_initial->pi,gen);
  //cout << "here is the setm 4 = "<< ip_initial->t[j].getm() << endl;

  drmu(ip_initial->t[j],ip_initial->xi,ip_initial->di,ip_initial->pi,gen);
  //cout << "here is the setm 5 = "<< ip_initial->t[j].getm() << endl;

//if (j==0)cout << "2ftemp = :"<<ip_initial->ftemp[0]<<endl;

  fit(ip_initial->t[j],ip_initial->xi,ip_initial->di,ip_initial->ftemp);
  //cout << "here is the setm 6 = "<< ip_initial->t[j].getm() << endl;

// if (j==0)cout << "3ftemp = :"<<ip_initial->ftemp[0]<<endl;

  for(size_t k=0;k<n;k++) ip_initial->allfit[k] += ip_initial->ftemp[k];
}

  if(vartype[ip_initial->di.p]==0)
    {
      //draw sigma
      double rss;  //residual sum of squares
      double restemp; //a residual
      rss=0.0;
      for(size_t k=0;k<n;k++)
	    {
  	  restemp=ip_initial->y[k]-ip_initial->allfit[k];
  	  rss += restemp*restemp;
  	  }
  //cout <<" restamp n rss = " << restemp << " , " << rss << endl;
      ip_initial->pi.sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    }
  else
    {
      for(size_t k=0;k<n;k++)
	      {
          if(ip_initial->z[k]==0)
  	      {
  	      ip_initial->y[k]=0.0;
    	      while(ip_initial->y[k]>=0.0)
    		    {
		           ip_initial->y[k]=gen.normal(ip_initial->allfit[k],1);
		        }
	        }
          else
    	    {
    	      ip_initial->y[k]=0.0;
    	      while(ip_initial->y[k]<=0.0)
    		    {
      		  ip_initial->y[k]=gen.normal(ip_initial->allfit[k],1);
      		  }
	        }
	      }
    }
  //cout<<"pi.sigma = "<<ip_initial->pi.sigma<<endl;
  return;
}

//--------------------------------
// main program
//--------------------------------

//' Multiply a number by two
//' @param x A single integer.
//' @export
//'
// [[Rcpp::export]]
int cpp_bart_y (NumericVector new_xroot, NumericVector new_yroot, int new_nd, int new_burn,
                int new_m, int new_nu, int new_kfac, int  new_nmissing,
                IntegerVector new_xmissroot, int new_bistart, NumericVector new_vartyperoot,
                NumericVector new_zroot, CharacterVector new_ffname, NumericVector new_lambda, int new_type)
{

  //cmd = paste(cmd,xroot,yroot,"0",nd,"0.0",
  //            burn,m,nu,kfac,nmissing,xmissroot,bistart,vartyperoot,zroot,ffname)
  // std::cout.precision(5);
  // std::cout  << std::fixed;


  cout<< "2017 build"<<endl;
  //cout << "************My code************\n";
  // double xtemp; //used to read in x and y without knowing how many there are.
  // double ytemp;


  //random number generation
  uint seed=99;
//  RNG * gen = new RNG(seed); //this one random number generator is used in all draws
  RNG gen(seed);

 // gen.zigset();

  //read in data
  //read y NOTE: we assume y is already centered at 0.

  //file to read y from
 // std::string new_yroot_conv = Rcpp::as<std::string>(new_yroot);
 // const char *somethingY = new_yroot_conv.c_str(); //Returns a pointer to the c-string representation of the string object's value.
  //std::ifstream yf(somethingY); //the constructor for an ifstream takes a const char*, not a string

 // cout<<"somethingY   =  " << somethingY<<endl;

 // cout << "AARTI HERE IS THE INSIDE OF THE Y FILE" << endl;
 //failed attempt to read in the file y when it was a ifstream object and thus was read into the cpp code
  // std::string line;
  // if (yf.is_open())
  // {
  //   while ( yf.good() )
  //   {
  //     line = getline(yf,line);
  //     cout << line << endl;
  //   }
  //   yf.close();
  // }
  //
  // else cout << "Unable to open file";
  // y.clear();
  y = Rcpp::as< std::vector<double> >(new_yroot);

  //...TO READ WHAT IS INSIDE Y vector
  // for (std::vector<double>::const_iterator i = y.begin(); i != y.end(); ++i)
  // {std::cout << *i << ' ';}
  //cout <<" ...TO READ WHAT IS INSIDE Y vector" <<endl;
  //printthisfunc(y);
  //cout << "DONE Y" <<endl;


 ////std::vector<double> y;
 // while(yf >> ytemp)
 //   {
 //    y.push_back(ytemp);
 //   } // void push_back (const value_type& val);

  n = y.size();
 // cout<<"\nValue of n is " << n << endl;
  if(n<1)
  {cout << "error y's size is <1\n";
    return 1;
  }

  // cout<< "Y included. Linear regression. Diffuse prior.";
  // cout << "\ny read in:\n";
  // cout << "n: " << n << endl;
  // cout << "y first and last:\n";
  // cout << y[0] << ", " << y[n-1] << endl;

//  cout << ">>>>>>>>>>>>>>>>>>> Y ENDS HERE <<<<<<<<<<<<<\n"<< endl;



  //read x
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on
  //file to read x from
 // std::ifstream xf(new_xroot.c_str());
  //file to read x from
  // std::string new_xroot_conv = Rcpp::as<std::string>(new_xroot);
  // const char *somethingX = new_xroot_conv.c_str(); //Returns a pointer to the c-string representation of the string object's value.
  // std::ifstream xf(somethingX); //the constructor for an ifstream takes a const char*, not a string
  //
  // //std::vector<double> x;
  // while(xf >> xtemp)
  //   {
  //   x.push_back(xtemp);
  //   }


  //********** Read the 3 stars  *** comment lines together!!
  // x.clear();
  x = Rcpp::as< std::vector<double> >(new_xroot); // ***aparently x is not a 2d vector.

  // cout <<" ...TO READ WHAT IS INSIDE X vector" <<endl;
  // printthisfunc(x);
  // cout << "DONE X" <<endl;


  // ***it is a 1-d vector.
 // int s = x[1].size(); // *** therefore this is wrong


 // *** So, this is importnat>>>That leaves, that ONLY the sequence of x matters.<<<<<<<<<<<<<<<<< not its number rows and columns

//   int s = x.size();
//  cout << "x's size = " << s << endl;
  p = x.size()/n; //  *** And p is regenerated anyways!! and we dont take p as = numb of col of incoming x,,, cause that is not possible in the vector form of  x here

  // *** and thats why when I sent transX, it is acceptable....


  // std::cout << "THIS IS X   \n" << endl;
  // for (std::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i)
  // {std::cout << *i << ' ';}

  if(x.size() != n*p){
    cout << "error: input x file has wrong number of values\n";
    return 1;
    }

  // cout << "\nx read in:\n";
  // cout << "p: " << p << endl;
  // cout << "first row: " <<  x[0] << " ...  " << x[p-1] << endl;
  // cout << "last row: " << x[(n-1)*p] << " ...  " << x[n*p-1] << endl;
 // cout << ">>>>>>>>>>>>>>>>>>> X ENDS HERE <<<<<<<<<<<<<\n" << endl;




  //optionally read in additional arguments
 // size_t burn = 100; //number of mcmc iterations called burn-in
 // size_t nd = 1000; //number of mcmc iterations
  //m=200; //number of trees in BART sum
  //lambda = 1.0; //this one really needs to be set
  //nu = 3.0;
  //size_t nvar=0; // # of covariates with missing values

  m = 200;
  nu = 3.0;
  kfac = 2.0;
  m = new_m;
  nu = new_nu;
  kfac = new_kfac;


  //const char *nd_p = new_nd.c_str()
  //const char *something = new_nd.c_str();

  // Because previously used 'argv' has pointers to access the array it stored the input argument values


  size_t burn = 100;
  size_t nd = 1000;
  size_t nvar=0;

   burn  = new_burn;
   nd = new_nd; // int atoi (const char * str); //lambda = std::atof(new_lambda);
   nvar = new_nmissing;

   bistart=0;
   binum=0;

  //std::ifstream mindf(new_xmissroot);  //file to read missing indicator from

  // std::string new_xmissroot_conv = Rcpp::as<std::string>(new_xmissroot);
  // const char *somethingXmiss = new_xmissroot_conv.c_str();
  // std::ifstream mindf(somethingXmiss);



  std::vector<int> ind_missing = Rcpp::as< std::vector<int> >(new_xmissroot);
  // int mtemp;
  // while(mindf>> mtemp)
  //  {
  //   ind_missing.push_back(mtemp);
  //  }

  // cout << "\nind_missing's size = " << ind_missing.size() << endl;
  //
  // cout <<"\nburn, nd, number of trees: " << burn << ", " << nd << ", " << m << endl;
  // cout <<"\nlambda =  " << lambda << "nu is = " << nu << " kfac is = " << kfac << endl;
  // cout <<"\nnvar: " << nvar<<endl;

  //new_bistart
  bistart = new_bistart;
  int* vartypeArray = new int[p+1];
  // std::string new_vartyperoot_conv = Rcpp::as<std::string>(new_vartyperoot);
  // const char *something_vartyperoot = new_vartyperoot_conv.c_str();
  // std::ifstream vart(something_vartyperoot);

  std::vector<int> vart2 = Rcpp::as< std::vector<int> >(new_vartyperoot);

  //int vartemp;
  int jj = 0;
  for (std::vector<int>::const_iterator i = vart2.begin(); i != vart2.end(); ++i)
  {
   // cout << "jj is here" << jj <<endl;
    vartypeArray[jj]=*i;
    if(jj>=bistart)
    {
      if(*i==1)
      {
        binum++;
      }
    }
    jj++;
  }
  // cout << " stupid binum is here" << binum<< endl;
  // while(vart >> vartemp)
  // {
  // }
  //  jj=0;
  //  while(vartype[jj]==0) {bistart++;jj++;}



  // std::ifstream bi(new_zroot);

  // std::string new_zroot_conv = Rcpp::as<std::string>(new_zroot);
  // const char *something_new_zroot = new_zroot_conv.c_str();
  // std::ifstream bi(something_new_zroot);
  //
  // double bitemp;
  // jj=0;
  // //std::vector<double> z;
  //
  // while(bi>>bitemp)
  // {
  //   z.push_back(bitemp);
  //   //   cout<<"z["<<jj<<"]"<<bitemp<<endl;
  //   jj++;
  // }

  // z.clear();
  z = Rcpp::as< std::vector<double> >(new_zroot);

//  cout << "z's size = " << z.size() << endl;
  // cout << " ...TO READ WHAT IS INSIDE Z vector" <<endl;
  // printthisfunc(z);
  // cout << "DONE Z"<<endl;


 // cout << "jj was dropped by  me" << endl;

  //string studentID[numStudents]; //wrong
  //string *studentID = new string[numStudents]; //right

  // Pointer Creation!!
  init *all_initial = new init[nvar+1]; //right // all_initial is an array name. it has nvar+1 values. Each value is of type init

  // all_initial = [ init 1, init 2, init 3, init4 ]

  //init all_initial[nvar+1]; //Commenting it for now!
  // for(size_t i=0;i<nvar+1;i++)
  // {
  //   (all_initial + i)->ftemp=new double [n];
  //   (all_initial + i)->r=new double [n];
  // }
  //


  for(size_t i=0;i<nvar+1;i++)
    {
     // cout << "HERE is where init_input wil be called for "<<i <<"th time" << endl;
      //cout<<"\n function"<<i+1<<endl;
      //cout << "p is" << p << endl;
      init_input(i,(all_initial + i),vartypeArray);
      // cout<<all_initial[i].r<<endl;
      // cout<<"first and last ys:\n"<<endl;
      // cout <<all_initial[i].y[0]<<"..."<<all_initial[i].y[n-1]<<endl;
      // cout << "\nfirst and last row of x"<<endl;
      // cout << "first row: " <<  all_initial[i].di.x[0] << " ...  " << all_initial[i].di.x[p-1-i] << endl;
      //  cout << "last row: " << all_initial[i].di.x[(n-1)*(p-i)] << " ...  " << all_initial[i].di.x[n*(p-i)-1] << endl;

    }

  // for (size_t i=0;i<n;i++){
  //   cout << "THIS place "<< (all_initial + 3)->allfit[i] <<endl;
  // }


   //storage for ouput
   //missing value


   //std::string name1(new_ffname);
   std::string new_ffname_conv = Rcpp::as<std::string>(new_ffname);
   //const char *something_new_ffname = new_ffname_conv.c_str();
   //std::ifstream name1(something_new_ffname);
   //cout<<name1;
   std::string filename="";
   filename.append(new_ffname_conv);
  // cout<<"here is the filename"<<filename <<endl;

   std::ofstream mif(filename.c_str());//file to save imputed missing value
   //std::ofstream miftemp("/Users/as82986/BART/miftemp.txt");
   //std::ofstream accpt("/home/dandan/Downloads/test/accpt.txt");//file to save imputed missing value
   //std::ofstream bsdf("/home/dandan/Downloads/test/bart-sd.txt"); //note that we write all burn+nd draws to this file.



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


   // //mcmc
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
  //VectorXd yy(y.data());
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > yy(y.data(),n,1);
  VectorXd ee(n);
  double* sigsq=new double[1];
  Map<Eigen::Matrix<double,1,1> > sigsq1(sigsq);
  MatrixXd cho(p+1,p+1);
  MatrixXd xxinv(n,n);
  //VectorXd rbeta(p+1);
  std::vector<double> rnor;
  for(size_t i=0;i<=p;i++) rnor.push_back(0.0);
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > rbeta(rnor.data(),p+1,1);
  double* meantemp1=new double[1];
  double* meantemp2=new double[1];
  Map<Eigen::Matrix<double,1,1> > meantemp11(meantemp1);
  Map<Eigen::Matrix<double,1,1> > meantemp21(meantemp2);
  VectorXd xxtemp(p+1);

  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > xx1( (all_initial + 0)->di.x,n,p);
  xx<<x1,xx1;
   // cout << "\nsize of xx "<<xx.size() << endl;
   // cout << "\nsize of yy "<<yy.size() << endl;
   // cout <<"memememmememem";
   // cout<<"\n this is xx="<<xx;
  beta=(xx.transpose()*xx).inverse()*xx.transpose()*yy;
  // cout<<" \n\nbeta"<<beta<<endl;
  // end change 1


  //double temp_a;
  di_temp.n=1;
  di_temp.y=0;
  std::vector<double> tempx;
  for(size_t i=0;i<=p;i++) tempx.push_back(0.0);
  Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> > xxtemp1(tempx.data(),p,1);
  xxtemp<<1,xxtemp1;
//cout << "\nMCMC:\n"<<endl;
//cout << "xxtemp1="<<xxtemp1<<endl;
//cout <<"xxtemp="<<xxtemp<<endl;
  clock_t tp;
  tp = clock();
  size_t l=0;
 //  int numberoftimes= nd+burn;
 //
 //   //(all_initial + k)->allfit[j]
 //   cout << "THIS place 22222222  = "<< (all_initial + 1)->allfit[1] <<endl;
 //
 //    // for (size_t i=0;i<n;i++){
 //    //   cout << "THIS place 3333333  = "<< (all_initial + 1)->allfit[i] <<endl;
 //    // }
 //
 // cout<<"\nthe numberof time the loop wil run is here "<< numberoftimes <<endl;


flag = 0;
int type = new_type;
cout << "Value of 'type' is =" << type<< endl;
//  cout << " \n ********************* THE BIGGEST LOOP IS HERE**************************" <<endl;
  for(size_t i=0;i<(nd+burn);i++)
    {
      if(i%100==0) cout << "i: " << i << endl;

      for(size_t j=1;j<nvar+1;j++){
	       lambda= (new_lambda[j]); if(i == 0){
          //cout << "j is = " << j << endl;
          //cout << "lambda = "<< j <<" = " <<lambda << endl;
          //cout << "THIS place xxxxxx  = "<< (all_initial + 1)->allfit[1] <<endl;
          }
	       //cout <<"calling MCMC " << "j = " <<j << "i = " <<i<<endl;
      	 mcmc((all_initial + j), gen,vartypeArray);
        }

      //new change2 -- posterior of beta and sigma

      // if(i == 0){
      //   cout <<" wewewewewewewewewee";
      //  // cout<<"\n this is xx =\n"<<xx;
      // }
      // xx<<x1,xx1;
      // if(i == 0){
      //   cout <<" yeyeyeyeyyeyeyyeyeyyeyeyeyyee";
      //  // cout<<"\n this is xx =\n"<<xx;
      // }

      // cout << "\nsize of xx 2 "<<xx.size();
      // cout << " size of yy 2 "<<yy.size();
      // cout <<" yeyeyeyeyyeyeyyeyeyyeyeyeyyee";
      // cout<<"\n this is xx =\n"<<xx;
      ee=yy-xx*beta;
      //cout<<"\nee is = "<<ee;
    	sigsq1=(ee.transpose()*ee);
    	sigsq[0]=sigsq[0]/gen.chi_square(double(n-p-1));
    	// cout<<"sig="<<sigsq[0]<<endl;
      xxinv=(xx.transpose()*xx).inverse()*sigsq[0];
    	cho=(xxinv.llt().matrixL());
    	// cout<<"cho="<<cho<<endl;
      gen.normal(rnor,0,1);
      //VectorXd rbeta(rnor.data());
      beta=xxinv/sigsq[0]*xx.transpose()*yy+cho*rbeta;
      //end change2

      // all_initial[0].t[0].pr();
      //  prxi(all_initial[0].xi);

      jj=0;
      //-------------------------------
      // impute missing values
      for(size_t j=0;j<n;j++)
      {
        for(size_t k=1;k<=nvar;k++) //change3//the row of x used to get prediction
        {
          if(ind_missing[j*(nvar+1)+k]==1)
          {
            for(size_t j1=0;j1<p;j1++)
      	    {
      	      tempx[j1]=(all_initial + 0)->di.x[j*p+j1];
      	    }
      		  //cout<<"missing"<<j<<"and"<<k<<endl;
        		//proposal
            if(vartypeArray[p-k]==0)
          	{
          		prop=gen.normal((all_initial + k)->allfit[j],(all_initial + k)->pi.sigma);
        		  //cout<<"proposal1 = "<<prop<<endl;
        		  //cout<<"mean = "<<(all_initial + k)->allfit[j]<<endl;
        		  //cout<<"sig = "<<(all_initial + k)->pi.sigma<<endl;
          		// calcualte proposal prob
        		  // pro_propdenom=pn(all_initial[k].y[j],all_initial[k].allfit[j],pow(all_initial[k].pi.sigma,2));
        		  pro_propdenom=1.0;
          		//pro_propnum=pn(prop,all_initial[k].allfit[j],pow(all_initial[k].pi.sigma,2));
        		  pro_propnum=1.0;
          	}
            else
            {
              prop=double(abs((all_initial + k)->z[j]-1));
              //cout<<"proposal2"<<prop<<endl;
              pro_propnum=fabs(prop-phi(-(all_initial + k)->allfit[j]));
              pro_propdenom=fabs(double((all_initial + k)->z[j])-phi(-(all_initial + k)->allfit[j]));
    		      //temp_a=phi(-all_initial[k].allfit[j]);
  		        //cout<<"phi"<<temp_a;
  		        //cout<<"zzzzzzzz"<<all_initial[k].z[j]<<"prop"<<prop<<endl;
  		        //cout<<"ppppppp"<<pro_propnum<<"ppppppp"<<pro_propdenom<<endl;
  		      }
      		  // change 4
      		  // cout << "xxtemp1="<<xxtemp1<<endl;
            xxtemp<<1,xxtemp1;
            // cout <<"xxtemp="<<xxtemp<<endl;
    		    meantemp11=(xxtemp.transpose()*beta);
    		    //  cout<<"meantemp11="<<meantemp11<<endl;
    		    // cout<<"beta="<<beta<<endl;
    		    // cout<<"meantemp1="<<meantemp1[0]<<endl;
          	pro_propdenom=pro_propdenom*pn((all_initial + 0)->y[j],meantemp1[0],sigsq[0]);
          	tempx[p-k]=prop;//change the element to proposal
            xxtemp<<1,xxtemp1;
            meantemp21=(xxtemp.transpose()*beta);
          	// cout << "xxtemp1="<<xxtemp1<<endl;
            // cout <<"xxtemp="<<xxtemp<<endl;
            // cout<<"meantemp2="<<meantemp2[0]<<endl;
            pro_propnum=pro_propnum*pn((all_initial + 0)->y[j],meantemp2[0],sigsq[0]);

            l=1; //change 5
            while(l<k)
    		    {
              if(vartypeArray[p-l]==0)
                {
    		          pro_propdenom=pro_propdenom*pn((all_initial + l)->y[j],(all_initial + l)->allfit[j],pow((all_initial + l)->pi.sigma,2));
                }
              else
                {
  			          pro_propdenom=pro_propdenom*fabs(double((all_initial + l)->z[j])-phi(-(all_initial + l)->allfit[j]));
            			// temp_a=phi(-all_initial[l].allfit[j]);
            			//   cout<<"phi"<<temp_a;
                }
              // cout<<"a"<<all_initial[l].y[j]<<"b"<<all_initial[l].allfit[j]<<"C"<<all_initial[l].pi.sigma<<endl;
              di_temp.p=p-l; di_temp.x=&tempx[0];
              di_temp.vartype=vartypeArray;
    		      ppredmean[l]=0.0;
    		      for(size_t q=0;q<m;q++)
        			{
        			  fit((all_initial + l)->t[q],(all_initial + l)->xi,di_temp,fpredtemp);
        			  ppredmean[l] += fpredtemp[0];
        			}
              if(vartypeArray[p-l]==0)
              {
    		        pro_propnum=pro_propnum*pn((all_initial + l)->y[j],ppredmean[l],pow((all_initial + l)->pi.sigma,2));
              }
              else
              {
  		          pro_propnum=pro_propnum*fabs(double((all_initial + l)->z[j])-phi(-ppredmean[l]));
              }
  		        // cout<<"d"<<ppredmean[l]<<endl;
  		        l++;
  		        // cout<<"denom"<<pro_propdenom<<"num"<<pro_propnum;
  		      }


            if(vartypeArray[p-k]==0)
            {
            	pro_prop=std::min(1.0,pro_propnum/pro_propdenom);
          	}
            else
            {
              pro_prop=pro_propnum/(pro_propnum+pro_propdenom);
      		    // if(j==6&&k==19){
      		    // cout<<"proposal_prob"<<pro_prop<<"prop"<<prop<<endl;}
          	}
      		  // cout<<"proposal_prob"<<pro_prop<<endl;
        		//if proposal accepted, update variables

            if(gen.uniform()<pro_prop)
            {if(vartypeArray[p-k]==0){(all_initial + k)->y[j]=prop;}
             else                    {(all_initial + k)->z[j]=int(prop);
                              			  if((all_initial + k)->z[j]==0)
                              			  {(all_initial + k)->y[j]=0.0;
                              			    while((all_initial + k)->y[j]>=0.0)
                              			    {(all_initial + k)->y[j]=gen.normal((all_initial + k)->allfit[j],1);}
                              			  }
                              			  else
                              			  {(all_initial + k)->y[j]=0.0;
                              			    while((all_initial + k)->y[j]<=0.0)
                              			    {(all_initial + k)->y[j]=gen.normal((all_initial + k)->allfit[j],1);}
                              			  }
                }
                l=0;
          		  while(l<k){(all_initial + l)->di.x[j*(p-l)-k+p]=prop;
        		               (all_initial + l)->allfit[j]=ppredmean[l];
                           l++;}
          		  mif<<prop<<endl; //cout<<"mif="<<prop<<endl;// miftemp<<prop<<endl;//   accpt<<1<<endl;
          		  }
            else
            {
              if(vartypeArray[p-k]==0)
                {mif<<all_initial[k].y[j] <<endl;//cout<<"mif1="<<all_initial[k].y[j]<<endl;// miftemp<<all_initial[k].y[j]<<endl;
                }
              else
                {mif<<all_initial[k].z[j] <<endl;//	cout<<"mif2="<<all_initial[k].z[j]<<endl; //miftemp<<all_initial[k].z[j]<<endl;
          		  }
                //   accpt<<0<<endl;
             }
          }
      }
    }
  }
  tp=clock()-tp;
  double thetime = (double)(tp)/(double)(CLOCKS_PER_SEC);
  cout << "time for loop: " << thetime << endl;
  cout << "Value of 'type' was =  "<< type <<endl;
  // std::ofstream timef("time.txt");
  // timef << thetime << endl;
  // cout << "deleting the Array" << endl;
  // delete gen;
   delete[] all_initial;
   return 0;
}
