#include "func_initialInput.h"


// using std::cout;
// using std::endl;

void init_input(size_t index, init &ip_initial,
                 int* vartype, size_t& m, double &kfac, int &bistart, int &binum, size_t &n, size_t &p,
                 std::vector<double> &z,std::vector<double> &y, std::vector<double> &x)
{
  std::vector<tree> t(m);
  ip_initial.t=t;
  ip_initial.di.x = new double [n*(p-index)];
  ip_initial.allfit = new double [n];
  //read in data
  //read y NOTE: we assume y is already centered at 0.
  double miny = INFINITY; //use range of y to calibrate prior for bottom node mu's
  double maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
  double ytemp;

  if(index==0)
  {
    if(vartype[p-index]==0)
    {
      //cout<<"continuous"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.y.push_back(y[i]);
      }
    }
    else
    {
     // cout<<"binary"<<endl;
      for(size_t i=0;i<n;i++)
      {
        ip_initial.y.push_back(z[i*binum+p-bistart]);
        ip_initial.z.push_back(int(y[i]));
      }
    }
    //cout<<ip_initial.y[0]<<endl;
    ip_initial.di.x=&x[0];
    //cout<<ip_initial.di.x[0]<<endl;
  }

  else
  {
    if(vartype[p-index]==0)
    {
      //cout<<"continuous"<<endl;
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
      //cout<<"binary"<<endl;
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
  //cout<<"\ndata import ok"<<endl;
  for(size_t i=0;i<n;i++) {
    ytemp=ip_initial.y[i];
    if(ytemp<miny) miny=ytemp;
    if(ytemp>maxy) maxy=ytemp;
    allys.sy += ytemp; // sum of y
    allys.sy2 += ytemp*ytemp; // sum of y^2
  }
  allys.n = n;
  // cout << "\ny read in:\n";
  // cout << "n: " << n << endl;
  // cout << "y first and last:\n";
  // cout << ip_initial.y[0] << ", " << ip_initial.y[n-1] << endl;
  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
  //cout << "ybar,shat: " << ybar << ", " << shat <<  endl;

  //read x
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  // cout << "\nx read in:\n";
  // cout << "p: " << p << endl;
  // cout << "first row: " <<  ip_initial.di.x[0] << " ...  " << ip_initial.di.x[p-1-index] << endl;
  // cout << "last row: " << ip_initial.di.x[(n-1)*(p-index)] << " ...  " << ip_initial.di.x[n*(p-index)-1] << endl;

  //x cutpoints
  size_t nc=100; //100 equally spaced cutpoints from min to max.
  makexinfo(p-index,n,&(ip_initial.di.x[0]),ip_initial.xi,nc,vartype);
  //cout<<ip_initial.di.x[0]<<endl;
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

    // cout << "\nalpha, beta: " << ip_initial.pi.alpha << ", " << ip_initial.pi.beta << endl;
    // cout << "sigma, tau: " << ip_initial.pi.sigma << ", " << ip_initial.pi.tau << endl;

  }
  else
  {
    ip_initial.pi.tau=3.0/(kfac*sqrt((double)m));
    ip_initial.pi.sigma=1.0;
  }

  //cout<<ip_initial.allfit[n-1]<<endl;
  ip_initial.di.n=n;

  //cout<<ip_initial.di.n<<endl;
  ip_initial.di.p=p-index;

  //cout<<ip_initial.di.p<<endl;
  ip_initial.di.y=ip_initial.r; //the y for each draw will be the residual

  //cout<<ip_initial.r<<endl;
  ip_initial.di.vartype=vartype;

  return;
}
