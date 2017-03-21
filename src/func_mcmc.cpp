
#include "func_mcmc.h"

void mcmc(init& ip_initial, RNG gen,int* vartype, size_t m, size_t n, double lambda, double nu)
{

  for(size_t j=0;j<m;j++)
  {

    fit(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.ftemp);

    for(size_t k=0;k<n;k++)
    {
      ip_initial.allfit[k] = ip_initial.allfit[k]-ip_initial.ftemp[k];
      ip_initial.r[k] = ip_initial.y[k]-ip_initial.allfit[k];
    }
    bd(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.pi,gen);

    drmu(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.pi,gen);

    fit(ip_initial.t[j],ip_initial.xi,ip_initial.di,ip_initial.ftemp);

    for(size_t k=0;k<n;k++) ip_initial.allfit[k] += ip_initial.ftemp[k];
  }

  if(vartype[ip_initial.di.p]==0)
  {
    double rss;
    double restemp;
    rss=0.0;
    for(size_t k=0;k<n;k++)
    {
      restemp=ip_initial.y[k]-ip_initial.allfit[k];
      rss += restemp*restemp;
    }
    ip_initial.pi.sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
  } else
  {
    for(size_t k=0;k<n;k++)
    { if(ip_initial.z[k]==0)
      {
      ip_initial.y[k]=0.0;
      while(ip_initial.y[k]>=0.0)
        {
         ip_initial.y[k]=gen.normal(ip_initial.allfit[k],1);
        }
      } else
      {
      ip_initial.y[k]=0.0;
      while(ip_initial.y[k]<=0.0)
        {
        ip_initial.y[k]=gen.normal(ip_initial.allfit[k],1);
        }
      }
    }
  }

  return;
 }
