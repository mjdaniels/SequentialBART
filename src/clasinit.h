#ifndef CLASSINIT_H
#define CLASSINIT_H

#include "info.h"
#include <vector>
#include "tree.h"


class init { //
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

#endif //CLASSINIT_H
