#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

#ifndef __particles__
#define __particles__ 

using namespace std;

class particles{

public:

    void set_values (vector<double>,vector<double>,vector<double>);
    vector<double> dir_product (vector<double>,vector<double>);
    vector<double> x;
    vector<double> vx;
    vector<double> mx;
    vector<double> num; // need the numbering of the particles in case there are crossings
};


#endif
