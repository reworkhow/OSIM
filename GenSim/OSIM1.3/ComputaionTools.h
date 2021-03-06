//
//  ComputaionTools.h
//  OSIM
//
//  Created by Hao Cheng on 4/28/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef __OSIM__ComputaionTools__
#define __OSIM__ComputaionTools__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <Eigen/Dense>
#include <random>
#include <time.h>


using namespace Eigen;
using namespace std;

//SAMPLE WITHOUT REPLACEMENT
static default_random_engine randGenTool;
VectorXi sample_without_replace(VectorXi &vec , int k);//sample k elements from vec without replacement


//TRANSFER STRING TO DOUBLE
double getDouble(std::string& Str) ;


//CALCULATE CORRELATION
class corrClass{
public:
    double sumA,sumB;
    double sumA2,sumB2,sumAB;
    double varA,varB,covAB,r,alpha,beta,num;
    double meanA,meanB;
    unsigned n;
    Eigen::ArrayXd vecA,vecB;
    
    corrClass(unsigned dim);
    void initialize(Eigen::VectorXd a, Eigen::VectorXd b);
    void getCorr(void);
};

//GET CURRENT TIME
string currentDateTime();



#endif /* defined(__OSIM__ComputaionTools__) */
