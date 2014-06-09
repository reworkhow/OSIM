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


using namespace Eigen;
using namespace std;

static default_random_engine randGenTool;

VectorXi sample_without_replace(VectorXi &vec , int k);//sample k elements from vec without replacement

double getDouble(std::string& Str) ;

class corrClass{
public:
    double sumA,sumB;
    double sumA2,sumB2,sumAB;
    double varA,varB,covAB,r,alpha,beta,num;
    double meanA,meanB;
    unsigned n;
    Eigen::ArrayXd vecA,vecB;
    
    corrClass(unsigned dim){
        n     = dim;
        sumA  = 0.0;
        sumB  = 0.0;
        sumA2 = 0.0;
        sumB2 = 0.0;
        sumAB = 0.0;
        varA  = 0.0;
        varB  = 0.0;
        covAB = 0.0;
        r     = 0.0;
    }
    
    void initialize(Eigen::VectorXd a, Eigen::VectorXd b){
        vecA  = a.array();
        vecB  = b.array();
        sumA  = vecA.sum();
        sumB  = vecB.sum();
        sumA2 = vecA.square().sum();
        sumB2 = vecB.square().sum();
        sumAB = (vecA*vecB).sum();
    }
    
    void getCorr(void){
        num = double(n);
        meanA = sumA/num;
        meanB = sumB/num;
        varA  = (sumA2-sumA*meanA)/(num-1);
        varB  = (sumB2-sumB*meanB)/(num-1);
        covAB = (sumAB-sumA*meanB)/(num-1);
        r     = covAB/sqrt(varA*varB);
        beta  = covAB/varB;
        alpha = meanA-beta*meanB;
    }
};


#endif /* defined(__OSIM__ComputaionTools__) */
