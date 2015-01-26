//
//  tools.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include "tools.h"

double getDouble(std::string& Str) {
    std::istringstream inputStrStream(Str.c_str());
    double val;
    inputStrStream >> val;
    return val;
}

VectorXi sample_without_replace(VectorXi &vec , int k){
    
    uniform_real_distribution<float> u(0,1);
    unsigned size=vec.size();
    int temp;
    
    for(int i=0;i<k;i++){
        float  random = u(randGenTool);
        int length=size-i;
        unsigned which = unsigned (random*length);
        if(which==length){which=which-1;};
        
        temp=vec[i];
        vec[i]=vec[which];
        vec[which]=temp;
    }
    
    return vec.head(k);
};

corrClass::corrClass(unsigned dim){
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

void corrClass::initialize(Eigen::VectorXf a, Eigen::VectorXf b){
    vecA  = a.array();
    vecB  = b.array();
    sumA  = vecA.sum();
    sumB  = vecB.sum();
    sumA2 = vecA.square().sum();
    sumB2 = vecB.square().sum();
    sumAB = (vecA*vecB).sum();
}

void corrClass::getCorr(void){
    num = float(n);
    meanA = sumA/num;
    meanB = sumB/num;
    varA  = (sumA2-sumA*meanA)/(num-1);
    varB  = (sumB2-sumB*meanB)/(num-1);
    covAB = (sumAB-sumA*meanB)/(num-1);
    r     = covAB/sqrt(varA*varB);
    beta  = covAB/varB;
    alpha = meanA-beta*meanB;
}

float corrClass::getVariance(Eigen::VectorXf a)
{
    num = float(n);
    vecA  = a.array();
    sumA  = vecA.sum();
    sumA2 = vecA.square().sum();
    meanA = sumA/num;
    varA  = (sumA2-sumA*meanA)/(num-1);
    return(varA);
}

string currentDateTime()
{
    time_t     now = time(0);
    struct tm  * timeinfo;
    
    time (&now);
    timeinfo = localtime(&now);
    
    return asctime(timeinfo);
}
