//
//  ComputaionTools.cpp
//  OSIM
//
//  Created by Hao Cheng on 4/28/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include "ComputaionTools.h"

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

