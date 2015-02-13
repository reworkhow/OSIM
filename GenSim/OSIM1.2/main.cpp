//
//  main.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "cohort.h"
#include "ComputaionTools.h"
#include "SimPop.h"
#include "global.h"


int main(int argc, const char * argv[])
{
    if(argc !=1 ){
        unsigned nChrm   =  atoi(argv[1]);
        unsigned nLoci   =  atoi(argv[2]);
        unsigned popSize =  atoi(argv[3]);
        unsigned nGen    =  atoi(argv[4]);
        double   mutRate =  atof(argv[5]);
        cout<<"I've got the parameters you type in. Thank you."<<endl<<endl;
    }
    
    unsigned nLoci   =  5100;
    unsigned nChrm   =     2;
    unsigned popSize =  1000;
    unsigned nGen    =     100;
    double   mutRate =   1e-5;

    
    SimPop osim(nChrm,nLoci,1.0,mutRate);
    osim.popSample(popSize,nGen,"",0);

    return 0;
}
