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
    
    
    unsigned nLoci   =  1000;
    unsigned nChrm   =     2;
    double chrLength =   0.1;
    unsigned popSize =  1000;
    unsigned nGen    =   100;
    double   mutRate =  1e-5;

    
    SimPop osim(nChrm,nLoci,chrLength,mutRate);
    SimPop osim2(nChrm,nLoci,chrLength,mutRate);
    osim.popFounders(popSize,"",0);
    osim2.popFounders(popSize,"",0);


    osim.popSample(popSize,nGen);
    osim2.popSample(popSize,nGen);
    
    osim.merge(osim2,popSize);
   
    osim.popSample(popSize/2,nGen);
    
    cout<<AnimalClass::founders.size()<<endl;
    cout<<AnimalClass::mutants.size()<<endl;
    cout<<AnimalClass::countChromosome<<endl;


    
    MatrixXf out;
    out=osim.getGenotypes();
    
    ofstream outFile("/Users/erxingfangshui/Documents/genotype_"+currentDateTime());
    outFile <<out;
    cout<<"DONE"<<endl;
    
    
    return 0;
}

