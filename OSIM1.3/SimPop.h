//
//  SimPop.h
//  OSIM1.01
//
//  Created by Hao Cheng on 5/11/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef OSIM1_01_SimPop_h
#define OSIM1_01_SimPop_h


#include "cohort.h"
#include "ped.h"


class SimPop {
public:

    cohort founders;
    cohort parents;
    cohort children;
    int myGen;
    
    vector<pedLine> pedTable;

    
    SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate);
  
    //GENERATE FOUNDERS
    void popFounders(unsigned size, string file, unsigned founder_size);
    
    //RANDOMLY MATING FOR nGen GENERATIONS
    void popSample(unsigned size, int nGen);
    //INPUT PEDIGREE
    void inputPedfile(string pedfile);
    //MATING AS PEDIGREE FOR ONE GENERATION
    void pedSample(string pedfile);
    //GET GENOTYPES
    MatrixXf getGenotypes();
    
    //merge with other populations
    void merge(SimPop, unsigned size);

};



#endif
