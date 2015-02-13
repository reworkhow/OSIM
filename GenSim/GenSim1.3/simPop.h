//
//  simPop.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef __GenSim1_3__simPop__
#define __GenSim1_3__simPop__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cohort.h"
#include "ped.h"
#include "parmMap.h"


class SimPop {
public:
    
    cohort founders;
    cohort parents;
    cohort children;
    int myGen;
    
    vector<pedLine> pedTable;
    
    
    SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate);         //initialization with constant numLoci, chrLength and random map postions
    SimPop(string paramFile, string mapFile);                                           //reading  parameter file and map position file
    
    //GENERATE FOUNDERS
    void popFounders(unsigned founderSize);
    void popFounders(unsigned founderSize, string haplotypeFile);
        
    //RANDOMLY MATING FOR nGen GENERATIONS
    void popSample(unsigned size, int nGen);
    //INPUT PEDIGREE
    void inputPedfile(string pedfile);
    //MATING AS PEDIGREE FOR ONE GENERATION
    void pedSample(string pedfile);
    //GET GENOTYPES
    MatrixXf getGenotypes();
    
    //crossing from other two lines
    void cross(SimPop &a,SimPop &b,unsigned size);
    
    //get a sub part of another population
    void sub(SimPop &a,unsigned size);
    
    
};


#endif /* defined(__GenSim1_3__simPop__) */
