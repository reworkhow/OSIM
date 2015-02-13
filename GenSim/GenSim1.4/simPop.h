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
#include <Eigen/Dense>
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
    
    //initialization with constant numLoci, chrLength and random map postions
    // or reading  parameter file and map position file
    SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate);
    SimPop(string parmFile);
    SimPop(string paramFile, string mapFile);
    
    //generate founders
    void popFounders(unsigned founderSize);
    void popFounders(unsigned founderSize, string haplotypeFile);
        
    //randomly mating for nGen generations
    void popSample(unsigned size, int nGen);
    
    //input pedigree
    void inputPedfile(string pedfile);
    
    //mating as pedigree for one generarion
    void pedSample(string pedfile);
    
    //get genotype
    MatrixXf getGenotypes();
    
    //get phenotype
    VectorXf getPhenotype(unsigned numQTL,double heritability);
    
    //crossing from other two lines
    void cross(SimPop &a,SimPop &b,unsigned size);
    
    //get a subset of current population
    void subset(unsigned size);
    
    //initialize genome infomation done or not
    //(genome information just need to be initialized once)
    static bool genomeInfoDone;
    void makeGenomeInfoDone()
    {
        genomeInfoDone=true;
    }
    
    //print out animal IDs fot current generation
    void printAnimalIDs();
};


#endif /* defined(__GenSim1_3__simPop__) */
