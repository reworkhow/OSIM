//
//  animal_class.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef __GenSim1_3__animal_class__
#define __GenSim1_3__animal_class__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include "genome_info.h"
#include "tools.h"

using namespace Eigen;
using namespace std;
using namespace boost;


struct chromosome{
    ArrayXf    haplotype;
    ArrayXi    ori;
    ArrayXf    pos;
};

struct mutantInfo{
    unsigned whichOri;
    float    whichStartPos;
    float    whichEndPos;
    float    mutPos;
    ArrayXf  val;
};

class AnimalClass {
    
public:
    AnimalClass(int father, int mother);
    
    AnimalClass(void){
        myId = countId++;
        unsigned numChromosomePair = G.get_num_chrom();
        GenomePat.resize(numChromosomePair);
        GenomeMat.resize(numChromosomePair);
    };
    
    vector<chromosome> GenomePat;
    vector<chromosome> GenomeMat;
    
    int myId, sireId, damId;
    
    static vector<AnimalClass*> founders;
    static vector<mutantInfo*> mutants;//Do I need a pointer here?
    
    //for founders
    void sampleFounder();
    void sampleFounder(vector<string> tokens1,vector<string> tokens2);
    void initFounderPosOri();
    void initFounderHaps();                  //random sample
    void inputFounderHaps(string filename);  //read from a haplotype file
    void inputFounderPatHaps(vector<string> tokens);
    void inputFounderMatHaps(vector<string> tokens);
    
    //for non-founders
    void sampleNonFounder(AnimalClass& father, AnimalClass& mother);
    void sampleMyPosOri(AnimalClass& father, AnimalClass& mother);
    void sampleMyMutation();
    void getMyHaps();
    ArrayXf getMyHapSeg(int i,int myOri,int start,int numcopy);
    
    void display();
    unsigned displayNumPos();
    ArrayXf myGenotype;
    void getMyGenotype();
    
    //static:declare and initialize outside
    static Genome_info G;
    static unsigned countChromosome;
    static unsigned countId;
    static default_random_engine randGen;
};



#endif /* defined(__GenSim1_3__animal_class__) */
