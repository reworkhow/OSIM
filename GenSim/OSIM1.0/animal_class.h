//
//  animal_class.h
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef __OSIM__animal_class__

#define __OSIM__animal_class__

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

using namespace Eigen;
using namespace std;
using namespace boost;


struct chromosome{
//    vector<unsigned> haplotype;
    ArrayXf    haplotype;
    ArrayXi    ori;
    ArrayXf    pos;
};

struct mutantInfo{
    unsigned whichOri;
    float whichStartPos;
    float whichEndPos;
    float mutPos;
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

    static Genome_info G;
    
	int myId, sireId, damId;
    
    static vector<AnimalClass*> founders;
    static vector<mutantInfo*> mutants;//Do I need a pointer here?
    //static vector<mutantInfo> mutants;//Do I need a pointer here?

    
    //int breed;
    void initFounderPosOri(); //for founders
    void initFounderHaps(); //random sample
    void inputFounderHaps(string filename);  //read from a haplotype file
	void sampleMyPosOri(AnimalClass& father, AnimalClass& mother);//for non-founders
    void sampleMyMutation();
    void getMyHaps();
    ArrayXf getMyHapSeg(int i,int myOri,int start,int numcopy);
    
    void sampleFounder();
    void sampleNonFounder(AnimalClass& father, AnimalClass& mother);

    void display();
    unsigned displayNumPos();
    ArrayXf myGenotype;
    void getMyGenotype();
    
    static unsigned countChromosome;//declare and initialize outside
    static unsigned countId;
    static default_random_engine randGen;
};



#endif /* defined(__OSIM__animal_class__) */
