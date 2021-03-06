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
#include "ComputaionTools.h"

using namespace Eigen;
using namespace std;
using namespace boost;


struct chromosome{
    ArrayXf    haplotype;
    ArrayXi    ori;
    ArrayXf    pos;
    ArrayXf    mut;
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
    
    //for founders
    void sampleFounder();
    void sampleFounder(vector<string> tokens1,vector<string> tokens2);
    void initFounderPosOriMut();
    void initFounderHaps();                  //random sample
    void inputFounderHaps(string filename);  //read from a haplotype file
    void inputFounderPatHaps(vector<string> tokens);
    void inputFounderMatHaps(vector<string> tokens);

	//for non-founders
    void sampleNonFounder(AnimalClass& father, AnimalClass& mother);
    void sampleMyPosOriMut(AnimalClass& father, AnimalClass& mother);
    void getMyHaps();
    
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



#endif /* defined(__OSIM__animal_class__) */
