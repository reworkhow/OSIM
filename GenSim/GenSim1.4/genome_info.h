//
//  genome_info.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef __GenSim1_3__genome_info__
#define __GenSim1_3__genome_info__

#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include <vector>
#include <map>

using namespace Eigen;
using namespace std;
using namespace boost;

class Locus_info
{
public:
    string locusType;       //QTLs or Markers
    float map_pos;
    //MatrixXd allele_freq;   //a row vector of alleles for each breed,thus a matrix for several breeds.
    double alleleFreq;
};

class Chromosome_info: public vector <Locus_info>
{
public:
    void     set_num_loci(unsigned n){resize(n);}
    unsigned get_num_loci(void){return size();}
    long double chr_length;
    ArrayXf MapPos;
    void     mkMapPosFromLocus_info();
};

class Genome_info: public vector<Chromosome_info>
{
public:
    Genome_info(){mapPosDone=false;}
    unsigned num_traits, num_breeds, num_chrom;
    double mutRate;
    void display(void);
    unsigned getTotalLoci();
    void set_num_chrom(unsigned n){resize(n);}
    unsigned get_num_chrom(void){return size();}
    
    bool mapPosDone;
    void mkMapPos()
    {
        for( auto &i : *this){
            i.mkMapPosFromLocus_info();
        }
        mapPosDone=true;
    };
    
    ArrayXd alleleFreqGenome;
};


#endif /* defined(__GenSim1_3__genome_info__) */
