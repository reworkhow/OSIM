//
//  genome_info.h
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef __OSIM__genome_info__
#define __OSIM__genome_info__

#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include <vector>
#include <map>

using namespace Eigen;
using namespace std;
using namespace boost;



class Locus_info {
public:
    string locusType;       //QTLs or Markers
    float map_pos;
    MatrixXd allele_freq;   //a row vector of allele for each breed,
                            //thus a matrix for several breeds.
    unsigned get_num_alleles(void){return allele_freq.cols();}
    void setNumAllelesNumBreeds(unsigned numAlleles, unsigned numBreeds);
};

class Chromosome_info: public  vector <Locus_info> {
public:
    void     set_num_loci(unsigned n){resize(n);}
    unsigned get_num_loci(void){return size();}
    long double chr_length;
    ArrayXf MapPos;
    void     mkMapPosFromLocus_info();
};


class Genome_info: public vector<Chromosome_info> {
public:
    Genome_info(){mapPosDone=false;}
    unsigned num_traits, num_breeds, num_chrom;
    double mutRate;
    void display(void);
    unsigned getTotalLoci();
    void set_num_chrom(unsigned n){resize(n);}
    unsigned get_num_chrom(void){return size();}
    
    bool mapPosDone;
    void mkMapPos(){
        for( auto &i : *this){
            i.mkMapPosFromLocus_info();
        }
        mapPosDone=true;
    };
};



#endif /* defined(__OSIM__genome_info__) */
