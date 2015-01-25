//
//  genome_info.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Dense>
#include "genome_info.h"

void Locus_info::setNumAllelesNumBreeds(unsigned numAlleles, unsigned numBreeds)
{
    allele_freq.resize(numBreeds,numAlleles);
}

void Chromosome_info::mkMapPosFromLocus_info()
{
    MapPos.resize(get_num_loci());
    for(unsigned i=0;i<MapPos.size();i++){
        MapPos[i]=(*this)[i].map_pos;
    }
}

unsigned Genome_info::getTotalLoci()
{
    unsigned totalLoci=0;
    for(auto &i: *this   ){
        totalLoci += i.get_num_loci();
    }
    
    return(totalLoci);
}

void Genome_info::display(void)
{
    unsigned n = get_num_chrom();
    for (unsigned i=0; i < n; i++)
    {
        cout << "Info for chromosome: " << i+1 << endl;
        unsigned n_loci = (*this)[i].get_num_loci();
        for (unsigned j=0; j < n_loci; j++)
        {
            cout << "Info for locus: " << j+1 << endl;
            cout << "Allele Frequencies" << endl;
            unsigned num_alleles = (*this)[i][j].get_num_alleles();
            for (unsigned br = 0; br < Genome_info::num_breeds; br++)
            {
                cout <<"Breed: " << br+1 << " ";
                for (unsigned k=0; k < num_alleles; k++)
                {
                    cout << setw(7) << setprecision(5) << (*this)[i][j].allele_freq(br,k);
                }
                cout << endl;
            }
        }
    }
    cout << "------------------------------------------------------" << endl;
}
