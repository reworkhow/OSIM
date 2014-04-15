//
//  genome_info.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include "genome_info.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Dense>
#include "genome_info.h"

//unsigned Genome_info::num_traits;
//unsigned Genome_info::num_breeds;
//unsigned Genome_info::num_chrom;
//unsigned Genome_info::nintval;
//unsigned Genome_info::nMarker;
//unsigned Genome_info::nQTL;
//double   Genome_info::pCrossOver;
//double   Genome_info::pMarkerMut;
//double   Genome_info::pQTLMut;
//unsigned Genome_info::nLoci;
//vector <unsigned> Genome_info::markerDensity;

void Locus_info::setNumAllelesNumBreeds(unsigned numAlleles, unsigned numBreeds)
{
	allele_freq.resize(numBreeds,numAlleles);
}

void Chromosome_info::mkMapPosFromLocus_info(){
    MapPos.resize(get_num_loci());
    for(unsigned i=0;i<MapPos.size();i++){
        MapPos[i]=(*this)[i].map_pos;
    }
}

unsigned Genome_info::getTotalLoci(){
    unsigned totalLoci=0;

    for(auto &i: *this   ){
        totalLoci += i.get_num_loci();
    }
    
    return(totalLoci);
}











//
//void Genome_info::display(void){
//    unsigned n = get_num_chrom();
//    for (unsigned i=0; i < n; i++){
//        cout << "Info for chromosome: " << i+1 << endl;
//        unsigned n_loci = (*this)[i].get_num_loci();
//        for (unsigned j=0; j < n_loci; j++){
//            cout << "Info for locus: " << j+1 << endl;
//            cout << "Allele Frequencies" << endl;
//            unsigned num_alleles = (*this)[i][j].get_num_alleles();
//            for (unsigned br = 0; br < Genome_info::num_breeds; br++){
//                cout <<"Breed: " << br+1 << " ";
//                for (unsigned k=0; k < num_alleles; k++){
//                    cout << setw(7) << setprecision(5) << (*this)[i][j].allele_freq(br,k);
//                }
//                cout << endl;
//            }
//        }
//    }
//    cout << "------------------------------------------------------" << endl;
//}
//
//
//
//void Gamete::setSize(Genome_info G){
//	std::cout<<"Resizing Gametes ";
//	unsigned numberChr = G.get_num_chrom();
//	resize(numberChr);
//	for (unsigned i=0;i<numberChr;i++){
//		unsigned numberLoci = G[i].get_num_loci();
//		(*this)[i].resize(numberLoci);
//	}
//	std::cout<<"  -->DONE."<<endl;
//}

