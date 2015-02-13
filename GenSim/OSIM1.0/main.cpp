//
//  main.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "cohort.h"

default_random_engine AnimalClass::randGen;
Genome_info AnimalClass::G;
unsigned AnimalClass::countChromosome=0;
unsigned AnimalClass::countId=0;
vector<AnimalClass*> AnimalClass::founders;
vector<mutantInfo*> AnimalClass::mutants;
//vector<mutantInfo> AnimalClass::mutants;


int main(int argc, const char * argv[])
{
    AnimalClass::G.num_breeds=1;
    //AnimalClass::G.set_num_chrom(50);
    AnimalClass::G.set_num_chrom(1);

    //AnimalClass::mutants.resize(AnimalClass::G.get_num_chrom());   STUPID
    AnimalClass::G.mutRate=1e-4;
    //AnimalClass::G.mutRate=0;


    
    unsigned numChromosomePair = AnimalClass::G.get_num_chrom();
    for(unsigned i=0;i<numChromosomePair;i++){
    AnimalClass::G[i].chr_length = 1;//with chr_length=1,mapdis is almost = recom rate.
    }
    //unsigned Ne=100;
    unsigned Ne=1500;

    
    for(auto &chromosome : AnimalClass::G){
        
        //unsigned numLoci=10;
        //unsigned numLoci=1000;
        unsigned numLoci=30000;

        chromosome.set_num_loci(numLoci);

        vector<float> MapPos;
        uniform_real_distribution<float> u(0,1);
        float incr = chromosome.chr_length/numLoci;
        float mappos = 0.000001;

        for (unsigned k=0; k<chromosome.size(); k++) {
            //MapPos.push_back(chromosome.chr_length*u(AnimalClass::randGen));
            MapPos.push_back(mappos);
            mappos=mappos+incr;
        }
        sort(MapPos.begin(),MapPos.end());
        
        unsigned i=0;
        for(auto &locus : chromosome){
            locus.setNumAllelesNumBreeds(2, AnimalClass::G.num_breeds);
            locus.locusType="Marker";
            locus.allele_freq<<0.5,0.5;
            locus.map_pos=MapPos[i];
            i++;
       }
    }

    AnimalClass::randGen.seed(314);
    cohort founderBoys,founderGirls;
    cohort boysOdd,girlsOdd;
    cohort boysEvn,girlsEvn;
    
    founderBoys.sampleFounders(Ne);
 
//    for(unsigned locus=0;locus<999;locus++){
//    cout<<(boysEvn[0])->GenomePat[0].haplotype[locus]<<endl;
//    }

    //founderGirls.sampleFounders(Ne);

    boysEvn.sampleChildren(Ne,founderBoys,founderBoys);
    //girlsEvn.sampleChildren(Ne,founderBoys,founderGirls);

    //unsigned nGen=1598;
    unsigned nGen=2998;

    for(unsigned gen=0;gen<nGen;gen++){
        
        if(gen%2){
            cout <<endl<< "This is generation "<< gen <<endl<<endl;
            boysEvn.sampleChildren(Ne,boysOdd,boysOdd);
            //boysEvn.display();
            //girlsEvn.sampleChildren(Ne,boysOdd,boysOdd);
            //girlsEvn.display();

            boysOdd.flush();
            //girlsOdd.flush();
            
            cout<< "Average number of pos in generation " <<gen<<" is ---> " << boysEvn.displaySumNumPos()/(2*Ne) <<endl;

        }else{
            cout << endl<<"This is generation "<< gen <<endl<<endl;
            boysOdd.sampleChildren(Ne,boysEvn,boysEvn);
            //boysOdd.display();
            //girlsOdd.sampleChildren(Ne,girlsEvn,girlsEvn);
            //girlsOdd.display();

            boysEvn.flush();
            //girlsEvn.flush();
            
            cout<< "Average number of pos in generation " <<gen<<" is ---> " << boysOdd.displaySumNumPos()/(2*Ne) <<endl;
        }
    
    }
    
    cout<<AnimalClass::mutants.size()<<endl<<endl;
//    for(int i=0;i<AnimalClass::mutants.size();i++){
//        std::cout << std::fixed << std::setw( 11 ) << std::setprecision(9) <<
//        AnimalClass::mutants[i]->whichStartPos<<"   "<<
//        AnimalClass::mutants[i]->whichEndPos<<"  "<<
//        AnimalClass::mutants[i]->mutPos<<"   "
//        <<((AnimalClass::mutants[i]->whichEndPos)>(AnimalClass::mutants[i]->mutPos))+10000000<<endl;
//    }
    

    boysEvn.getHaps();
    boysEvn.getNPMatrix();
    //cout<<boysEvn.NPMatrix<<endl;
    MatrixXf out= boysEvn.NPMatrix;
    ArrayXf genfreq=out.colwise().mean()/2;

    cout<<endl;
    //cout<<genfreq<<endl;
    std::cout << "DONE!\n";
    
//    ofstream datafile;
//    string genfreqfile= "/Users/erxingfangshui/Dropbox/sharedFolderRohan/OSIM/genfreq_1ch_0.1M_10000SNP_" + to_string(nGen+2) + ".txt";
//    datafile.open(genfreqfile);
//    datafile<<genfreq<<'\n';
//    datafile.close();
    
//    ofstream datafile2;
//    string NPMatrixfile= "/Users/erxingfangshui/Dropbox/sharedFolderRohan/OSIM/NPMatrix_1ch_100M_100000SNP_" + to_string(nGen+2) + ".txt";
//    datafile2.open(NPMatrixfile);
//    datafile2<<out<<'\n';
//    datafile2.close();
    
    
    //system("Rscript /Users/erxingfangshui/Dropbox/sharedFolderRohan/OSIM/plotGenFreq.R");
    return 0;
}

