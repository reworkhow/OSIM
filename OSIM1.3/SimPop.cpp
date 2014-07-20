//
//  SimPop.cpp
//  OSIM1.2
//
//  Created by Hao Cheng on 7/18/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include "SimPop.h"


SimPop::SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate){
    
    AnimalClass::randGen.seed(314);
    AnimalClass::G.num_breeds=1;
    AnimalClass::G.set_num_chrom(numChr);
    AnimalClass::G.mutRate=mutRate;
    
    unsigned numChromosomePair = AnimalClass::G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
        AnimalClass::G[i].chr_length = chrLength;
    }
    
    for(auto &chromosome : AnimalClass::G){//randomly create map positions
        
        chromosome.set_num_loci(numLoci);
        
        vector<float> MapPos;
        uniform_real_distribution<float> u(0,1);
        float incr = chromosome.chr_length/numLoci;
        float mappos = incr/10;
        
        for (unsigned k=0; k<chromosome.size(); k++) {
            MapPos.push_back(mappos);
            mappos=mappos+incr;
        }
        
        unsigned i=0;
        for(auto &locus : chromosome){
            locus.setNumAllelesNumBreeds(2, AnimalClass::G.num_breeds);
            locus.locusType="Marker";
            locus.allele_freq<<0.5,0.5;
            locus.map_pos=MapPos[i];
            i++;
        }
    }
    
    if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
}

void SimPop::popFounders(unsigned size, string file, unsigned founder_size){
    
    if (founders.size()==0) {
        if(file!=""){
            founders.sampleFounders(founder_size,file);
        }else{
            founders.sampleFounders(size);
        }
    }
}

void SimPop::popSample(unsigned size, int nGen){
    
        nGen--;
        parents.sampleChildren(size,founders,founders);
        myGen = 2;
        nGen--;
    
    for (int i=1; i<=nGen; i++){
        cout << "Generation " << myGen << " ";
        children.sampleChildren(size,parents,parents);
        parents.copy(children);
        myGen++;
    }
}


void SimPop::inputPedfile(string pedfile){
    
    ifstream pedFile(pedfile);
    if(!pedFile){
        cout<<"Cannot open "<<pedfile<<endl;
        exit(-1);
    }
    
    unsigned individual, father, mother;
    while(pedFile>>individual>>father>>mother){
        pedLine ped;
        ped.individual= individual;
        ped.father    = father;
        ped.mother    = mother;
        pedTable.push_back(ped);
    }
}

void SimPop::pedSample(string pedfile){
    
    inputPedfile(pedfile);
    children.sampleChildren_pedigree(pedTable,parents,parents);
    parents.copy(children);
    myGen++;
}


void SimPop::merge(SimPop, unsigned size){
    
    children.sampleChildren(size,parents,SimPop::parents);
    parents.copy(children);

}

MatrixXf SimPop::getGenotypes(){
    parents.getHaps();
    parents.getNPMatrix();
    return parents.NPMatrix;
}




///Reading in mappos
//    for(auto &chromosome : AnimalClass::G){
//
//        unsigned numLoci=10;
//
//        chromosome.set_num_loci(numLoci);
//
//        vector<float> MapPos;
//        float mappos = 0;
//
//        for (unsigned k=0; k<chromosome.size(); k++) {
//            //MapPos.push_back(chromosome.chr_length*u(AnimalClass::randGen));
//
//            ifstream datafile;
//            datafile.open(file.c_str());
//            if(!datafile) {
//                cerr << "Couldn't open data file: " << file << endl;
//                exit (-1);
//            }
//
//            std::string inputStr;
//
//            unsigned i=0;
//            while (getline(datafile,inputStr)){
//                mappos = getDouble(inputStr);
//                MapPos.push_back(mappos);
//            }
//
//            datafile.clear();
//            datafile.close();
//
//        }
//
//        sort(MapPos.begin(),MapPos.end());
//
//        unsigned i=0;
//        for(auto &locus : chromosome){
//            locus.setNumAllelesNumBreeds(2, AnimalClass::G.num_breeds);
//            locus.locusType="Marker";
//            locus.allele_freq<<0.5,0.5;
//            locus.map_pos=MapPos[i];
//            i++;
//       }
//    }

