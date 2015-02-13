//
//  simPop.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include "simPop.h"

SimPop::SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate){
    
    AnimalClass::randGen.seed(314);
    AnimalClass::G.num_breeds=1;
    AnimalClass::G.set_num_chrom(numChr);
    AnimalClass::G.mutRate=mutRate;
    
    unsigned numChromosomePair = AnimalClass::G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
        AnimalClass::G[i].chr_length = chrLength;
    }
    
    
    
    for(auto &chromosome : AnimalClass::G){//evenly create map positions
            
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


SimPop::SimPop(string paramFile, string mapFile){

    AnimalClass::G.num_breeds=1;
    //read ParaMap to get parameteres for  genome information
    ParmMap parameters;
    parameters.inputParms(paramFile);
    
    vector<double> numChr_v      =   parameters["nChrm"];
    vector<double> nLoci_v       =   parameters["nLoci"];
    vector<double> chrLength_v   =   parameters["chrLength"];
    vector<double> mutRate_v     =   parameters["mutRate"];

    int numChr= int(numChr_v[0]);
    AnimalClass::G.set_num_chrom(numChr);
    double mutRate = mutRate_v[0];
    AnimalClass::G.mutRate=mutRate;
    
    int numLoci;
    double chrLength;
    
    int start=0;
    
    //reaed mapPos file to get map positions of all loci
    ifstream datafile;
    datafile.open(mapFile.c_str());
    if(!datafile) {
        cerr << "Couldn't open data file: " << mapFile << endl;
        exit (-1);
    }

    std::string inputStr;
    vector<string> tokens;
    
    for(auto &chromosome : AnimalClass::G){
        
        numLoci=nLoci_v[start];
        chrLength=chrLength_v[start];
        chromosome.set_num_loci(numLoci);
        chromosome.chr_length = chrLength;
        start++;
        
        double mapPos;
        
        int lineNumber= start;
        
        //for(int i=0; i<lineNumber;i++){
        getline(datafile,inputStr); //The outer loop is not needed; Every time the getlline is run, it goes to nexr line in the file as long as the file has not been cleared
        //}
        

        
        boost::split(tokens, inputStr, boost::is_any_of(" "));
       
        unsigned i=0;
        for(auto &locus : chromosome){
            locus.setNumAllelesNumBreeds(2, AnimalClass::G.num_breeds);
            locus.locusType="Marker";
            locus.allele_freq<<0.5,0.5;
            locus.map_pos=getDouble(tokens[i]);
            i++;
        }
        
    }
    
    if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
    datafile.clear();
    datafile.close();

}


void SimPop::popFounders(unsigned founderSize){
    
    if (founders.size()==0) {
        founders.sampleFounders(founderSize);
    }
}

void SimPop::popFounders(unsigned founderSize,string haplotypeGile){
    
    if (founders.size()==0) {
            founders.sampleFounders(founderSize,haplotypeGile);
        }
}



void SimPop::popSample(unsigned size, int nGen){
    
    //this just run once after founders are created
    myGen=1;//need double check
    
    if(founders.size()!= 0){
        nGen--;
        cout << "Generation 1 ---> ";
        parents.sampleChildren(size,founders,founders);
        myGen = 2;
        nGen--;
        founders.clear();
    }
    
    for (int i=1; i<=nGen; i++){
        cout << "Generation " << myGen << " ---> ";
        children.sampleChildren(size,parents,parents);
        //children.showIds();
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


void SimPop::cross(SimPop &a,SimPop &b,unsigned size){
    children.sampleChildren(size,a.parents,b.parents);
    parents.copy(children);
}


void SimPop::sub(SimPop &a,unsigned size){
    children.sampleChildren(size,a.parents,a.parents);
    parents.copy(children);
}

MatrixXf SimPop::getGenotypes(){
    parents.getHaps();
    parents.getNPMatrix();
    return parents.NPMatrix;
}
