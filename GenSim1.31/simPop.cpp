//
//  simPop.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include "simPop.h"

SimPop::SimPop(unsigned numChr, unsigned numLoci,double chrLength, double mutRate)
{
    if (genomeInfoDone==false)
    {
        AnimalClass::G.num_breeds=1;
        AnimalClass::G.set_num_chrom(numChr);
        AnimalClass::G.mutRate=mutRate;
    
        unsigned numChromosomePair = AnimalClass::G.get_num_chrom();
    
        for(unsigned i=0;i<numChromosomePair;i++){
            AnimalClass::G[i].chr_length = chrLength;
        }
    
        for(auto &chromosome : AnimalClass::G)//evenly create map positions
        {
            chromosome.set_num_loci(numLoci);
            vector<float> MapPos;
            uniform_real_distribution<float> u(0,1);
            float incr = chromosome.chr_length/numLoci;
            float mappos = incr/10;
        
            for (unsigned k=0; k<chromosome.size(); k++)
            {
                MapPos.push_back(mappos);
                mappos=mappos+incr;
            }
        
            unsigned i=0;
            for(auto &locus : chromosome){
                locus.setNumAllelesNumBreeds(2,AnimalClass::G.num_breeds);
                locus.locusType="Marker";
                locus.allele_freq<<0.5,0.5;
                locus.map_pos=MapPos[i];
                i++;
            }
        }
        
        if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
        //AnimalClass::G.display();
        makeGenomeInfoDone();
    }
}

SimPop::SimPop(string paramFile)
{
    if (genomeInfoDone==false)
    {
        AnimalClass::G.num_breeds=1;
        
        //read ParaMap to get parameteres for  genome information
        ParmMap parameters;
        parameters.inputParms(paramFile);
        
        //set seed
        vector<double> seed_v      =   parameters["seed"];
        int seed= int(seed_v[0]);
        AnimalClass::randGen.seed(seed);

        //set genome infomation
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
        
        for(auto &chromosome : AnimalClass::G)
        {
            
            numLoci=nLoci_v[start];
            chrLength=chrLength_v[start];
            chromosome.set_num_loci(numLoci);
            chromosome.chr_length = chrLength;
            start++;
            
            vector<float> MapPos;
            uniform_real_distribution<float> u(0,1);
            float incr = chromosome.chr_length/numLoci;
            float mappos = incr/10;
            
            for (unsigned k=0; k<chromosome.size(); k++)
            {
                MapPos.push_back(mappos);
                mappos=mappos+incr;
            }
            
            unsigned i=0;
            for(auto &locus : chromosome)
            {
                locus.setNumAllelesNumBreeds(2, AnimalClass::G.num_breeds);
                locus.locusType="Marker";
                locus.allele_freq<<0.5,0.5;
                locus.map_pos=MapPos[i];
                i++;
            }
            
        }
        
        if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
        //AnimalClass::G.display();
        makeGenomeInfoDone();        
    }
    
}


SimPop::SimPop(string paramFile, string mapFile){
    
    if (genomeInfoDone==false)
    {
        AnimalClass::G.num_breeds=1;
        
        //read ParaMap to get parameteres for  genome information
        ParmMap parameters;
        parameters.inputParms(paramFile);
        
        //set seed
        vector<double> seed_v      =   parameters["seed"];
        int seed= int(seed_v[0]);
        cout<< "This is the seed  " <<seed<<endl;

        
        AnimalClass::randGen.seed(seed);

        //set genome infomation
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
    
        //read mapPos file to get map positions of all loci
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
        
            getline(datafile,inputStr);//read one line and go to next
        
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
        //AnimalClass::G.display();
        makeGenomeInfoDone();
    
        datafile.clear();
        datafile.close();
    }
}

void SimPop::popFounders(unsigned founderSize){
    founders.sampleFounders(founderSize);
    myGen=1;
}

void SimPop::popFounders(unsigned founderSize,string haplotypeFile){
    founders.sampleFounders(founderSize,haplotypeFile);
    myGen=1;
}

void SimPop::popSample(unsigned size, int nGen){
        
    if(founders.size()!= 0){
        cout << "Generation 1 ---> ";
        parents.sampleChildren(size,founders,founders);
        myGen = 2;
        founders.clear();
  
        nGen--;
        //printAnimalIDs();
    }
    
    for (int i=1; i<=nGen; i++){
        cout << "Generation " << myGen << " ---> ";
        children.sampleChildren(size,parents,parents);
        parents.copy(children);
        myGen++;
        
        //printAnimalIDs();
    }
}

void SimPop::inputPedfile(string pedfile)
{
    ifstream pedFile(pedfile);
    if(!pedFile)
    {
        cout << "Cannot open " << pedfile << endl;
        exit(-1);
    }
    unsigned individual, father, mother;
    while(pedFile>>individual>>father>>mother)
    {
        pedLine ped;
        ped.individual= individual;
        ped.father    = father;
        ped.mother    = mother;
        pedTable.push_back(ped);
    }
}

void SimPop::pedSample(string pedfile)
{
    inputPedfile(pedfile);

    if(founders.size()!= 0)
    {
        cout << "Generation " << myGen << " ---> ";
        parents.sampleChildrenWithPedigree(pedTable,founders,founders);
        founders.clear();
    }
    else
    {
        cout << "Generation " << myGen << " ---> ";
        children.sampleChildrenWithPedigree(pedTable,parents,parents);
        parents.copy(children);
    }
    
    myGen++;
    //printAnimalIDs();
}

void SimPop::cross(SimPop &a,SimPop &b,unsigned size)
{
    children.sampleChildren(size,a.parents,b.parents);
    parents.copy(children);
    cout << "Generation " << myGen << " ---> "
         << "cross two populations";
    myGen++;
    //printAnimalIDs();
}

void SimPop::subset(unsigned size)
{
    myGen--;
    //cohort::iterator it;
    for(auto it=parents.begin()+size;it!=parents.end();it++)
    {
        delete(*it); //delete memory pointed by cohort
    }

    parents.resize(size);

    cout<< "Generation " << myGen << " ---> "
        << "get subsets of current populations "
        << "-->"<<size << " individuals"<<endl;;
    myGen++;
    //printAnimalIDs();
}

MatrixXf SimPop::getGenotypes()
{
    parents.getHaps();
    parents.getNPMatrix();
    return parents.NPMatrix;
}

void SimPop::printAnimalIDs()
{
    parents.showIds();
}
