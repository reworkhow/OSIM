//
//  cohort.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/19/14.
//  Copyright (c) 2014 Hao. All rights reserved.S
//

#include "cohort.h"

void cohort::sampleFounders(unsigned numAnimals)
{
    std::cout<<"Sampling "<<numAnimals<<" animals into base population."<<endl;
    for (unsigned i=0; i<numAnimals; i++){
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleFounder();
        this->push_back(animal);//push back to cohort vector
        AnimalClass::founders.push_back(animal);
    }
}

void cohort::sampleFounders(unsigned numAnimals, string file)
{
    std::cout<<"Sampling "<<numAnimals<<" animals from known genotypes into base population."<<endl;

    ifstream datafile;
    datafile.open(file.c_str());
    if(!datafile) {
    cerr << "Couldn't open data file: " << file << endl;
    exit (-1);
    }
    std::string inputStr;
    vector<string> tokens;
    vector<string> tokens_temp;
    unsigned lineNumber = 0;
    unsigned numCols, n;

    int nMarkers;
    
    while (getline(datafile,inputStr)){
        
        boost::split(tokens, inputStr, boost::is_any_of(" "));
        lineNumber++;
        
        if (lineNumber==1){
            nMarkers = tokens.size();
        }else{
            n = tokens.size();
            if (n != nMarkers) {
                cerr << "Line "  << lineNumber << " of data file "
                << file << "  has "   << n <<" columns; "
                << nMarkers << " expected " << endl;
                exit (-1);
            }
        }
        
        if(lineNumber%2){
            tokens_temp=tokens;
        }else{
            AnimalClass *animal = new AnimalClass(); //'new' return the address
            animal->sampleFounder(tokens_temp,tokens);
            this->push_back(animal);//push back to cohort vector
            AnimalClass::founders.push_back(animal);
        }
        
    }
    
    datafile.close();
    
}



void cohort::sampleChildren(unsigned numAnimals,cohort &fathers, cohort &mothers)
{
    
    std::cout<<"Sampling "<<numAnimals<<" children into next generation."<<endl;
    
    for (unsigned i=0; i<numAnimals; i++){
        AnimalClass *father = fathers.getRandomInd();
        AnimalClass *mother = mothers.getRandomInd();
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleNonFounder(*father,*mother);
        this->push_back(animal);//push back to cohort vector
        
    }

}

void cohort::sampleChildren_pedigree(vector<pedLine> pedTable,cohort &fathers, cohort &mothers){
    
    for(unsigned i=0;i<pedTable.size();i++){
        
        unsigned fatherId=pedTable[i].father;
        unsigned motherId=pedTable[i].mother;
        
        AnimalClass *father = fathers[fatherId];
        AnimalClass *mother = mothers[motherId];
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleNonFounder(*father,*mother);
        this->push_back(animal);//push back to cohort vector

    }
}

AnimalClass* cohort::getRandomInd(void){
    uniform_real_distribution<float> u(0,1);
    float  random = u(AnimalClass::randGen);
	//unsigned size = this->size();
    unsigned long size = this->size();

	unsigned long i = (unsigned long)(random*size);
    if(i==size){i=i-1;};
	return (*this)[i];
}

void cohort::display()
{
    for(int i=0;i<this->size();i++){
        (*this)[i]->display();
    }
}

void cohort::flush()
{
    cohort::iterator it;
    for(it=begin();it!=end();it++){
        delete(*it); //delete memory pointed by cohort
    }
    clear();//clean cohort
}

unsigned cohort::displaySumNumPos()
{
    unsigned sumPosSize=0;

    for(int i=0;i<this->size();i++){
        sumPosSize+=(*this)[i]-> displayNumPos();
    }
    return(sumPosSize);

}

void cohort::getHaps()
{
    if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
    int k=0;
    for(auto &i: *this   ){
        i->getMyHaps();
        //cout<<"DONE with individual "<< k <<endl;
        k++;
    }
}

void cohort::getNPMatrix()
{
    unsigned nrows= this->size();
    unsigned nclumns=AnimalClass::G.getTotalLoci();
    NPMatrix.resize(nrows,nclumns);
    
    unsigned j=0;
    for(auto &i: *this   ){
        i->getMyGenotype();
        NPMatrix.row(j)=i->myGenotype;
        j++;
    }
}

void cohort::copy(cohort &from){
    this->flush();
    for (auto i : from){
        this->push_back(i);
    }
    from.clear();//GOOD, will not remove the memory for content in animals class, just remove those animla class
}



    









