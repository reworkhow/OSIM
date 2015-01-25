//
//  cohort.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include "cohort.h"

void cohort::sampleFounders(unsigned numAnimals)
{
    std::cout << "Sampling " << numAnimals << " animals randomly into base population." <<endl;
    
    for (unsigned i=0; i<numAnimals; i++)
    {
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleFounder();
        this->push_back(animal);//push back to cohort vector
        AnimalClass::founders.push_back(animal);
    }
}

void cohort::sampleFounders(unsigned numAnimals, string hapFile)
{
    std::cout << "Sampling "<< numAnimals << " animals with known genotypes into base population." <<endl;
    
    ifstream datafile;
    datafile.open(hapFile);
    if(!datafile)
    {
        cerr << "Couldn't open data file: " << hapFile << endl;
        exit (-1);
    }
    std::string inputStr;
    vector<string> tokens;
    vector<string> tokens2;
    unsigned lineNumber = 0;
    
    while (getline(datafile,inputStr)){
        
        lineNumber++;
        
        if(lineNumber%2==1)
        {
            boost::split(tokens, inputStr, boost::is_any_of(" "));
        }
        else
        {
            boost::split(tokens2, inputStr, boost::is_any_of(" "));
            AnimalClass *animal = new AnimalClass();
            animal->sampleFounder(tokens,tokens2);
            this->push_back(animal);
            AnimalClass::founders.push_back(animal);
        }
    }
    datafile.close();
}

void cohort::sampleChildren(unsigned numAnimals,cohort &fathers, cohort &mothers)
{
    
    std::cout<<"Sampling "<<numAnimals<<" children randomly."<<endl;
    
    for (unsigned i=0; i<numAnimals; i++)
    {
        AnimalClass *father = fathers.getRandomInd();
        AnimalClass *mother = mothers.getRandomInd();
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleNonFounder(*father,*mother);
        this->push_back(animal);//push back to cohort vector
    }
}

void cohort::sampleChildrenWithPedigree(vector<pedLine> pedTable,cohort &fathers, cohort &mothers)
{
    std::cout<<"Sampling "<<pedTable.size()<<" children with pedigree."<<endl;

    for(unsigned i=0;i<pedTable.size();i++)
    {
        unsigned fatherId=pedTable[i].father;
        unsigned motherId=pedTable[i].mother;
        
        unsigned fatherCount=0;
        unsigned motherCount=0;
        
        for(int i=0;i<fathers.size()&&fathers[i]->myId!=fatherId;i++)
        {
            fatherCount++;
        }
        for(int i=0;i<mothers.size()&&mothers[i]->myId!=motherId;i++)
        {
            motherCount++;
        }
        
        AnimalClass *father = fathers[fatherCount];
        AnimalClass *mother = mothers[motherCount];
        AnimalClass *animal = new AnimalClass();//'new' return the address
        animal->sampleNonFounder(*father,*mother);
        this->push_back(animal);//push back to cohort vector
    }
}

AnimalClass* cohort::getRandomInd(void)
{
    uniform_real_distribution<float> u(0,1);
    float  random = u(AnimalClass::randGen);
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
    for(it=begin();it!=end();it++)
    {
        delete(*it); //delete memory pointed by cohort
    }
    clear();//clean cohort
}

void cohort::showIds()
{
    cout<<"Animal IDs are ";
    for(int i=0;i<this->size();i++)
    {
        cout<<"("
            <<(*this)[i]-> myId     <<", "
            <<(*this)[i]-> sireId   <<", "
            <<(*this)[i]-> damId    <<")\t";
    }
    cout<<endl;
}


void cohort::getHaps()
{
    if(!AnimalClass::G.mapPosDone) AnimalClass::G.mkMapPos();
    for(auto &i: *this)
    {
        i->getMyHaps();
    }
}

void cohort::getNPMatrix()
{
    unsigned nrows= this->size();
    unsigned nclumns=AnimalClass::G.getTotalLoci();
    NPMatrix.resize(nrows,nclumns);
    
    unsigned j=0;
    for(auto &i: *this)
    {
        i->getMyGenotype();
        NPMatrix.row(j)=i->myGenotype;
        j++;
    }
}

void cohort::copy(cohort &from){
    this->flush();
    for (auto i : from)
    {
        this->push_back(i);
    }
    from.clear();//GOOD, will not remove the memory for content in animals class, just remove those animla class
}

unsigned cohort::displaySumNumPos()
{
    unsigned sumPosSize=0;
    for(int i=0;i<this->size();i++)
    {
        sumPosSize+=(*this)[i]-> displayNumPos();
    }
    return(sumPosSize);
}


