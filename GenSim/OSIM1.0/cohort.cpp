//
//  cohort.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/19/14.
//  Copyright (c) 2014 Hao. All rights reserved.
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
        cout<<"DONE with individual "<< k <<endl;
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








