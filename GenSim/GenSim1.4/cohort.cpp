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

void cohort::sampleFounders(unsigned numAnimals, string hapFile) //also get allele frequencies
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
    ArrayXf genFreq;
    
    
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
    
    
    if(lineNumber/2!=numAnimals)
    {
        cerr << "the haplotype file for founder animals is wrong"<<endl;
        exit(-1);
    }
    
    //calculate allele frequencies for all loci 
    unsigned numChromosomePair = AnimalClass::G.get_num_chrom();
    unsigned totalLoci = AnimalClass::G.getTotalLoci();
    AnimalClass::G.alleleFreqGenome.resize(totalLoci);
    double allelFreqTemp;
    unsigned nLoci=0;
    
    for(unsigned i=0;i<numChromosomePair;i++)
    {
        unsigned numLoci= AnimalClass::G[i].get_num_loci();
        for(unsigned j=0;j<numLoci;j++){
            allelFreqTemp = AnimalClass::G[i][j].alleleFreq/(numAnimals*2);
            AnimalClass::G[i][j].alleleFreq = allelFreqTemp;
            AnimalClass::G.alleleFreqGenome[nLoci] = allelFreqTemp;
            nLoci++;
        }
    }

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

void cohort::getNPMatrixQTL(unsigned numQTL)
{
    unsigned nrows= this->size();
    unsigned nclumns=numQTL;
    NPMatrixQTL.resize(nrows,nclumns);

    unsigned totalLoci = AnimalClass::G.getTotalLoci();
    uniform_real_distribution<float> u(0,1);
    
    for (unsigned i=0; i<numQTL; i++) {
        
        float  random = u(AnimalClass::randGen);
        unsigned long nQTL = (unsigned long)(random*totalLoci);
        if(nQTL==totalLoci){nQTL=totalLoci-1;};
        if (AnimalClass::G.alleleFreqGenome[nQTL]==0.0) {
            i--;
        }
        
        NPMatrixQTL.col(i)=NPMatrix.col(nQTL);
    }
    
}

void cohort::getPhenotype(unsigned numQTL,double heritability)
{
    unsigned nInd= this->size();

    breedingValue.resize(nInd);
    phenotype.resize(nInd);
    
    //simulate QTL effects
    VectorXf qtlEffect;
    qtlEffect.resize(numQTL);

    normal_distribution<> norm1(0,1);
    for (unsigned i=0; i<numQTL; i++) {
        float  random = norm1(AnimalClass::randGen);
        qtlEffect(i)=random;
    }
    
    //get npQTL matrix
    getNPMatrixQTL(numQTL);
    
    //get breeding value
    breedingValue.resize(nInd);
    breedingValue=NPMatrixQTL*qtlEffect;
    
    //get residual
    corrClass cor(nInd);
    float geneticVariance = cor.getVariance(breedingValue);
    float residualVariance = geneticVariance/heritability-geneticVariance;
    
    normal_distribution<> norm2(0,sqrt(residualVariance));
    VectorXf residual(nInd);
    for (unsigned j=0; j<nInd; j++)
    {
        residual(j) = norm2(AnimalClass::randGen);
    }

    //get phenotype
    phenotype = breedingValue + residual;
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
        //cout << (*this)[i]-> displayNumPos() <<endl;
    }
    
    cout<<"Average of length of ori vector is "<<sumPosSize/(this->size())<<endl;
    return(sumPosSize);
}


