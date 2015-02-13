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
#include "ComputaionTools.h"
#include "SimPop.h"
#include "SimPop_Gus.h"

default_random_engine AnimalClass::randGen;
Genome_info AnimalClass::G;
unsigned AnimalClass::countChromosome=0;
unsigned AnimalClass::countId=0;
vector<AnimalClass*> AnimalClass::founders;
vector<mutantInfo*> AnimalClass::mutants;
//vector<mutantInfo> AnimalClass::mutants;


int main(int argc, const char * argv[])
{
    if(argc !=1 ){
        unsigned nChrm   =  atoi(argv[1]);
        unsigned nLoci   =  atoi(argv[2]);
        unsigned popSize =  atoi(argv[3]);
        unsigned nGen    =  atoi(argv[4]);
        double   mutRate =  atof(argv[5]);
        cout<<"I've got the parameters you type in. Thank you."<<endl<<endl;
    }
    
    unsigned nLoci   =  5100;
    unsigned nChrm   =     1;
    unsigned popSize =  1000;
    unsigned nGen    =     1;
    double   mutRate =   1e-5;

    
//    SimPop osim(nChrm,nLoci,mutRate);
//    osim.sample(popSize,nGen);
    
    SimPop_Gus osim(nChrm,nLoci,mutRate);
    
    //string file= "/Users/erxingfangshui/Dropbox/sharedFolderRohan/OSIM1.01/OutputOSIM_5m.txt";
    string file= "";

    
    ofstream datafile;
    string Varfile= "/Users/erxingfangshui/Dropbox/sharedFolderRohan/OSIM1.01/varFile.txt";
    datafile.open(Varfile);
    
    MatrixXf XPX,XX,XQ;
    unsigned nMarkers=5000;
    unsigned nQTL=nLoci-nMarkers;
    
    XPX.setZero(nLoci,nLoci);
    XQ.setZero(nMarkers,nQTL);
    XX.setZero(nMarkers,nMarkers);

    
    
    //sample loci as markers or QTL
    VectorXi vec;
    vec.resize(nLoci);
    for(int i=0;i<nLoci;i++){
        vec[i]= i;
    }
    VectorXi shuffle,whichQTL,whichMarker;
    shuffle.resize(nLoci);
    shuffle=sample_without_replace(vec,nLoci);//vec is changed after this
    whichQTL   =shuffle.head(nQTL);
    whichMarker=shuffle.tail(nMarkers);
    
    
    //sample QTL effect
    VectorXf effects;
    effects.resize(nQTL);
    
    normal_distribution<float> rnorm(0.0,2.0);
    
    for(unsigned i=0;i<nQTL; i++){
        //effects[i]=rnorm(AnimalClass::randGen);
        effects(i)=1;

    }
    
    MatrixXf out;

    for(unsigned sample=0;sample<1;sample++){
        
        osim.sample_Gus(popSize,nGen,file,1092);
        cout<<"This is the  "<<sample+1<<"  X"<<popSize<<" individuals"<<endl;

        out=osim.getGenotypes();
        
        osim.parents.flush();
        
        VectorXf colMean=out.colwise().mean();
        out.rowwise() -= colMean.transpose();
        
        XPX =XPX+(out.transpose()*out-popSize*XPX)/((sample+1)*popSize);
        //XPX=out.transpose()*out;
        
        //extract Marker matrix Sigma_XX from
        int rowi,columnj;
        
        for(int i=0;i<nMarkers;i++){
            rowi=whichMarker(i);
            for(int j=0;j<nMarkers;j++){
                columnj=whichMarker(j);
                XX(i,j)=XPX(rowi,columnj);
            }
        }
        cout<< sample << " extracting XX is done" <<endl;
        
        //extract QTL matrix Sigma_XQ from XPX matrix
        for(int i=0;i<nMarkers;i++){
            rowi=whichMarker(i);
            for(int j=0;j<nQTL;j++){
                columnj=whichQTL(j);
                XQ(i,j)=XPX(rowi,columnj);
            }
        }
        cout<< sample << " extracting XQ is done" <<endl;
        
        //calculate variance
        double var;
        var=effects.transpose()*XQ.transpose()*XX.partialPivLu().inverse()*XQ*effects;
        cout<<var<<endl;
        var=effects.transpose()*XQ.transpose()*XQ*effects;

        
        cout<<var<<endl;
        datafile<<var<<'\n';
        
    }
    
    datafile.close();
    std::cout << "DONE!\n";

    return 0;
}

