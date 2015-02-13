//
//  main.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "cohort.h"
#include "tools.h"
#include "simPop.h"
#include "global.h"
#include "parmMap.h"

int main(int argc, const char * argv[])
{
    
    ///user-defined parameters and map positions
    unsigned popSize =  10;
    unsigned nGen    =  3;
    unsigned founderSize=10;
  
    string genomeFile="/Users/erxingfangshui/Dropbox/GenSim/GenSim1.3_xcode/genomeInfo.txt";
    string mapFile="/Users/erxingfangshui/Dropbox/GenSim/GenSim1.3_xcode/mapPos.txt";
    string haplotype="/Users/erxingfangshui/Dropbox/GenSim/GenSim1.3_xcode/haplotype.txt";
    string pedigree="/Users/erxingfangshui/Dropbox/GenSim/GenSim1.3_xcode/pedigree.txt";
    
    SimPop osim(genomeFile,mapFile);
    osim.popFounders(founderSize,haplotype);
   
    SimPop osim1(genomeFile,mapFile);
    SimPop osim2(genomeFile,mapFile);
    SimPop osim3(genomeFile,mapFile);
    
    //osim1.popFounders(popSize,haplotype);
    //osim1.popSample(popSize,nGen);
    
    ///osim2.sub(osim1,popSize);
    ///osim3.sub(osim1,popSize);
    
    ///osim2.popSample(popSize,nGen);
    ///osim3.popSample(popSize,nGen);
    
    //osim.parents.copy(osim.founders);
    //osim.pedSample(pedigree);
    osim.popSample(popSize,nGen);
    //osim.subset(5);
    //osim.popSample(20,2);
    

    //cout<<osim.parents.phenotype<<endl;
    //osim.parents.showIds();

    
    MatrixXf outGeno;
    outGeno=osim.getGenotypes();
    ofstream outFile("/Users/erxingfangshui/Desktop/genotype3.txt");
    outFile << outGeno;
    cout<<"All individuals have been genotyped."<<endl;
    cout<<endl<<endl;
    
    
    VectorXf outPhen=osim.getPhenotype(10, 0.5);
    cout<<outPhen<<endl;
 
    /*
    ///constant nLoci, chrLength and random map positions
    unsigned nLoci   =  1;
    unsigned nChrm   =  1;
    double chrLength =  1;
    unsigned popSize =  1000;
    unsigned nGen    =  100;
    double   mutRate =  0;
    
    SimPop osim(nChrm,nLoci,chrLength,mutRate);
    osim.popFounders(popSize);
    */
    /*
    for (int i=0;i<10;i++ )
    {
        osim.popSample(popSize, nGen);
        osim.parents.displaySumNumPos();
    }
    */
 

    /*
    MatrixXf out;
    out=osim.getGenotypes();
    ofstream outFile("/Users/erxingfangshui/genotype.example");
    outFile << out;
    cout<<"DONE"<<endl;
    return 0;
    */
    
    
    ///ACKNOWLEDGEMENT
    ///printf(R"EOF( /////////////////////////////////////////
    ///                                                 //////
    ///                                                 //////
    ///   _____             ____    ____                //////
    ///  / ____|           / ____|_   _|                //////
    /// | |  __  ___ _ __ | (___   | |  _ __ ___        //////
    /// | | |_ |/ _ \ '_ \ \___ \  | | | '_ ` _ \       //////
    /// | |__| |  __/ | | |____) |_| |_| | | | | |      //////
    ///  \_____|\___|_| |_|_____/|_____|_| |_| |_|      //////
    ///                                                 //////
    ///       )EOF");                                   //////
    ///                                                 //////
    ///cout<<"Thank you for using GenSim1.3. "<<endl;   //////
    //////////////////////////////////////////////////////////
    cout<<endl<<endl<<endl<<endl;
    cout<< "   _____             ____    ___"<<endl;
    cout<< "  / ____|           / ____|_   _|"<<endl;
    cout<< " | |  __  ___ _ __ | (___   | |  _ __ ___"<<endl;
    cout<< " | | |_ |/ _ \\ '_ \\ \\___ \\  | | | '_ ` _ \\ " <<endl;
    cout<< " | |__| |  __/ | | |____) |_| |_| | | | | |"<<endl;
    cout<< "  \\_____|\\___|_| |_|_____/|_____|_| |_| |_|"<<endl;
    cout<< "\n";
    
    cout<<"  Thank you for using GenSim1.3. "<<endl<<endl<<endl;;
    return 0;
}