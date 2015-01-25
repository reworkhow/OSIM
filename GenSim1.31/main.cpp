
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
    unsigned nGen    =  10;
  
    string path="/Users/erxingfangshui/Dropbox/GenSim/GenSim1.4/data/";
    string genomeFile=path+"genomeInfo.txt";
    string mapFile=path+"mapPos.txt";
    string haplotype=path+"haplotype.txt";

    SimPop osim(genomeFile,mapFile);
    osim.popFounders(popSize,haplotype);
    
    osim.popSample(popSize,nGen);
    MatrixXf out;
    out=osim.getGenotypes();
    ofstream outFile(path+"genotype.txt");
    outFile << out;
    
    return 0;

}
