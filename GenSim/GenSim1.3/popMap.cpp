//
//  popMap.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/14/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#include "popMap.h"

void PopMap::inputParms(string fileName){
 
    std::ifstream datafile(fileName);
    if (!datafile) {
        std::cerr << "Couldn't open parmFile: " <<fileName << std::endl;
        exit (-1);
    }
    
    string parmName, inputStr,sep;
    double parmValue;
    vector<double> parmValueVector;
    PopMap::iterator mapit;
    
    while (getline(datafile,inputStr)) {
        
        vector<string> tokens;
        boost::split(tokens, inputStr, boost::is_any_of(" "));
        
        unsigned numArg = tokens.size();
        
        parmName = tokens[0];
        
        int numPar = tokens.size();
        if(numPar==3){
            
        }
        
        
        for (unsigned i=1;i<numArg;i++){
            parmValue = getDouble(tokens[i]);
            parmValueVector.push_back(parmValue);
        }
        
        vector<double> parmV;
        parmV = parmValueVector;
 //       (*this)[parmName] = parmV;
        parmValueVector.clear();
    }
    
    datafile.clear();
    datafile.close();
}

