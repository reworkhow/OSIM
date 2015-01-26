//
//  parmMap.cpp
//  GenSim1.3
//
//  Created by Hao Cheng on 12/13/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//  learn parmMap from MATVEC


#include "parmMap.h"

void ParmMap::inputParms(string fileName){
    std::ifstream datafile(fileName);
    if (!datafile) {
        std::cerr << "Couldn't open parmFile: " <<fileName << std::endl;
        exit (-1);
    }

    string parmName, inputStr,sep;
    double parmValue;
    vector<double> parmValueVector;
    ParmMap::iterator mapit;
    
    while (getline(datafile,inputStr)) {
        
        vector<string> tokens;
        boost::split(tokens, inputStr, boost::is_any_of(" "));
        
        unsigned numArg = unsigned(tokens.size());
        
        parmName = tokens[0];
        
        for (unsigned i=1;i<numArg;i++){
            parmValue = stod(tokens[i]);
            parmValueVector.push_back(parmValue);
        }
        
        vector<double> parmV;
        parmV = parmValueVector;
        (*this)[parmName] = parmV;
        parmValueVector.clear();
    }
    datafile.clear();
    datafile.close();
}

void ParmMap::display(void) {
        
    ParmMap::iterator mapit;
    for (mapit=this->begin(); mapit!= this->end(); mapit++)
    {
        std::cout<<mapit->first<<"\t";
        for (auto it : (mapit->second))
        {
            cout<< it <<"\t";
        }
        cout<<endl;
    }
}