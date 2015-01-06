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

            unsigned numArg = tokens.size();
            
            parmName = tokens[0];
            
            for (unsigned i=1;i<numArg;i++){
                parmValue = getDouble(tokens[i]);
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
        for (mapit=this->begin(); mapit!= this->end(); mapit++) {
            std::cout<<mapit->first<<"\t";
            
            for (auto it : (mapit->second)){
                cout<< it <<"\t";
            }
            cout<<endl;
        }
    }



 /*   double ParmMap::getDoubleValue(std::string parmName) {
        ParmMap::iterator mapit = find(parmName);
        if (mapit== end()) {
            std::cerr << "Unrecognized parameter: " << parmName << std::endl;
            throw exception("ParmMap::getDoubleValue: Error\n");
        }
        std::string parmValue = mapit->second;
        std::istringstream istr(parmValue);
        double doubleValue;
        istr >> doubleValue;
        return doubleValue;
    }
    
    unsigned ParmMap::getUnsignedValue(std::string parmName) {
        ParmMap::iterator mapit = find(parmName);
        if (mapit== end()) {
            std::cerr << "Unrecognized parameter: " << parmName << std::endl;
            throw exception("ParmMap::getUnsignedValue: Error\n");;
        }
        std::string parmValue = mapit->second;
        std::istringstream istr(parmValue);
        unsigned unsignedValue;
        istr >> unsignedValue;
        return unsignedValue;
    }
    
    std::string ParmMap::getStringValue(std::string parmName) {
        ParmMap::iterator mapit = find(parmName);
        if (mapit== end()) {
            std::cerr << "Unrecognized parameter: " << parmName << std::endl;
            throw exception("ParmMap::getStringValue: Error\n");;
        }
        std::string parmValue = mapit->second;
        
        // remove Mac carriage return from argument...
        size_t foundloc = parmValue.find('\r');
        if (foundloc!=string::npos) {
            std::string tpv;
            tpv = parmValue.substr(0,foundloc);
            parmValue = tpv;
        }		
        
        return parmValue;
    }
    
    
    char* ParmMap::getCharPtr(std::string parmName) {
        ParmMap::iterator mapit = find(parmName);
        if (mapit== end()) {
            std::cerr << "Unrecognized parameter: " << parmName << std::endl;
            throw exception("ParmMap::getCharPtr: Error\n");;
        }
        char* parmValue = (char*)(mapit->second).c_str();
        return parmValue;
    }*/
    
    
