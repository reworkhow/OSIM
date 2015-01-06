//
//  popMap.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/14/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef __GenSim1_3__popMap__
#define __GenSim1_3__popMap__

#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <exception>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include "tools.h"


using namespace std;
using namespace boost;


class PopMap: public std::map<std::string,vector<string>> {
public:
    void inputParms(std::string fileName);
    double getDoubleValue(std::string parmName);
    unsigned getUnsignedValue(std::string parmName);
    std::string getStringValue(std::string parmName);
    char* getCharPtr(std::string parmName);
    void display(void);
    void display(std::string fileName);
};

#endif /* defined(__GenSim1_3__popMap__) */
