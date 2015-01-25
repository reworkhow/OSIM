//
//  parmMap.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/13/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//  learn parmMap from MATVEC

#ifndef __GenSim1_3__parmMap__
#define __GenSim1_3__parmMap__

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

class ParmMap: public std::map<std::string,vector<double> >
{
    public:
        void inputParms(std::string fileName);
        void display(void);
        void display(std::string fileName);
};

#endif /* defined(__GenSim1_3__parmMap__) */
