//
//  SimPop_Gus.h
//  OSIM1.01
//
//  Created by Hao Cheng on 5/11/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef OSIM1_01_SimPop_Gus_h
#define OSIM1_01_SimPop_Gus_h

#include "SimPop.h"

class SimPop_Gus: public SimPop{
    
    public:
    SimPop_Gus(unsigned numChr, unsigned numLoci,double mutRate):SimPop(numChr, numLoci, mutRate){};

    void sample_Gus(unsigned size, int nGen, string file, unsigned founder_size){
        
        if (founders.size()==0) {
            if(file!=""){
                founders.sampleFounders(founder_size,file);
            }else{
                founders.sampleFounders(size);
            }
        }
        
        parents.sampleChildren(size,founders,founders);
        
        }
    
    MatrixXf getGenotypes(){
        parents.getHaps();
        parents.getNPMatrix();
        return parents.NPMatrix;
    }
    
    void resultClear(){
        
    }

};


#endif
