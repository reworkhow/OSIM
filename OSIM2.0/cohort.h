//
//  cohort.h
//  OSIM
//
//  Created by Hao Cheng on 3/19/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#ifndef __OSIM__cohort__
#define __OSIM__cohort__

#include <iostream>
#include "animal_class.h"


class cohort:public vector<AnimalClass*> {
public:
    void sampleFounders(unsigned numAnimals);
    void sampleFounders(unsigned numAnimals, string file);
    void sampleChildren(unsigned numAnimals,cohort &fathers, cohort &mothers);
    AnimalClass* getRandomInd(void);
    void display();
    void flush();
    unsigned displaySumNumPos();//average pos for different generations
    void getHaps();
    void copy(cohort &from);
    
    MatrixXf NPMatrix;
    void getNPMatrix();
};




#endif /* defined(__OSIM__cohort__) */
