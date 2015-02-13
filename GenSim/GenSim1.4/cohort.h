//
//  cohort.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef __GenSim1_3__cohort__
#define __GenSim1_3__cohort__


#include <iostream>
#include "animal_class.h"
#include "ped.h"


class cohort:public vector<AnimalClass*> {
public:
    void sampleFounders(unsigned numAnimals);
    void sampleFounders(unsigned numAnimals, string hapFile);
    void sampleChildren(unsigned numAnimals,cohort &fathers, cohort &mothers);
    void sampleChildrenWithPedigree(vector<pedLine> pedTable,cohort &fathers, cohort &mothers);
    
    AnimalClass* getRandomInd(void);
    void display();
    void showIds();
    void flush();
    void getHaps();
    void copy(cohort &from);
    
    MatrixXf NPMatrix;
    MatrixXf NPMatrixQTL;
    void getNPMatrix();
    void getNPMatrixQTL(unsigned numQTL);
    void getPhenotype(unsigned numQTL,double heretability);
    
    VectorXf phenotype;
    VectorXf breedingValue;
    
    unsigned displaySumNumPos();//average pos for different generations
};

#endif /* defined(__GenSim1_3__cohort__) */
