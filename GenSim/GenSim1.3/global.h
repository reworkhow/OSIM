//
//  global.h
//  GenSim1.3
//
//  Created by Hao Cheng on 12/12/14.
//  Copyright (c) 2014 Hao Cheng. All rights reserved.
//

#ifndef GenSim1_3_global_h
#define GenSim1_3_global_h

default_random_engine AnimalClass::randGen;
Genome_info AnimalClass::G;
unsigned AnimalClass::countChromosome=0;
unsigned AnimalClass::countId=0;
vector<AnimalClass*> AnimalClass::founders;
vector<mutantInfo*> AnimalClass::mutants;


#endif
