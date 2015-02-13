//
//  global.h
//  
//
//  Created by Rohan L Fernando on 7/9/14.
//
//

#ifndef _global_h
#define _global_h

default_random_engine AnimalClass::randGen;
Genome_info AnimalClass::G;
unsigned AnimalClass::countChromosome=0;
unsigned AnimalClass::countId=0;
vector<AnimalClass*> AnimalClass::founders;
vector<mutantInfo*> AnimalClass::mutants;


#endif
