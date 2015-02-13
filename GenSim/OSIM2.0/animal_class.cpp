//
//  animal_class.cpp
//  OSIM
//
//  Created by Hao Cheng on 3/17/14.
//  Copyright (c) 2014 Hao. All rights reserved.
//

#include "animal_class.h"
#include "genome_info.h"


//make(create) an animal as a founder
void AnimalClass::sampleFounder()
{
    sireId=0;
    damId=0;
    initFounderPosOriMut();
    initFounderHaps();
    
}

void AnimalClass::sampleFounder(vector<string> tokens1,vector<string> tokens2)
{
    sireId=0;
    damId=0;
    initFounderPosOriMut();
    inputFounderPatHaps(tokens1);
    inputFounderMatHaps(tokens2);
    
}


//initialize arrays for pos and ori in founders
void AnimalClass::initFounderPosOriMut()
{
    unsigned numChromosomePair = G.get_num_chrom();

    for(unsigned i=0;i<numChromosomePair;i++){
        GenomePat[i].ori.resize(1);
        GenomePat[i].ori[0]=countChromosome;
        GenomePat[i].pos.resize(1);
        GenomePat[i].pos[0]=0;
        GenomePat[i].mut.resize(0);


        GenomeMat[i].ori.resize(1);
        GenomeMat[i].ori[0]=countChromosome+1;
        GenomeMat[i].pos.resize(1);
        GenomeMat[i].pos[0]=0;
        GenomeMat[i].mut.resize(0);

    }
    countChromosome += 2;
}


//initialize paternal and maternal haplotype arrays in one founder randomly
void AnimalClass::initFounderHaps()
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
    
        unsigned numLoci=G[i].get_num_loci();
        GenomePat[i].haplotype.resize(numLoci);
        GenomeMat[i].haplotype.resize(numLoci);
        
        for(unsigned j=0;j<numLoci;j++){
            binomial_distribution<int> Binom(1,G[i][j].allele_freq(0));
            GenomePat[i].haplotype[j]= Binom(randGen);
            GenomeMat[i].haplotype[j]= Binom(randGen);
        }
        
    }

}

//initialize the paternal haplotype array in one founder from a file
void AnimalClass::inputFounderPatHaps(vector<string> tokens)
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
        
        unsigned numLoci=G[i].get_num_loci();
        GenomePat[i].haplotype.resize(numLoci);
        
        for(unsigned j=0;j<numLoci;j++){
            GenomePat[i].haplotype[j]= getDouble(tokens[j]);
        }
    }
}

//initialize the maternal haplotype array in one founder from a file
void AnimalClass::inputFounderMatHaps(vector<string> tokens)
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
        
        unsigned numLoci=G[i].get_num_loci();
        GenomeMat[i].haplotype.resize(numLoci);
        
        for(unsigned j=0;j<numLoci;j++){
            GenomeMat[i].haplotype[j]= getDouble(tokens[j]);
        }
        
    }
    
}


//get new Ori and Pos arrays for the child after recombination
void AnimalClass::sampleMyPosOriMut(AnimalClass& father,
                               AnimalClass& mother)
{
    if(!G.mapPosDone) G.mkMapPos();
    
    sireId  = father.myId;
    damId   = mother.myId;
	
	chromosome *currentFatherChrom;
	chromosome *currentMotherChrom;
	
	unsigned numChromosomePair = G.get_num_chrom();
    
    //construct paternal chromosome for myId
    ArrayXf tempPos;
    ArrayXi tempOri;
    ArrayXf tempMut;
    tempPos.resize(100000);
    tempOri.resize(100000);
    tempMut.resize(100000);
    
    for(unsigned i=0;i<numChromosomePair;i++){
        double chrLength = G[i].chr_length;
        unsigned binomialN = chrLength*3 + 1;
        vector<float> rPos;
        binomial_distribution<int> Binom(binomialN,chrLength/binomialN);
        int numCrossover    =   Binom(randGen);
        
        uniform_real_distribution<float> u(0,1);
        for (unsigned k=0; k<numCrossover; k++) {
            rPos.push_back(chrLength*u(randGen));
        }
        rPos.push_back(chrLength);
        sort(rPos.begin(),rPos.end());
        

        unsigned   startPosParent=0;
        unsigned   startPosMe=0;
        unsigned   startPosMe_mut=0;
        unsigned   startPosParent_mut=0;


        currentFatherChrom = (u(randGen)<0.5)?&father.GenomePat[i]:&father.GenomeMat[i];

        for(unsigned j=0;j<rPos.size();j++){

            unsigned numCopy=(currentFatherChrom->pos < rPos[j]).count()-startPosParent;
            tempPos.block(startPosMe,0,numCopy,1)=currentFatherChrom->pos.block(startPosParent,0,numCopy,1);
            tempOri.block(startPosMe,0,numCopy,1)=currentFatherChrom->ori.block(startPosParent,0,numCopy,1);
            
            unsigned numCopy_mut=(currentFatherChrom->mut < rPos[j]).count()-startPosParent_mut;
            tempMut.block(startPosMe_mut,0,numCopy_mut,1)=currentFatherChrom->mut.block(startPosParent_mut,0,numCopy_mut,1);

            currentFatherChrom= (currentFatherChrom==&father.GenomePat[i])?&father.GenomeMat[i]:&father.GenomePat[i];
            
            startPosParent=(currentFatherChrom->pos < rPos[j]).count();
            startPosMe+=numCopy;
            tempPos(startPosMe)=rPos[j];
            tempOri(startPosMe)=currentFatherChrom->ori(startPosParent-1);
            startPosMe++;

            startPosParent_mut=(currentFatherChrom->mut < rPos[j]).count();
            startPosMe_mut+=numCopy_mut;

        } //the length of tempPos should be startPosMe now. Tricky here.
        
       
        unsigned keep=0;
        GenomePat[i].pos.resize(startPosMe);
        GenomePat[i].ori.resize(startPosMe);
        GenomePat[i].pos[0]=tempPos[0];
        GenomePat[i].ori[0]=tempOri[0];

        for(unsigned m=1;m < startPosMe-1; m++){

            if(tempOri[m]!=tempOri[m-1]){
            keep=keep+1;
            GenomePat[i].pos[keep]=tempPos[m];
            GenomePat[i].ori[keep]=tempOri[m];
            }
        }
        
        unsigned PosOriSize=keep+1;
        GenomePat[i].pos.conservativeResize(PosOriSize);
        GenomePat[i].ori.conservativeResize(PosOriSize);

        tempMut.conservativeResize(startPosMe_mut);
        
        //sample mutaions
        
        unsigned nLoci=G[i].get_num_loci();//num of loci is euqal to num of mapPos
        binomial_distribution<int> Binom_mut(nLoci,G.mutRate);
        int nMut=Binom_mut(randGen);
        
        vector<float> mutPos;
        for (unsigned k=0; k<nMut; k++) {
            unsigned which  = unsigned(nLoci*u(randGen));
            if(which!=nLoci){
                mutPos.push_back(G[i].MapPos(which));
            }
        }
        
        nMut=mutPos.size();
        sort(mutPos.begin(),mutPos.end());
        
        GenomePat[i].mut.resize(startPosMe_mut+nMut);
        
        unsigned start_mut=0;
        unsigned start_temp=0;
        
        if(nMut){
            
            for (unsigned k=0; k<nMut; k++) {
                
                unsigned numcopy=(tempMut < mutPos[k]).count()-start_temp;
                GenomePat[i].mut.block(start_mut,0,numcopy,1)=tempMut.block(start_temp,0,numcopy,1);
                start_mut=start_mut+numcopy;
                start_temp=start_temp+numcopy;
                GenomePat[i].mut[start_mut]= mutPos[k];     //get map position of mutation loci
                start_mut++;
                
            }
        }
        
        GenomePat[i].mut.conservativeResize(start_mut);

    }
    ///
    ///construct maternal chromosome for myId
    tempMut.resize(100000);
    for(unsigned i=0;i<numChromosomePair;i++){
        double chrLength = G[i].chr_length;
        unsigned binomialN = chrLength*3 + 1;
        vector<float> rPos;
        binomial_distribution<int> Binom(binomialN,chrLength/binomialN);
        int numCrossover    =   Binom(randGen);
        
        uniform_real_distribution<float> u(0,1);
        for (unsigned k=0; k<numCrossover; k++) {
            float rPos_Hao=chrLength*u(randGen);
                rPos.push_back(rPos_Hao);
        }

        rPos.push_back(chrLength);
        sort(rPos.begin(),rPos.end());
        
        unsigned   startPosParent=0;
        unsigned   startPosMe=0;
        unsigned   startPosMe_mut=0;
        unsigned   startPosParent_mut=0;
   
        currentMotherChrom = (u(randGen)<0.5)?&mother.GenomePat[i]:&mother.GenomeMat[i];
        
        for(unsigned j=0;j<rPos.size();j++){
            
            unsigned numCopy=(currentMotherChrom->pos < rPos[j]).count()-startPosParent;
            tempPos.block(startPosMe,0,numCopy,1)=currentMotherChrom->pos.block(startPosParent,0,numCopy,1);
            tempOri.block(startPosMe,0,numCopy,1)=currentMotherChrom->ori.block(startPosParent,0,numCopy,1);
            
            unsigned numCopy_mut=(currentMotherChrom->mut < rPos[j]).count()-startPosParent_mut;
            tempMut.block(startPosMe_mut,0,numCopy_mut,1)=currentMotherChrom->mut.block(startPosParent_mut,0,numCopy_mut,1);
            
            currentMotherChrom= (currentMotherChrom==&mother.GenomePat[i])?&mother.GenomeMat[i]:&mother.GenomePat[i];
            
            startPosParent=(currentMotherChrom->pos < rPos[j]).count();
            startPosMe+=numCopy;
            tempPos(startPosMe)=rPos[j];
            tempOri(startPosMe)=currentMotherChrom->ori(startPosParent-1);
            startPosMe++;

            startPosParent_mut = (currentMotherChrom->mut < rPos[j]).count();
            startPosMe_mut    += numCopy_mut;
            
            
        }
        
        unsigned keep=0;
        GenomeMat[i].pos.resize(startPosMe);
        GenomeMat[i].ori.resize(startPosMe);
        GenomeMat[i].pos[0]=tempPos[0];
        GenomeMat[i].ori[0]=tempOri[0];
        
        
        for(unsigned m=1;m < startPosMe-1; m++){
            
            if(tempOri[m]!=tempOri[m-1]){
                keep=keep+1;
                GenomeMat[i].pos[keep]=tempPos[m];
                GenomeMat[i].ori[keep]=tempOri[m];
            }
            
        }
        
        unsigned PosOriSize=keep+1;
        GenomeMat[i].pos.conservativeResize(PosOriSize);
        GenomeMat[i].ori.conservativeResize(PosOriSize);
        
        tempMut.conservativeResize(startPosMe_mut);
        
        //sample mutaions
        
        unsigned nLoci=G[i].get_num_loci();//num of loci is euqal to num of mapPos
        binomial_distribution<int> Binom_mut(nLoci,G.mutRate);
        int nMut=Binom_mut(randGen);
        
        vector<float> mutPos;
        for (unsigned k=0; k<nMut; k++) {
            unsigned which  = unsigned(nLoci*u(randGen));//which can be equal to nLoci
            if(which!=nLoci){
                mutPos.push_back(G[i].MapPos(which));
            }
        }
        
        nMut=mutPos.size();
        
        sort(mutPos.begin(),mutPos.end());
        
        GenomeMat[i].mut.resize(startPosMe_mut+nMut);

        unsigned start_mut=0;
        unsigned start_temp=0;
        
        if(nMut){
            
            for (unsigned k=0; k<nMut; k++) {
                
                unsigned Hao=(tempMut < mutPos[k]).count();
                    unsigned numcopy=(tempMut < mutPos[k]).count()-start_temp;
                    GenomeMat[i].mut.block(start_mut,0,numcopy,1)=tempMut.block(start_temp,0,numcopy,1);
                    start_mut=start_mut+numcopy;
                    start_temp=start_temp+numcopy;
                    GenomeMat[i].mut[start_mut]=mutPos[k];     //get map position of mutation loci
                    start_mut++;
            }
        }

        GenomeMat[i].mut.conservativeResize(start_mut);

    }
}



//make(create) an animal as a child (non-founder)
void AnimalClass::sampleNonFounder(AnimalClass& father, AnimalClass& mother)
{
    sireId=father.myId;
    damId=mother.myId;
    sampleMyPosOriMut(father,mother);
};



//get the haplotype array in one animal based on pos/ori arrays and mutation information for this animal
void AnimalClass::getMyHaps() //&
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++){
        
        unsigned numLoci=G[i].get_num_loci();
        GenomePat[i].haplotype.resize(numLoci);
        GenomeMat[i].haplotype.resize(numLoci);

        ///For paternal halotype
        unsigned numOriPat = GenomePat[i].ori.size();
        
        
        for(unsigned segment=0;segment<numOriPat;segment++){
            unsigned   numcopy;
            unsigned   start=(G[i].MapPos<GenomePat[i].pos[segment]).count();
            int        myOri=GenomePat[i].ori[segment];
            
            if(segment!= numOriPat-1){
                numcopy=(G[i].MapPos<GenomePat[i].pos[segment+1]).count()-start;
            }else{
                numcopy=(G[i].MapPos<G[i].chr_length).count()-start;
            }

            unsigned whichFounder=unsigned(myOri/2);//i.e. which Ori
            chromosome* GenomePatOrMatInThisFounder= (myOri%2)? &AnimalClass::founders[whichFounder]->GenomeMat[i]: &AnimalClass::founders[whichFounder]->GenomePat[i];        //i.e. GenomePat or GenomeMat in this founder "whichFounder"
            
            GenomePat[i].haplotype.block(start,0,numcopy,1)=GenomePatOrMatInThisFounder->haplotype.block(start,0,numcopy,1);
        }


        for(unsigned mut_i=0;mut_i<GenomePat[i].mut.size();mut_i++){
            
            unsigned mutLocus=(G[i].MapPos < GenomePat[i].mut(mut_i)).count();
            
            if(GenomePat[i].haplotype(mutLocus)==1){
                GenomePat[i].haplotype(mutLocus)=0;
            }else{
                GenomePat[i].haplotype(mutLocus)=1;
            }
        }

        ///For maternal haplotype
        unsigned numOriMat = GenomeMat[i].ori.size();
        
        for(unsigned segment=0;segment<numOriMat;segment++){
            
            unsigned   numcopy;
            unsigned   start=(G[i].MapPos<GenomeMat[i].pos[segment]).count();
            int        myOri= GenomeMat[i].ori[segment];
            
            if(segment!=numOriMat-1){
                numcopy=(G[i].MapPos<GenomeMat[i].pos[segment+1]).count()-start;
            }else{
                numcopy=(G[i].MapPos<G[i].chr_length).count()-start;
            }
            
            unsigned whichFounder=unsigned(myOri/2);//i.e. which Ori
            chromosome* GenomePatOrMatInThisFounder= (myOri%2)? &AnimalClass::founders[whichFounder]->GenomeMat[i]: &AnimalClass::founders[whichFounder]->GenomePat[i];        //i.e. GenomePat or GenomeMat in this founder "whichFounder"
            
            GenomeMat[i].haplotype.block(start,0,numcopy,1)=GenomePatOrMatInThisFounder->haplotype.block(start,0,numcopy,1);
        }
        
        
        for(unsigned mut_i=0;mut_i<GenomeMat[i].mut.size();mut_i++){
            
            unsigned mutLocus=(G[i].MapPos < GenomeMat[i].mut(mut_i)).count();
            
            if(GenomeMat[i].haplotype(mutLocus)==1){
                GenomeMat[i].haplotype(mutLocus)=0;
            }else{
                GenomeMat[i].haplotype(mutLocus)=1;
            }
        }


    }
}


//after getMyhaps,my genotype is got by summing paternal and maternal haplotypes
void AnimalClass::getMyGenotype()
{
    unsigned numChromosomePair = G.get_num_chrom();
    unsigned numTotalLoci=0;
    
    for(unsigned i=0;i<numChromosomePair;i++){
        numTotalLoci+=G[i].get_num_loci();
    }
    myGenotype.resize(numTotalLoci);

    unsigned start=0;
    for(unsigned i=0;i<numChromosomePair;i++){
        
        unsigned numLoci=G[i].get_num_loci();
        myGenotype.block(start,0,numLoci,1)= GenomePat[i].haplotype+GenomeMat[i].haplotype;
        start+=numLoci;
    }
    
}


//display the ori/pos arrays for this animal
void AnimalClass::display()
{
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    unsigned numChromosomePair = G.get_num_chrom();
    for(unsigned i=0;i<numChromosomePair;i++){
        cout<<"These are positions on chromosome "<<i<<" for pat"<<endl;
        cout<< GenomePat[i].pos.format(CommaInitFmt)<<endl;
        cout<<"These are origins on chromosome "<<i<<" for pat"<<endl;
        cout<< GenomePat[i].ori.format(CommaInitFmt)<<endl;
        cout<<"These are positions on chromosome "<<i<<" for mat"<<endl;
        cout<< GenomeMat[i].pos.format(CommaInitFmt)<<endl;
        cout<<"These are origins on chromosome "<<i<<" for mat"<<endl;
        cout<< GenomeMat[i].ori.format(CommaInitFmt)<<endl;
    }
};

//display the number of ori segments in one animal
unsigned AnimalClass::displayNumPos()
{
    unsigned sumSizePat=0;
    unsigned sumSizeMat=0;
    for(unsigned i=0;i<G.get_num_chrom();i++){
        sumSizePat+=GenomePat[i].pos.size();
        sumSizeMat+=GenomeMat[i].pos.size();
    }
    cout<<sumSizeMat+sumSizePat<<endl;
    return((sumSizePat+sumSizeMat));
};














    
