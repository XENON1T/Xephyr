#include "ToyGenerator.h"

ToyGenerator::ToyGenerator(TString sampleName, TString outDir):errorHandler("ToyGenerator"){

    treeName = sampleName;
    dir = outDir;
    averageCalEvnt  = -9;
    averageDataEvnt = -9;
    likeHood = NULL;

}


void ToyGenerator::setSeed(int seed){

    randomizeMyass.SetSeed(seed);
}


void ToyGenerator::generateCalibration(int N){

}


void ToyGenerator::generateData(double mu, int N){

}

void ToyGenerator::randomizeNuissanceParameter(){

}
