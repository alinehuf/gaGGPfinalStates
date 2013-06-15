//
//  ga.h
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013

#ifndef gaGGPfinalState_ga_h
#define gaGGPfinalState_ga_h

typedef struct { // convenient to sort weights of the best ever genome
  double weight;
  char * fact;
} result;

typedef struct {
  double * vGenes;
  double dFitness;
  double dAdjFitness;
  double dToSpawn;
} genome;

typedef struct {
  genome * vIndiv;
  int iNumIndiv;
  genome sBestEver;
  int iBestEverEpoch;
  double dMinFit, dAvgFit, dMaxFit;
  int iMaxId; // to replace best ever fitness if needed
  double dGlobalFitness;
  int iMaxGenerations;
  int iEpoch;
  float fMutationRate;
  float fCrossoverRate;
} population;

void ga(testBase * data, double cross, double mut, int player, int popsize,
        int genmax);
population * createPopulation(testBase * data, double cross, double mut,
                              int popsize, int genmax);
void freePopulation(population * pop);
genome * createRandomGenomes(int sizePop, int sizeGen);
genome * createGenomes(int sizePop, int sizeGen);
void freeGenomes(genome * indivs, int sizePop);
void getAllFitnesses(population * pop, testBase * data, int player);
double getFitness(genome * indiv, testBase * data, int player);
void adjustAllFitnesses(population * pop);
void dumpWeightAndFitness(population * pop, int sizeGen);
void getNextGeneration(population * pop, testBase * data);
void getNextGeneration2(population * pop, testBase * data);
int selectIndiv(population * pop, int totalToSpawn);
void mutateGenomes(population * pop, testBase * data);
void dumpPopulationState(population * pop, testBase * data);
int cmpResultsByWeight(const result * a, const result * b);

#endif
