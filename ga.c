//
//  ga.c
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013
//
// this file contents the Genetic Algorithm itself.
// ga() creates a population and runs for genmax generations.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "data.h"
#include "utils.h"
#include "ga.h"

extern int Verbose;
extern int UltraVerbose;
extern int Dump;
extern FILE * DumpFile;
extern int Log;
extern FILE * LogFile;
extern int GlobalLog;
extern FILE * GlobalLogFile;


/* the genetic algorithm 
 * popsize individual evolved for genmax generations
 * the fitnes is estimated based on data
 */
void ga(testBase * data, double cross, double mut, int player, int popsize,
        int genmax) {
  //create popultation
  population * pop = createPopulation(data, cross, mut, popsize, genmax);
  // genetic algorithm
  int i;
  for (i = 0; i < pop->iMaxGenerations; i++) {
    pop->iEpoch++;
    getAllFitnesses(pop, data, player);
    getNextGeneration2(pop, data);
    mutateGenomes(pop, data);
    dumpPopulationState(pop, data);
    if (GlobalLog) {
      if (i == 0) {
        fprintf(GlobalLogFile, "start fitness - max: %f - avg: %f - min %f\n",
                pop->dMaxFit, pop->dAvgFit, pop-> dMinFit);
      } else if (i == pop->iMaxGenerations - 1) {
        fprintf(GlobalLogFile, "end   fitness - max: %f - avg: %f - min %f",
                pop->dMaxFit, pop->dAvgFit, pop-> dMinFit);
        fprintf(GlobalLogFile, "\t\t\tbest ever fitness %f found at epoch %d\n",
                pop->sBestEver.dFitness, pop->iBestEverEpoch);
      }
    }
  }
  // free population
  freePopulation(pop);
}


/* creates the population
 */
population * createPopulation(testBase * data, double cross, double mut,
                                                      int popsize, int genmax) {
  // calloc: to initialize all other values ​​to 0
  population * pop = (population *) calloc(1, sizeof(*pop));
  // params
  pop->fCrossoverRate = cross;
  pop->fMutationRate = mut;
  pop->iNumIndiv = popsize; // must be even number
  if (pop->iNumIndiv % 2 == 1) pop->iNumIndiv++;
  pop->iMaxGenerations = genmax;
  // intialize individuals
  pop->vIndiv = createRandomGenomes(pop->iNumIndiv, data->iNumFacts);
  // init best ever individual
  pop->sBestEver.vGenes =
    (double *) calloc(data->iNumFacts,sizeof(*(pop->sBestEver.vGenes)));
  return pop;
}

/* creates sizePop random genomes/individuals
 */
genome * createRandomGenomes(int sizePop, int sizeGen) {
  int i, j;
  genome * indivs = (genome *) malloc(sizePop * sizeof(*indivs));
  srand((unsigned int) time(NULL));
  for (i = 0; i < sizePop; i++) {
    indivs[i].vGenes = (double *) malloc(sizeGen * sizeof(*(indivs[i].vGenes)));
    for (j = 0; j < sizeGen; j++) {
      indivs[i].vGenes[j] = randDouble(-1, 1);
    }
    indivs[i].dFitness = 0;
    indivs[i].dAdjFitness = 0;
  }
  return indivs;
}

/* creates sizePop genomes initialized with zero weigth
 */
genome * createGenomes(int sizePop, int sizeGen) {
  int i;
  genome * indivs = (genome *) calloc(sizePop, sizeof(*indivs));
  for (i = 0; i < sizePop; i++)
    indivs[i].vGenes = (double *) calloc(sizeGen, sizeof(*(indivs[i].vGenes)));
  return indivs;
}

/* free all the genomes of a population
 */
void freeGenomes(genome * indivs, int sizePop) {
  int i;
  for (i = 0; i < sizePop; i++) {
    free(indivs[i].vGenes);
  }
  free(indivs);
}

/* free a population
 */
void freePopulation(population * pop) {
  free(pop->sBestEver.vGenes);
  freeGenomes(pop->vIndiv, pop->iNumIndiv);
  free(pop);
}

/* get the fitness of all the individuals and adjust the fitness span 
 * according to Golberg advice
 */
void getAllFitnesses(population * pop, testBase * data, int player) {
  int i;
  pop->dMinFit = INT_MAX;
  pop->dMaxFit = INT_MIN;
  pop->dAvgFit = 0;
  if (UltraVerbose) puts("Evaluation of each individual :");
  for (i = 0; i < pop->iNumIndiv; i++) {
    pop->dAvgFit += getFitness(&(pop->vIndiv[i]), data, player);
    if (pop->vIndiv[i].dFitness > pop->dMaxFit) {
      pop->dMaxFit = pop->vIndiv[i].dFitness;
      pop->iMaxId = i;
    }
    else if (pop->vIndiv[i].dFitness < pop->dMinFit)
      pop->dMinFit = pop->vIndiv[i].dFitness;
  }
  if (UltraVerbose) printf("\n");
  pop->dAvgFit /= pop->iNumIndiv;
  // remember the best ever fitness
  if (pop->vIndiv[pop->iMaxId].dFitness > pop->sBestEver.dFitness) {
    for(i = 0; i < data->iNumFacts; i++)
      pop->sBestEver.vGenes[i] = pop->vIndiv[pop->iMaxId].vGenes[i];
    pop->sBestEver.dFitness = pop->vIndiv[pop->iMaxId].dFitness;
    pop->iBestEverEpoch = pop->iEpoch;
  }
  adjustAllFitnesses(pop);
  if (UltraVerbose) dumpWeightAndFitness(pop, data->iNumFacts);
}

/* get the fitness of one individual
 */
double getFitness(genome * indiv, testBase * data, int player) {
  int i, j;
  indiv->dFitness = 0;
  double evaluation, error;
  if (UltraVerbose) {
    puts("individual weights :");
    for (i = 0; i < data->iNumFacts; i++)
      printf("%f ", indiv->vGenes[i]);
    puts("");
  }
  for (i = 0; i < data->iNumFinalStates; i++) {
    evaluation = 0;
    if (UltraVerbose) {
      printf("final state %d - list of facts :\n", i);
    }
    for (j = 0; j < data->vAllFinalStates[i].iNumFacts; j++) {
      if (UltraVerbose)
        printf("(%s) * %f\n",
               data->vAllFacts[data->vAllFinalStates[i].vFactIdxs[j]],
               indiv->vGenes[data->vAllFinalStates[i].vFactIdxs[j]] );
      evaluation += indiv->vGenes[data->vAllFinalStates[i].vFactIdxs[j]];
    }
    if (UltraVerbose) {
      printf("evaluation = %f/%d = ", evaluation,
                                      data->vAllFinalStates[i].iNumFacts);
    }
    //evaluation /= data->vAllFinalStates[i].iNumFacts; // inside [-1, 1]
    if (evaluation > 1) evaluation = 1;    // inside [-inf, 1]
    if (evaluation < -1) evaluation = -1;  // inside [-1, 1]
    if (UltraVerbose)  printf("%f => ", evaluation);
    //printf("%f\n", evaluation);
    evaluation = 50 * (evaluation + 1);   // inside [0, 100] lineaire
    //evaluation = 100 * sigmoid(evaluation); // inside [0, 100]
    
    if (UltraVerbose)
      printf("result %f instead of %d - ",
             evaluation, data->vAllFinalStates[i].vScore[player]);
    error = (double) data->vAllFinalStates[i].vScore[player] - evaluation;
    error *= error;
    error /= 100;
    //if (error < 0) error = -error;
    if (UltraVerbose) printf("error = %f\n", error);
    indiv->dFitness += 100 - error;            // fitness for this evaluation
  }
  // final fitness : average of all the final state evaluations
  indiv->dFitness /= data->iNumFinalStates;
  if (UltraVerbose) printf("fitness = %f\n", indiv->dFitness);
  return indiv->dFitness;
}

// linear transformation of Fitness curve
/* Goldberg p79 scaling routine :
 * "for typical population (n = 50 to 100) a Cmult = 1.2 to 2 
 * has been used successfully"
 */
void adjustAllFitnesses(population * pop) {
  const double fitMult = 2;      // Fitness Multiple
  double delta, slope, offset ;  // slope and offset for linear equation
  int i;
  // if the fitness of the population has already converged too quickly, 
  // avoid rescale => risk of dividing by zero (delta=0)
  if (pop->dMaxFit - pop->dAvgFit < 1) {
    for (i = 0 ; i < pop->iNumIndiv ; i++) {
      pop->vIndiv[i].dAdjFitness = pop->vIndiv[i].dFitness;
    }
  } else {
    // compute coefficient for linear scale
    // we want to rescale max to have new_max = CM * avg and new_avg = avg
    // but new_min must be remain a positive number
    if (pop->dMinFit > fitMult * pop->dAvgFit - pop->dMaxFit) { // non-negative ?
      delta = pop->dMaxFit - pop->dAvgFit;     // normal scale
      slope = pop->dAvgFit / delta;
      offset = pop->dAvgFit * (pop->dMaxFit - fitMult * pop->dAvgFit) / delta;
    } else {
      delta = pop->dAvgFit - pop->dMinFit;     // scale as much as possible
      slope = pop->dAvgFit / delta;
      offset = pop->dAvgFit * -pop->dMinFit / delta;
    }
    if (UltraVerbose)
      printf("Adjusted Fitnesses : slope = %f offset = %f\n\n", slope, offset);
    // y = ax + b, translated adaptation value
    pop->dGlobalFitness = 0;
    for (i = 0 ; i < pop->iNumIndiv ; i++) {
      pop->vIndiv[i].dAdjFitness = slope * (pop->vIndiv[i].dFitness) + offset;
      pop->dGlobalFitness += pop->vIndiv[i].dAdjFitness;
    }    
  }
}

// stochastic sampling without replacement (round wheel)
void getNextGeneration(population * pop, testBase * data) {
  // pool containing the individuals for reproduction
  // * 2 because the number of parents on the wheel can be superior of
  // population size
  genome ** pool = (genome* *) malloc(2 * pop->iNumIndiv * sizeof(*pool));
  int poolIdx = 0;
  int i, j, numToSpawn;
  if (UltraVerbose) puts("Quantity of offspring for each individual :");
  for (i = 0 ; i < pop->iNumIndiv ; i++) {
    numToSpawn = (int) (0.5 + pop->vIndiv[i].dAdjFitness * pop->iNumIndiv /
                        pop->dGlobalFitness);
    if (UltraVerbose)
      printf("%d - spawn=%f => %d\n", i, pop->vIndiv[i].dAdjFitness *
                              pop->iNumIndiv / pop->dGlobalFitness, numToSpawn);
    for (j = 0; j < numToSpawn; j++) {
      if (poolIdx == 2 * pop->iNumIndiv) break;
      pool[poolIdx++] = &(pop->vIndiv[i]);
    }
  }

  // selects parents and create the new population
  genome * newPop = createGenomes(pop->iNumIndiv, data->iNumFacts);
  int mum, dad, crossPoint, newPopIdx = 0;
  //  int clone, poolSize = poolIdx;
  while(newPopIdx < pop->iNumIndiv) {
    dad = randInt(0, poolIdx);
    mum = randInt(0, poolIdx);
    // probability of crossover
    if (randDouble(0, 1) < (double) pop->fCrossoverRate) { // crossover
      crossPoint = randInt(0, data->iNumFacts);
    } else {                                               // cloning
      crossPoint = 0;
    }
    // do the crossover
    for (i = 0; i < crossPoint; i++) {
      newPop[newPopIdx].vGenes[i] = pool[dad]->vGenes[i];
      newPop[newPopIdx + 1].vGenes[i] = pool[mum]->vGenes[i];
    }
    for (i = crossPoint; i < data->iNumFacts; i++) {
      newPop[newPopIdx].vGenes[i] = pool[mum]->vGenes[i];
      newPop[newPopIdx + 1].vGenes[i] = pool[dad]->vGenes[i];
    }
    // next offspring
    newPopIdx += 2;
  }
  free(pool);
  freeGenomes(pop->vIndiv, pop->iNumIndiv);
  pop->vIndiv = newPop;
}

// stochastic sampling with replacement (round wheel)
// converges slower and approaches closer to the optimal solution
void getNextGeneration2(population * pop, testBase * data) {
  int i;
  // calculates the quantity to spawn for each individual
  if (UltraVerbose) puts("Quantity of offspring for each individual :");
  for (i = 0 ; i < pop->iNumIndiv ; i++) {
    pop->vIndiv[i].dToSpawn = pop->vIndiv[i].dAdjFitness * pop->iNumIndiv /
                                                            pop->dGlobalFitness;
    if (UltraVerbose) printf("%d - spawn=%f\n", i, pop->vIndiv[i].dToSpawn);
  }

  // selects parents and create the new population
  genome * newPop = createGenomes(pop->iNumIndiv, data->iNumFacts);
  int mum, dad, crossPoint, newPopIdx = 0;
  while(newPopIdx < pop->iNumIndiv) {
    dad = selectIndiv(pop, pop->iNumIndiv - newPopIdx);
    mum = selectIndiv(pop, pop->iNumIndiv - newPopIdx - 1);

    // probability of crossover
    if (randDouble(0, 1) < (double) pop->fCrossoverRate) { // crossover
      crossPoint = randInt(0, data->iNumFacts);
    } else {                                               // cloning
      crossPoint = 0;
    }
    // do the crossover
    for (i = 0; i < crossPoint; i++) {
      newPop[newPopIdx].vGenes[i] = pop->vIndiv[dad].vGenes[i];
      newPop[newPopIdx + 1].vGenes[i] = pop->vIndiv[mum].vGenes[i];
    }
    for (i = crossPoint; i < data->iNumFacts; i++) {
      newPop[newPopIdx].vGenes[i] = pop->vIndiv[mum].vGenes[i];
      newPop[newPopIdx + 1].vGenes[i] = pop->vIndiv[dad].vGenes[i];
    }
    // next offspring
    newPopIdx += 2;
  }
  freeGenomes(pop->vIndiv, pop->iNumIndiv);
  pop->vIndiv = newPop;
}

int selectIndiv(population * pop, int totalToSpawn) {
  double wheelslot = randDouble(0, totalToSpawn - 1);
  double position = 0;
  int i = 0;
  
  while (position <= wheelslot && i < pop->iNumIndiv)
    position += pop->vIndiv[i++].dToSpawn;

  i--;
  pop->vIndiv[i].dToSpawn--;
  if ( pop->vIndiv[i].dToSpawn < 0) pop->vIndiv[i].dToSpawn = 0;
  return i;
}

void mutateGenomes(population * pop, testBase * data) {
  int i, numMutation, span, rdm1, rdm2;
  double tmp;
  for (i = 0 ; i < pop->iNumIndiv ; i++) {
    numMutation = randInt(0, data->iNumFacts) * pop->fMutationRate;
    while(numMutation--) {
      rdm1 = randInt(0, data->iNumFacts);
//      if (randInt(0, 3) == 0) {
//        pop->vIndiv[i].vGenes[randInt(0, data->iNumFacts)] = randDouble(-1,1);
//      } else
      if (randInt(0, 2) == 0) {
        span = randDouble(0.1, 0.5); // add a value in [-x, x], x E [0.1, 0.5]
        //span = 0.05;                 // add a value in [-0.05, 0.05]
        pop->vIndiv[i].vGenes[rdm1] += randDouble(-span, span);
        if (pop->vIndiv[i].vGenes[rdm1] > 1.0)
          pop->vIndiv[i].vGenes[rdm1] = 1.0;
        if (pop->vIndiv[i].vGenes[rdm1] < -1.0)
          pop->vIndiv[i].vGenes[rdm1] = -1.0;
//      } else if (randInt(0, 1) == 0) { // sign inversion
//        pop->vIndiv[i].vGenes[rdm1] = -pop->vIndiv[i].vGenes[rdm1];
      } else {                         // inversion of the weight of two facts
        rdm2 = randInt(0, data->iNumFacts);
        tmp = pop->vIndiv[i].vGenes[rdm1];
        pop->vIndiv[i].vGenes[rdm1] = pop->vIndiv[i].vGenes[rdm2];
        pop->vIndiv[i].vGenes[rdm2] = tmp;
      }
    }
  }
}

void dumpWeightAndFitness(population * pop, int sizeGen) {
  int i,j;
  printf("Fitmess and weights : \n");
  for (i = 0; i < pop->iNumIndiv; i++) {
    printf("Fitness=%6.2f AdjFitness=%6.2f weights :\n",
           pop->vIndiv[i].dFitness, pop->vIndiv[i].dAdjFitness);
    for (j = 0; j < sizeGen; j++) {
      printf("%f ", pop->vIndiv[i].vGenes[j]);
    }
    if (pop->vIndiv[i].dFitness == pop->vIndiv[pop->iMaxId].dFitness)
      printf(" *");
    printf("\n");
  }
  printf("\n");
}

void dumpPopulationState(population * pop, testBase * data) {
  // log file
  if (Log) {
    fprintf(LogFile, "------------ EPOCH %d ------------\n", pop->iEpoch);
    fprintf(LogFile, "fitness - max: %f - avg: %f - min %f\n",
            pop->dMaxFit, pop->dAvgFit, pop-> dMinFit);
    fprintf(LogFile, "best ever fitness %f found at epoch %d\n",
            pop->sBestEver.dFitness, pop->iBestEverEpoch);
  }
  // dump file
  if (Dump)
    fprintf(DumpFile, "%f\t%f\t%f\n", pop->dMaxFit, pop->dAvgFit, pop->dMinFit);
  // standard output
  if (Verbose) {
    printf("------------ EPOCH %d ------------\n", pop->iEpoch);
    printf("fitness - max: %f - avg: %f - min %f\n",
           pop->dMaxFit, pop->dAvgFit, pop-> dMinFit);  
    printf("best ever fitness %f found at epoch %d\n", pop->sBestEver.dFitness,
                                                       pop->iBestEverEpoch);
  }
  // details of best ever individual
  if (Verbose || Log) {
    if (pop->iEpoch % 100 == 1) {
      result * results = (result *) malloc(data->iNumFacts * sizeof(*results));
      int i;
      for (i = 0; i < data->iNumFacts; i++) {
        results[i].fact = data->vAllFacts[i];
        results[i].weight = pop->sBestEver.vGenes[i];
      }
      qsort(results, data->iNumFacts, sizeof(*results),
                       (int (*) (const void *, const void *)) cmpResultsByWeight);
      for (i = 0; i < data->iNumFacts; i++) {
        if (Log)
          fprintf(LogFile, "(%s) => %f\n", results[i].fact, results[i].weight);
        if (Verbose)
          printf("(%s) => %f\n", results[i].fact, results[i].weight);
      }
      free(results);
    }
    if (Verbose) printf("\n");
  }
}

int cmpResultsByWeight(const result * a, const result * b) {
  if (a->weight == b->weight) return 0;
  else if (a->weight < b->weight) return -1;
  else return 1;
}

