//
//  data.h
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013

#ifndef gaGGPfinalState_loaddata_h
#define gaGGPfinalState_loaddata_h

#define MAX_LINE_SIZE 4096
#define INIT_VECT_SIZE 1000
#define MAX_PLAYERS_NUMBER 10

typedef enum {FALSE = 0, TRUE} bool;

typedef struct {
  int * vFactIdxs;
  int * vScore;
  int iNumFacts;
} finalState;

typedef struct {
  int iNumPlayers;
  char ** vAllFacts;
  int iNumFacts;
  int iMaxFacts;
  finalState * vAllFinalStates;
  int iNumFinalStates;
  int iMaxFinalStates;
} testBase;

/* prototypes */
testBase * loadData(const char * filename);
void freeData(testBase * data);
int readGroup(char * text, int textIdx, char * buffer);
int * readScore(char * scoreBuffer, int * nbPlayers);
int getNumberFacts(char * buffer);
int getFactPointer(char * buffer, testBase * data);
int addFact(char * buffer, testBase * data);
void addFinalState(finalState * fstate, testBase * data);
void dumpFacts(testBase * data);
void dumpFinalStates(testBase * data);

#endif
