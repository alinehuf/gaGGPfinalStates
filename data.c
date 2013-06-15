//
//  data.c
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013
//
// this files creates a structure for storing the data describing all the final
// states of a game. It reads the data file ans return a testBase struct.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

testBase * loadData(const char * filename) {
  FILE * fileref;
  if ((fileref = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "error : unable to open dump file %s", filename);
    exit(EXIT_FAILURE);
  }
  // test database structure
  // calloc: to initialize all other values ​​to 0
  testBase * data = (testBase *) calloc(1, sizeof(*data));
  data->iMaxFacts = INIT_VECT_SIZE;
  data->vAllFacts =
                  (char* *) malloc(INIT_VECT_SIZE * sizeof(*(data->vAllFacts)));
  data->iMaxFinalStates = INIT_VECT_SIZE;
  data->vAllFinalStates =
       (finalState *) malloc(INIT_VECT_SIZE * sizeof(*(data->vAllFinalStates)));

  int i;
  int idxFact = -1;
  char lineBuffer[MAX_LINE_SIZE];     // a line of the data file
  int charIdx;                        // id of the current char in the line
  char scoreBuffer[MAX_LINE_SIZE];    // buffer for the string of scores
  char buffer[MAX_LINE_SIZE];         // general buffer (for differents facts)

  while(fgets(lineBuffer, MAX_LINE_SIZE, fileref) != NULL) {
    // new final state
    finalState fstate;
    
    // read score and player strings
    charIdx = readGroup(lineBuffer, 0, scoreBuffer);      // score
    
    // get the score
    fstate.vScore = readScore(scoreBuffer, &data->iNumPlayers);

    // get the facts => final state of the game
    // get the number of facts
    fstate.iNumFacts = getNumberFacts(&lineBuffer[charIdx]);
    // allocate memory
    fstate.vFactIdxs =
                 (int *) malloc(fstate.iNumFacts * sizeof(*(fstate.vFactIdxs)));
    i = 0;
    while ((charIdx = readGroup(lineBuffer, charIdx, buffer)) != -1) {
      // fact already in the base
      if ((idxFact = getFactPointer(buffer, data)) != -1)
        fstate.vFactIdxs[i++] = idxFact;
      // else add it
      else fstate.vFactIdxs[i++] = addFact(buffer, data);
    }
    addFinalState(&fstate, data);
  }
  fclose(fileref);
  return data;
}

void freeData(testBase * data) {
  int i;
  for (i = 0; i < data->iNumFacts; i++)
    free(data->vAllFacts[i]);
  free(data->vAllFacts);
  for (i = 0; i < data->iNumFinalStates; i++) {
    free(data->vAllFinalStates[i].vScore);
    free(data->vAllFinalStates[i].vFactIdxs);
  }
  free(data->vAllFinalStates);
  free(data);
}

int readGroup(char * text, int textIdx, char * buffer) {
  int bufferIdx = 0;
  // buffer is NULL
  if (buffer == 0) {
    fprintf(stderr, "error in readGroup() : buffer is NULL");
    exit(EXIT_FAILURE);
  }
  // text is NULL
  if (text == NULL || text[textIdx] == '\0') {
    buffer[0] = '\0';
    return -1;
  }
  // find the opening parenthesis
  while (text[textIdx] != '\0' && text[textIdx] != '(')
    textIdx++;
  while (text[textIdx] != '\0' && text[textIdx] == '(')
    textIdx++;
  // copy the text until the closing parenthesis
  while (text[textIdx] != '\0' && text[textIdx] != ')')
    buffer[bufferIdx++] = text[textIdx++];
  // nothing in the buffer
  if (bufferIdx == 0) return -1;
  // else add end string char
  buffer[bufferIdx] = '\0';
  return textIdx;
}

int * readScore(char * scoreBuffer, int * iNumPlayers) {
  int i;
  int intBuffer[MAX_PLAYERS_NUMBER];
  int numScores = 0;
  char * number;

  while((number = strsep(&scoreBuffer, " \t")) != NULL) {
    if ((i = sscanf(number, "%d", &intBuffer[numScores++])) < 1) {
      numScores--;
      break;
    }
  }
  
  if (*iNumPlayers != 0 && numScores != *iNumPlayers) {
    fprintf(stderr, "error in readScore(): different number of scores");
    exit(EXIT_FAILURE);
  }
  *iNumPlayers = numScores;
  if (*iNumPlayers == 0)
    fprintf(stderr, "error in readScore(): no score found");
  int * scores = malloc(numScores * sizeof(*scores));
  memcpy(scores, intBuffer, numScores * sizeof(*scores));
  return scores;
}

int getNumberFacts(char * buffer) {
  int i = 0;
  int nb = 0;
  while (buffer[i] != 0)
    if (buffer[i++] == '(') nb++;
  return nb;
}

int getFactPointer(char * buffer, testBase * data) {
  int i;
  for (i = 0; i < data->iNumFacts; i++)
    if (strcmp(buffer, data->vAllFacts[i]) == 0)
      return i;
  return -1;
}

int  addFact(char * buffer, testBase * data) {
  if (data->iNumFacts == data->iMaxFacts) {
    data->iMaxFacts += INIT_VECT_SIZE;
    data->vAllFacts = (char* *) realloc(data->vAllFacts,
                                  data->iMaxFacts * sizeof(*(data->vAllFacts)));
  }
  data->vAllFacts[data->iNumFacts] = strdup(buffer);
  data->iNumFacts++;
  return data->iNumFacts - 1;
}

void addFinalState(finalState * fstate, testBase * data) {
  if (data->iNumFinalStates == data->iMaxFinalStates) {
    data->iMaxFinalStates += INIT_VECT_SIZE;
    data->vAllFinalStates = (finalState *) realloc(data->vAllFinalStates,
                      data->iMaxFinalStates * sizeof(*(data->vAllFinalStates)));
  }
  data->vAllFinalStates[data->iNumFinalStates++] = *fstate;
}

void dumpFacts(testBase * data) {
  int i;
  printf("%d differents facts :\n", data->iNumFacts);
  for (i = 0; i < data->iNumFacts; i++)
    printf("%d - (%s)\n", i, data->vAllFacts[i]);
  printf("\n");
}

void dumpFinalStates(testBase * data) {
  int i, j;
  puts("All final states :");
  for (i = 0; i < data->iNumFinalStates; i++) {
    for (j = 0; j < data->iNumPlayers; j++)
      printf("%d ", data->vAllFinalStates[i].vScore[j]);
    printf("=> ");
    for (j = 0; j < data->vAllFinalStates[i].iNumFacts; j++)
      printf("(%s) ", data->vAllFacts[data->vAllFinalStates[i].vFactIdxs[j]]);
    printf("\n");
  }
  printf("\n");
}
