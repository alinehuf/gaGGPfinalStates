//
//  main.c
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013
//
// this file manage the differents tests
// it takes the content of game_dir and run the GA for each data file
// the results are saved into the log_dir directory.
// GlobalLogFile contents a resume of each run of the GA for each file with
//               each sets of parameters
// dump files contents the minimal, average and xaximum fitness for each
//      generation in a run
// log files gives more détails for each generation and the différents weight
//     for the facts evety 100 generations

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include "data.h"
#include "ga.h"


/* GLOBALES */
int Verbose = 0;
int UltraVerbose = 0;
int Dump = 1;
FILE * DumpFile = NULL;
int Log = 1;
FILE * LogFile = NULL;
int GlobalLog = 1;
FILE * GlobalLogFile = NULL;
//const char * games_dir = "data/";
//const char * log_dir = "log/";
const char * games_dir = "/Users/dex/FAC/M1_MEMOIRE/gaGGPfinalState/gaGGPfinalState/data/";
const char * log_dir = "/Users/dex/FAC/M1_MEMOIRE/gaGGPfinalState/gaGGPfinalState/log/";


void createGlobalLogFile(const char * dir);
void createDumpFile(const char * dir, const char * gameName,
                    double c,double m, int p);
void createLogFile(const char * dir, const char * gameName,
                   double c, double m, int p);


int main(int argc, const char * argv[]) {
  testBase * data;

  /* Goldberg p71 :
   * pmutaion = 0.0333 (inversely proportional to the population size)
   * pcross = 0.6
   * popsize = 30
   */
  double pcross[] = {0.1, 0.3, 0.5, 0.7};   // crossover
  int numCross = 4;
  double pmut[] = {0.1, 0.2, 0.3, 0.4};   // mutation
  int numMut = 4;
  double numIndiv = 100; // number of individuals
  double numGen = 1001;  // number of generations

  // to identify the current test bench
  printf("%d\n", getpid());

  // opens games directory and run the genetic algorithm on each
  struct dirent *fileInfo;
  DIR *rep;
  rep = opendir(games_dir);
  
  // open global log file
  if (GlobalLog) createGlobalLogFile(log_dir);
  
  while ((fileInfo = readdir(rep))) {
    // not a file name
    if (fileInfo->d_name[0] == '.') continue;

    // load game file data
    char gamePath[512];
    sprintf(gamePath, "%s%s", games_dir, fileInfo->d_name);
    data = loadData(gamePath);

    // cut file extension
    *(strrchr (fileInfo->d_name, '.')) = '\0';
    printf("%s loaded\n", fileInfo->d_name);

    // print game data, list of differents facts and all finals states
    if (UltraVerbose) dumpFacts(data);
    if (UltraVerbose) dumpFinalStates(data);

    // skip a few lines before the new game
    if (GlobalLog) fprintf(GlobalLogFile,  "\n\n");
              
    int i, j, player;
    for (player = 0; player < data->iNumPlayers; player++) {
      for (i = 0 ; i < numMut ; i++) {
        for (j = 0 ; j < numCross ; j++) {
          printf("%d ", i * numCross + j);
          
          fflush(stdout);

          // open fitness dump file
          if (Dump) createDumpFile(log_dir, fileInfo->d_name, pcross[j], pmut[i], player);
          // open log file
          if (Log) createLogFile(log_dir, fileInfo->d_name, pcross[j], pmut[i], player);
          // write file name in global log file
          if (GlobalLog)
            fprintf(GlobalLogFile,  "--------------- %s - pcross = %4.2f - "
                                    "pmut = %4.2f - eval score %d ---------------\n",
                                    fileInfo->d_name, pcross[j], pmut[i], player);
          // run the genetic algorithm
          ga(data, pcross[j], pmut[i], player, numIndiv, numGen);
          // close dump file
          if (Dump) fclose(DumpFile);
          // close log file
          if (Log) fclose(LogFile);
          fflush(GlobalLogFile);
        }
      }
      if (GlobalLog) fprintf(GlobalLogFile, "\n");
      puts("done.");
    }
    
    // free game data
    freeData(data);
  }
  // close global log file
  if (Log) fclose(GlobalLogFile);
  closedir(rep);

  return 0;
}

void createDumpFile(const char * dir, const char * gameName, double c, double m, int p){
    char filePath[512];
    sprintf(filePath, "%s%d_%s_dump_c%4.2f_m%4.2f_p%d.txt",
            dir, getpid(), gameName, c, m, p);
    if ((DumpFile = fopen(filePath, "w")) == NULL) {
      fprintf(stderr, "error in main() : can't open dump file\n");
      exit(EXIT_FAILURE);
    }
    fprintf(DumpFile, "%s\n--------------- fitness ---------------\n"
                      "max\t\tavg\t\tmin\n\n", gameName);
}

void createLogFile(const char * dir, const char * gameName, double c, double m, int p){
  char filePath[512];
  sprintf(filePath, "%s%d_%s_log_c%4.2f_m%4.2f_P%d.txt",
          dir, getpid(), gameName, c, m, p);
  if ((LogFile = fopen(filePath, "w")) == NULL) {
    fprintf(stderr, "error in main() : can't open log file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(LogFile, "--------------- %s ---------------\n", gameName);
}

void createGlobalLogFile(const char * dir) {
  char filePath[512];
  sprintf(filePath, "%s%d__GLOBAL-LOG.txt", dir, getpid());
  if ((GlobalLogFile = fopen(filePath, "w")) == NULL) {
    fprintf(stderr, "error in main() : can't open global log file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(GlobalLogFile, "------------ test bench %d ------------\n", getpid());
}

