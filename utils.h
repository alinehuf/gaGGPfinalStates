//
//  utils.h
//  gaGGPfinalState
//
//  Created by AlineHUF 05/04/2013
//
// some tools to generate random numbers and the sigmoid function

#ifndef gaGGPfinalState_utils_h
#define gaGGPfinalState_utils_h

#include <math.h>

/* give the result of a sigmoid fonction applied on value
 */
//inline
double sigmoid(double value) {
  const double e = 2.71828182845904523536;
  // -value * 3 * e => x E [-1, 1] => y E [0, 1]
  // else with -value => x E [-1, 1] => y E [0.275, 0.725]
  return 1 / (1 + exp(-value * 3 * e));
}

/* generates a random double between vmin and vmax (not included)
 * requires that srand(SEED) is called before
 */
//inline
double randDouble(double vmin, double vmax) {
  return (double) rand() / RAND_MAX * (vmax-vmin) + vmin;
}

/* generates a random int between vmin and vmax (not included)
 * requires that srand(SEED) is called before
 */
//inline
int randInt(int vmin, int vmax) {
  if (vmax == vmin) return vmin;
  return rand() % (vmax-vmin) + vmin;
}


#endif
