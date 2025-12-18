/*
 * symbandlu.h
 */

#include <matrix.h>

void band(matrix at, vector f, int N, int m);
void symbandlu(matrix at, int N, int m);
void symbandsolve(matrix at, vector f, int N, int m);
