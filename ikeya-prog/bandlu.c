/*
 * bandlu.c
 * bandtest1.c を参考
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matrix.h>

void bandlu(matrix a, int n, int m)
{
  int i, j, k, mm;
  double q;
  for (i = 0; i < n; i++) {
    mm = i + m; if (mm > n-1) mm = n-1;
    for (j = i + 1; j <= mm; j++) {
        q = a[j][i] / a[i][i];
        for (k = i + 1; k <= mm; k++)
            a[j][k] -= q * a[i][k];
        a[j][i] = q; 
    }
  }
}

void bandsolve(matrix a, vector b, int n, int m)
{
  int i, j, mm;
  for (i = 0; i < n; i++) {
    mm = i + m; if (mm > n-1) mm = n-1;
    for (j = i + 1; j <= mm; j++)
      b[j] -= a[j][i] * b[i];
  }

  b[n-1] /= a[n-1][n-1];
  for (i = n - 2; i >= 0; i--) {
    mm = i + m; if (mm > n-1) mm = n-1;
    for (j = i + 1; j <= mm; j++)
      b[i] -= a[i][j] * b[j];
    b[i] /= a[i][i];
  }
}

