/*
 * testptrilu.c --- ptrilu.c のテスト
 */

#include <stdio.h>
#include <matrix.h>
#include "ptrilu.h"

int main()
{
  int i, n;
  vector al, ad, au, ab, ar, b, x;
  n = 8;
  al = new_vector(n); ad = new_vector(n); au = new_vector(n);
  ab = new_vector(n); ar = new_vector(n);  b = new_vector(n);
  x = new_vector(n);
  for (i = 0; i < n; i++) {
    al[i] = - 1; ad[i] = 4; au[i] = - 1; ar[i] = ab[i] = 0;
  }
  ar[0] = al[0]; ab[0] = au[n-1];
  for (i = 0; i < n; i++)
    x[i] = i+1;
  for (i = 0; i < n; i++) {
    if (i == 0)
      b[i] = ar[0] * x[n-1] + ad[i] * x[i] + au[i] * x[i+1];
    else if (i == n-1)
      b[i] = al[i] * x[i-1] + ad[i] * x[i] + ab[0] * x[0];
    else
      b[i] = al[i] * x[i-1] + ad[i] * x[i] + au[i] * x[i+1];
  }
  ptrilu0(n, al, ad, au, ab, ar);
  printf("L, D, U, R\n");
  for (i = 0; i < n; i++)
    if (i == 0)
      printf("        %7.4f %7.4f %7.4f\n", ad[i],au[i], ar[i]);
    else if (i < n-2)
      printf("%7.4f %7.4f %7.4f %7.4f\n",al[i],ad[i],au[i], ar[i]);
    else
      printf("%7.4f %7.4f %7.4f\n",al[i],ad[i],au[i]);
  printf("B\n");
  for (i = 0; i < n-2; i++)
    printf("%7.4f ", ab[i]);
  printf("\n");
  ptrisol0(n, al, ad, au, ab, ar, b);
  for (i = 0; i < n; i++) {
    printf("x[%d]=%g, sol[%d]=%g\n", i, x[i], i, b[i]);
  }
  return 0;
}
