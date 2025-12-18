/*
 * find_Bessel_roots.c -- 二分法 (bisection method) で方程式 Jn(x)=0 を解く
 *
 *   正の解を小さい方から指定した個数だけ求める
 *   コンパイルは gcc -o find_Bessel_roots find_Bessel_roots.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print(int N, int n, double *a)
{
  int i;
  printf("  k           μ%d,k                 J%d(μ%dk)\n", n, n, n);
  for (i = 1; i <= N; i++)
    printf("%3d %25.18e %15e %f\n", i, a[i], jn(n,a[i]), a[i]-a[i-1]);
}

int main(int argc, char **argv)
{
  int i, j, maxitr = 100;
  double a, b, c, eps;
  double fa, fb, fc;
  int n;
  int num_of_roots;
  double *roots;

  if (argc == 3) {
    n = atoi(argv[1]);
    num_of_roots = atoi(argv[2]);
  }
  else {
    printf("Jn の最初の N 個の零点を求めます。\n");
    printf("n="); scanf("%d", &n);
    printf("N="); scanf("%d", &num_of_roots);
  }
  if ((roots = malloc(sizeof(double) * (num_of_roots + 1))) == NULL) {
    fprintf(stderr, "cannot allocate memory\n");
    exit(1);
  }
  a = 0.1; b = 0.2; eps = 5e-16;
  roots[0] = 0;
  fa = jn(n, a); fb = jn(n, b);
  for (i = 1; i <= num_of_roots; i++) {
    while (fa * fb > 0) {
      a = b; fa = fb;
      b += 0.1; fb = jn(n, b);
    }
    for (j = 0; j <= maxitr && (b - a) / a > eps; j++) {
      c = (a + b) / 2; fc = jn(n, c);
      if (fc == 0.0)
        break;
      else if (fa * fc <= 0.0) {
        /* 左側 [a,c] に根がある */
        b = c; fb = fc;
      }
      else {
        /* 左側 [a,c] には根がないかもしれない。[c,b] にあるはず */
        a = c; fa = fc;
      }
    }
    if (j >= maxitr) {
      fprintf(stderr, "%d番目で失敗\n", i+1);
      print(i, n, roots);
      exit(1);
    }
    roots[i] = a;
    a = b; b = a + 0.1;
    fa = jn(n, a); fb = jn(n, b);
  }
  print(num_of_roots, n, roots);
  return 0;
}
