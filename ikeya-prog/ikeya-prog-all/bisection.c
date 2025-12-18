/*
 * bisection.c -- 二分法 (bisection method) で方程式 f(x)=0 を解く
 *   コンパイル: gcc -o bisection bisection.c -lm
 *   実行: ./bisection
 *   入力は "0 1 1e-15" など。ここで 1e-15 は 1かける 10 の -15 乗を意味する。
 */

#include <stdio.h>
#include <math.h>

double f(double x)
{
  return cos(x) - x;
}

int main()
{
  int i, maxitr = 100;
  double a, b, c, eps;
  double fa, fb, fc;

  printf(" 探す区間の左端, 右端, 許容精度ε: ");
  scanf("%lf %lf %lf", &a, &b, &eps);

  fa = f(a);
  fb = f(b);
  if (fa * fb > 0.0) {
    printf(" f(α) f(β) > 0 なのであきらめます。\n");
    exit(0);
  }
  else {
    printf("f(%20.15f)=%9.2e, f(%20.15f)=%9.2e\n", a, fa, b, fb);
    for (i = 0; i < maxitr; i++) {
      c = (a + b) / 2;
      fc = f(c);
      if (fc == 0.0)
	break;
      else if (fa * fc <= 0.0) {
	/* 左側 [a,c] に根がある */
	b = c;
	fb = fc;
      }
      else {
	/* 左側 [a,c] には根がないかもしれない。[c,b] にあるはず */
	a = c;
	fa = fc;
      }
      printf ("f(%20.15f)=%9.2e, f(%20.15f)=%9.2e\n", a, fa, b, fb);
      if ((b - a) <= eps)
	break;
    }
    printf ("二分法による近似解=%20.15f\n", c);
  }
  return 0;
}
