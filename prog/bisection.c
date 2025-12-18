/*
 * bisection.c -- 二分法 (bisection method) で方程式 f(x)=0 を解く
 *   コンパイルは cc -o bisection bisection.c -lm
 *   東海林先生の cco コマンドがあれば cco bisection.c
 *   ccxでコンパイル
 *   いずれも bisection という実行ファイルができる。
 *   Cーc g でgoto line
 */

#include <stdio.h>
#include <math.h>

int main()
{
  int i, maxitr = 100;
  double alpha, beta, a, b, c, eps;
  double fa, fb, fc;
  double f(double);
  double j1(double x);

  printf(" 探す区間の左端α, 右端β, 許容精度ε=");
  scanf("%lf %lf %lf", &alpha, &beta, &eps);

  a = alpha; b = beta;

  fa = f(a); fb = f(b);
  if (fa * fb > 0.0) {
    printf(" f(α) f(β) > 0 なのであきらめます。\n");
    exit(0);
  }
  else {
    for (i = 0; i < maxitr; i++) {
      c = (a + b) / 2; fc = f(c);
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
      printf("f(%24.15e)=%9.2e, f(%24.15e)=%9.2e\n", a, fa, b, fb);
      if ((b - a) <= eps)
        break;
    }
    printf("f(%24.15e)=%9.2e\n", c, fc);
  }
  return 0;
}
/*    ０次のBessel関数の零点を求めたいときは j0(x)
 *    １次のBessel関数の零点を求めたいときは j1(x)
 *     とする。
 */ 
double f(double x)
{
  return j0(x);
}







