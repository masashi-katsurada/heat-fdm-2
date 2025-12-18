/*
 * newton.c -- Newton 法で方程式 f(x)=0 を解く
 *  コンパイル gcc -o newton newton.c -lm
 *  いずれも実行可能ファイルの名前は newton で、実行は ./newton
 */

#include <stdio.h>
#include <math.h>

int main ()
{
  int i, maxitr = 100;
  double f(double), dfdx(double), x, dx, eps;

  printf(" 初期値 x0, 許容精度ε=");
  scanf("%lf%lf", &x, &eps);

  for (i = 0; i < maxitr; i++) {
    dx = - f(x) / dfdx(x);
    x += dx;
    printf("f(%20.15f)=%9.2e\n", x, f(x));
    if (fabs(dx) <= eps)
      break;
  }
  return 0;
}

double f(double x)
{
  return cos(x) - x;
}

/* 関数 f の導関数 (df/dx のつもりで名前をつけた) */
double dfdx(double x)
{
  return - sin(x) - 1.0;
}
