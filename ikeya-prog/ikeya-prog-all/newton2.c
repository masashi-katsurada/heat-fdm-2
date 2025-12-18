/*
 * newton.c -- Newton 法で方程式 f(x)=0 を解く
 *  コンパイル gcc -o newton newton.c -lm
 *  いずれも実行可能ファイルの名前は newton で、実行は ./newton
 */

#include <stdio.h>
#include <math.h>
#include<stdlib.h>

double r;

int main ()
{
  int i, maxitr = 100;
  double f(double), dfdx(double), x, dx, eps;

  printf(" 初期値 x0, 既知定数γ, 許容精度ε=");
  scanf("%lf%lf%lf", &x, &r, &eps);

  for (i = 0; i < maxitr; i++) {
    dx = - f(x) / dfdx(x);
    x += dx;
    printf("f(%20.15f)=%9.2e\n", x, f(x));
    if (fabs(dx) <= eps)
      break;
  }
  printf("λ=%20.15f\n)", x*x);
 
  return 0;
}

double f(double x)
{
  return tan(x) - (2*x*r)/(x*x - r*r);
}

/* 関数 f の導関数 (df/dx のつもりで名前をつけた) */
double dfdx(double x)
{
  return 1/(cos(x)*cos(x)) + (2*r*r*r + 2*x*x*r)/((x*x - r*r)*(x*x - r*r));
}
