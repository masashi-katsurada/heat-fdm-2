/*
 * heat2d-disk-i.c -- solve the heat equation in a two dimensional disk
 * by semi-implicit finite difference method.
 *
 *      To compile this program on WS's in 6701,
 *              gcc -o heat2d-disk-i heat2d-disk-i.c ptrilu.o call_gnuplot.o -lmatrix -lglscd -lX11 -lm
 *      or try
 *              ccmg heat2d-disk-i.c ptrilu.o call_gnuplot.o
 */

#include <stdio.h>
#include <math.h>

/* to use matrix, new_matrix() */
#include "matrix.h"
/* to use gnuplot */
#include "call_gnuplot.h"
/* to use ptrilu() */
#include "ptrilu.h"

double u0(double, double);
double maxnorm(int, int, matrix);
double exactu(double, double, double);

double pi;

int main()
{
  int Nr, Nt, i, j, n, skip, nMax;
  double hr, ht, ri, theta_j;
  double lambda_r, lambda_t;
  double tau, tau_limit, tau_limit2, Tmax, t, dt, s, ex, M;
  matrix u, newu;
  vector *al, *ad, *au, *ab, *ar, b;
  char label[100];
  /* 次の変数を 0 にすると原点だけで誤差を計算する。
   * 0 以外の値 (例えば 1) にすると円盤全体で誤差を計算する。
   */
  int zentai = 1;

  pi = 4.0 * atan(1.0);

  /* 分割の細かさを指定する */
  printf("Nr, Nt: ");
  scanf("%d %d", &Nr, &Nt);
  /* 差分解を記憶する変数の準備 */
  if ((u = new_matrix(Nr + 1, Nt + 1)) == NULL) {
    fprintf(stderr, "配列 u を確保できませんでした。");
    exit(1);
  }
  if ((newu = new_matrix(Nr + 1, Nt + 1)) == NULL) {
    fprintf(stderr, "配列 newu を確保できませんでした。");
    exit(1);
  }
  /* 連立1次方程式の係数行列を記憶するための変数の準備 */
  al = malloc(sizeof(vector *) * Nr);
  ad = malloc(sizeof(vector *) * Nr);
  au = malloc(sizeof(vector *) * Nr);
  ab = malloc(sizeof(vector *) * Nr);
  ar = malloc(sizeof(vector *) * Nr);
  for (i = 1; i < Nr; i++) {
    if ((al[i] = new_vector(Nt)) == NULL)
      perror("cannot allocate matirx.");
    if ((ad[i] = new_vector(Nt)) == NULL)
      perror("cannot allocate matirx.");
    if ((au[i] = new_vector(Nt)) == NULL)
      perror("cannot allocate matirx.");
    if ((ab[i] = new_vector(Nt)) == NULL)
      perror("cannot allocate matirx.");
    if ((ar[i] = new_vector(Nt)) == NULL)
      perror("cannot allocate matirx.");
  }
  if ((b = new_vector(Nt + 1)) == NULL)
    perror("cannot allocate matirx.");

  printf("Tmax: ");
  scanf("%lf", &Tmax);

  hr = 1.0 / Nr;
  ht = 2 * pi / Nt;

  tau_limit = 0.5 * (hr * hr * ht * ht) / (1 + ht * ht);
  tau_limit2 = 0.25 * hr * hr;
  printf("τ(陽解法の場合の上限=%g, λr≦1/4のための上限=%g): ",
         tau_limit, tau_limit2);
  scanf("%lf", &tau);

  lambda_r = tau / (hr * hr);
  lambda_t = tau / (ht * ht);

  printf("Δt: ");
  scanf("%lf", &dt);
  if (dt < tau) {
    dt = tau;
    printf("Δt=%g\n", dt);
  }
  skip = rint(dt / tau);

  open_gnuplot();

  for (i = 0; i <= Nr; i++) {
    ri = i * hr;
    for (j = 0; j <= Nt; j++) {
      theta_j = j * ht;
      u[i][j] = u0(ri, theta_j);
    }
  }
  /* 係数行列に値をセットして、LU 分解しておく */
  for (i = 1; i < Nr; i++) {
    double ri, ri2, d, od;
    ri = i * hr;
    ri2 = ri * ri;
    d = 1 + 2 * lambda_t / ri2;
    od = -lambda_t / ri2;
    for (j = 0; j < Nt; j++) {
      ad[i][j] = d;
      al[i][j] = au[i][j] = od;
      ab[i][j] = ar[i][j] = 0.0;
    }
    ab[i][0] = ab[i][Nt - 2] = ar[i][0] = ar[i][Nt - 2] = od;
    ptrilu0(Nt, al[i], ad[i], au[i], ab[i], ar[i]);
  }

  /* 時間に関するループ */
  nMax = rint(Tmax / tau);
  for (n = 1; n <= nMax; n++) {
    for (i = 1; i < Nr; i++) {
      double alpha, beta;
      /* 右辺を作る */
      alpha = lambda_r * (1.0 + 0.5 / i);
      beta = lambda_r * (1.0 - 0.5 / i);
      for (j = 0; j < Nt; j++) {
        b[j] = (1 - 2 * lambda_r) * u[i][j] +
          alpha * u[i + 1][j] + beta * u[i - 1][j];
      }
      /* 連立1次方程式を解く */
      ptrisol0(Nt, al[i], ad[i], au[i], ab[i], ar[i], b);
      /* コピーする */
      for (j = 0; j < Nt; j++)
        newu[i][j] = b[j];
      newu[i][Nt] = b[0];
    }
    /* 同次 Dirichlet 境界条件 */
    for (j = 0; j <= Nt; j++)
      newu[Nr][j] = 0.0;

    /* 原点 */
    s = 0.0;
    for (j = 0; j < Nt; j++)
      s += u[1][j];

    newu[0][0] = 4 * lambda_r * (s / Nt - u[0][0]) + u[0][0];
    for (j = 1; j <= Nt; j++)
      newu[0][j] = newu[0][0];

    /* 更新 */
    for (i = 0; i <= Nr; i++)
      for (j = 0; j <= Nt; j++)
        u[i][j] = newu[i][j];

    t = n * tau;
    /* 誤差を測る */
    if (zentai) {
      M = 0.0;
      for (i = 0; i <= Nr; i++) {
        ri = i * hr;
        for (j = 0; j <= Nt; j++) {
          double e;
          theta_j = j * ht;
          e = fabs(exactu(ri, theta_j, t) - u[i][j]);
          if (e > M)
            M = e;
        }
      }
      printf("n=%d, norm=%g, t=%g, 誤差=%g\n",
             n, maxnorm(Nr + 1, Nt + 1, u), t, M);
    }
    else {
      ex = exactu(0.0, 0.0, t);
      M = fabs(u[0][0] - ex);
      printf("n=%d, norm=%g, t=%g, exactu=%g , 誤差=%g\n",
             n, maxnorm(Nr + 1, Nt + 1, u), t, ex, M);
    }

    if (n % skip == 0) {
      t = n * tau;
      sprintf(label, "t=%g", t);
      disk(Nr, Nt, u, label);
    }
  }

  close_gnuplot();

  return 0;
}

#define mu01 (2.404825557695773)

/* 初期値 */
double u0(double r, double theta)
{
  return j0(mu01 * r);
}

/* 厳密解 */
double exactu(double r, double theta, double t)
{
  return exp(- mu01 * mu01 * t) * j0(mu01 * r);
}

double maxnorm(int m, int n, matrix u)
{
  int i, j, i0, j0;
  double tmpmax, absu;
  i0 = 0;
  j0 = 0;
  tmpmax = fabs(u[0][0]);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      if ((absu = fabs(u[i][j])) > tmpmax) {
        tmpmax = absu;
        i0 = i;
        j0 = j;
      }
  printf("(i,j)=(%d,%d)", i0, j0);
  return tmpmax;
}
