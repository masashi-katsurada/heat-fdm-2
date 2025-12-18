/*
 * heat2d-disk-e.c --- 円盤における熱方程式を陽解法で解く
 *
 *      To compile this program on WS's in 6701,
 *              gcc -c call_gnuplot.c
 *              ccmg heat2d-disk-e2.c call_gnuplot.o
 *      You can get call_gnuplot.c and call_gnuplot.h on the page
 *              http://www.math.meiji.ac.jp/~mk/program/
 *
 *      Sample input

oyabun% ./heat2d-e2
Nr, Nt: 20 80
Tmax: 1
τ(≦7.66336e-06): 5e-6
λ=0.00281057 になりました。
Δt(≧5e-06): 5e-3

 */

#include <stdio.h>
#include <math.h>

/* to use matrix, new_matrix() */
#include <matrix.h>

#include "call_gnuplot.h"

double u0(double, double);
double exactu(double, double, double);
double maxnorm(int, int, matrix);

int main()
{
  int Nr, Nt, i, j, n, skip, nMax, id;
  double pi, hr, ht, ri, ri2, theta_j, M, ex;
  double lambda_r, lambda_t, lambda;
  double tau, Tmax, t, dt, s;
  matrix u, newu;
  char label[100], fname[100];
  /* 次の変数を 0 にすると原点だけで誤差を計算する。
   * 0 以外の値 (例えば 1) にすると円盤全体で誤差を計算する。
   */
  int zentai = 1;

  pi = 4.0 * atan(1.0);
  /* 分割数を決定 */
  printf("Nr, Nt: ");
  scanf("%d %d", &Nr, &Nt);
  if ((u = new_matrix(Nr + 1, Nt + 1)) == NULL) {
    fprintf(stderr, "配列 u を確保できませんでした。");
    exit(1);
  }
  if ((newu = new_matrix(Nr + 1, Nt + 1)) == NULL) {
    fprintf(stderr, "配列 newu を確保できませんでした。");
    exit(1);
  }
  /* 最終時刻の決定 */
  printf("Tmax: ");
  scanf("%lf", &Tmax);

  /* 刻み幅 */
  hr = 1.0 / Nr;
  ht = 2 * pi / Nt;

  /* 時間刻み幅の決定 */
  printf("τ(≦%g): ",
         0.5 * (hr * hr * ht * ht) / (1 + ht * ht));
  scanf("%lf", &tau);

  /* λr, λθ */
  lambda_r = tau / (hr * hr);
  lambda_t = tau / (ht * ht);
  lambda = lambda_r + lambda_t;

  printf("λ=%g になりました。\n", lambda);

  /* 結果を出力する時間間隔を決定 */
  printf("Δt(≧%g): ", tau);
  scanf("%lf", &dt);
  if (dt < tau) {
    dt = tau;
    printf("Δ=%g\n", dt);
  }
  skip = rint(dt / tau);

  /* GNUPLOT の準備 */
  open_gnuplot2();

  /* 初期値の設定 */
  for (i = 0; i <= Nr; i++) {
    ri = i * hr;
    for (j = 0; j <= Nt; j++) {
      theta_j = j * ht;
      u[i][j] = u0(ri, theta_j);
    }
  }
  id = 0;
  sprintf(label, "t=%g", 0.0);
  sprintf(fname, "disk%05d.jpg", id++);
  disk2(Nr, Nt, u, label, fname);

  /* ループ */
  nMax = rint(Tmax / tau);
  for (n = 1; n <= nMax; n++) {
    for (i = 1; i < Nr; i++) {
      ri = i * hr;
      ri2 = ri * ri;
      for (j = 0; j < Nt; j++) {
        int jm1 = j - 1;
        if (jm1 == -1)
          jm1 = Nt - 1;
        newu[i][j] = (1 - 2 * lambda_r - 2 * lambda_t / (ri2)) * u[i][j]
          + lambda_r * (u[i + 1][j] + u[i - 1][j])
          + lambda_r * hr / (2 * ri) * (u[i + 1][j] - u[i - 1][j])
          + lambda_t / (ri2) * (u[i][j + 1] + u[i][jm1]);
      }
    }
    for (i = 1; i < Nr; i++)
      newu[i][Nt] = newu[i][0];
    /* Dirichlet 境界条件 */
    for (j = 0; j <= Nt; j++)
      newu[Nr][j] = 0.0;

    /* 原点での値 */
    s = 0;
    for (j = 0; j < Nt; j++)
      s += u[1][j];

    newu[0][0] = lambda_r * (4 * (s / Nt - u[0][0])) + u[0][0];
    for (j = 1; j <= Nt; j++)
      newu[0][j] = newu[0][0];

    /* 値の更新 */
    for (i = 0; i <= Nr; i++)
      for (j = 0; j <= Nt; j++)
        u[i][j] = newu[i][j];

    t = tau * n;

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

    /* Δtの整数倍の時刻ではグラフを描く */
    if (n % skip == 0) {
      t = tau * n;
      sprintf(label, "t=%g", t);
      sprintf(fname, "disk%05d.jpg", id++);
      disk2(Nr, Nt, u, label, fname);
    }
  }

  close_gnuplot();

  return 0;
}

/* 厳密解を計算する関数 */

#define mu01    (2.404825557695773)

double exactu(double r, double theta, double t)
{
  return exp(- mu01 * mu01 * t) * j0(mu01 * r);
}

/* 初期データ */
double u0(double r, double theta)
{
  return j0(mu01 * r);
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
