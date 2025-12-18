/*
 * heat2d-disk-adi.c --- ADI法により二次元円盤の熱方程式を解くプログラム
 *
 *  参考: http://www.math.meiji.ac.jp/~mk/labo/text/heat-fdm-2.pdf
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ptrilu.h"
#include "trid-lu.h"
/* to use gnuplot */
#define USEGNUPLOT
#ifdef USEGNUPLOT
#include "call_gnuplot.h"
#endif

#define BESSEL

/* 初期値 */
double u0(double, double);
/* 境界値 */
double boundary_data(double, double, double);
/* 厳密解 */
double exactu(double r, double theta, double t);
/* ノルム計算 */
double maxnorm(int m, int n, matrix u);
/* 円周率 */
double PI;

int main()
{
  int n, i, j, nMax, NN, Nr, Nt;
  double tau, t, Tmax, hr, ht, dt;
  double ri, rr, coef_a, coef_b, coef_c, lambda_r, lambda_t, s;
  vector b, Bd, Bl, Bu;
  matrix Ad, Al, Au, Ab, Ar, U, UU;
  int zentai = 1, output = 0;
  double MaxError;
  double tau_limit,tau_limit2;
  int skip;
  char label[100];

  PI = 4.0 * atan(1.0);

  /* いつまで計算するか */
  printf("Tmax: "); scanf("%lf", &Tmax);

  /* r, θ各方向の分割数 */
  printf("Nr, Nθ: ");
  scanf("%d%d", &Nr, &Nt);

  /* 刻み幅 */
  hr = 1.0 / Nr;
  ht = 2 * PI / Nt;

  tau_limit = 0.5 * (hr * hr * ht * ht) / (1 + ht * ht);
  tau_limit2 = 0.25 * hr * hr;
  printf("τ(陽解法の場合の上限=%g, λr≦1/4のための上限=%g): ",
         tau_limit, tau_limit2);
  printf("τ: "); scanf("%lf", &tau);

  printf("結果を出力する時間間隔: "); scanf("%lf", &dt);
  skip = rint(dt / tau);

  /* λr,λθ */
  lambda_r = tau / (2 * hr * hr);
  lambda_t = tau / (2 * ht * ht);

  /* メモリの確保 */

  /* 連立1次方程式を解くための作業ベクトル */
  NN = (Nt > Nr) ? Nt : Nr;
  b = new_vector(NN + 1);
  /* 差分解 */
  U = new_matrix(Nr + 1, Nt + 1);
 /* 差分解(作業用) */
  UU = new_matrix(Nr + 1, Nt + 1);
  /* 連立方程式の係数行列 B */
  Bl = new_vector(Nr + 1);
  Bd = new_vector(Nr + 1);
  Bu = new_vector(Nr + 1);
  /* 連立方程式の行列 Ai */
  Al = new_matrix(Nr + 1, Nt+1);
  Ad = new_matrix(Nr + 1, Nt+1);
  Au = new_matrix(Nr + 1, Nt+1);
  Ar = new_matrix(Nr + 1, Nt+1);
  Ab = new_matrix(Nr + 1, Nt+1);
  /* 確保できたか確認 */
  if ((b == NULL) || (U == NULL) || (UU == NULL)
      || (Bd == NULL) || (Bl == NULL) || (Bu == NULL)
      || (Ad == NULL) || (Al == NULL) || (Au == NULL)
      || (Ar == NULL) || (Ab == NULL)) {
    fprintf(stderr, "memory allocation error\n");
    exit(1);
  }

  /* r方向の係数行列 B を作る */
  for (i = 1; i < Nr; i++) {
    Bd[i] = 1.0 + 2 * lambda_r;
    Bl[i] = - lambda_r * (1.0 - 1.0 / (2 * i));
    Bu[i] = - lambda_r * (1.0 + 1.0 / (2 * i));
  }
  /* 係数行列 B をLU分解する */
  trilu(Nr - 1, Bl + 1, Bd + 1, Bu + 1);

  /* θ方向の係数行列 Ai を作り、LU分解する */
  for (i = 1; i < Nr; i++) {
    ri = i * hr;
    rr = ri * ri;
    coef_a = 1 + 2 * lambda_t / rr;
    coef_b = - lambda_t / rr;
    for (j = 0; j < Nt; j++) {
      Al[i][j] = coef_b; Ad[i][j] = coef_a; Au[i][j] = coef_b;
      Ar[i][j] = 0;
      Ab[i][j] = 0;
    }
    Ar[i][0] = coef_b;
    Ab[i][0] = coef_b;
    ptrilu0(Nt, Al[i], Ad[i], Au[i], Ab[i], Ar[i]);
  }

  /* 初期値 */
  for (i = 0; i <= Nr; i++)
    for (j = 0; j <= Nt; j++)
      U[i][j] = u0(i * hr, j * ht);

  /* GNUPLOT を起動 */
#ifdef USEGNUPLOT
  open_gnuplot();
#endif

  /* 時間サイクルをまわす */
  nMax = rint(Tmax / tau);
  for (n = 0; n < nMax; n++) {

    /* 第n段から第n+1/2段 */
    /* 原点での値を計算 (陽解法) */
    s = 0;
    for (j = 0; j < Nt; j++)
      s += U[1][j];
    UU[0][0] = 4 * lambda_r * (s / Nt - U[0][0]) + U[0][0];
    for (j = 1; j <= Nt; j++)
      UU[0][j] = UU[0][0];
    /* r方向を陰的に解く */
    for (j = 0; j < Nt; j++) {
      /*境界の代入 */
      UU[Nr][j] = boundary_data(1.0, j * ht, (n + 0.5) * tau);
      /* 差分方程式の右辺既知項を計算 (係数を事前に計算しておくと良い？) */
      for (i = 1; i < Nr; i++) {
        ri = i * hr;
        rr = ri * ri;
        coef_a = 1.0 - 2 * lambda_t / rr;
        coef_b = lambda_t / rr;
        if (j != 0)
          b[i] = coef_a * U[i][j] + coef_b * (U[i][j - 1] + U[i][j + 1]);
        else
          b[i] = coef_a * U[i][j] + coef_b * (U[i][Nt - 1] + U[i][j + 1]);
      }
      /* 右辺に移項する処理 */
      b[1] += lambda_r * (1 - 1.0 / 2) * UU[0][j];
      b[Nr - 1] += lambda_r * (1 + 1.0 / (2 * (Nr - 1))) * UU[Nr][j];
      /* r方向方程式を解く */
      trisol(Nr - 1, Bl + 1, Bd + 1, Bu + 1, b + 1);
      /* UU に計算結果を残す */
      for (i = 1; i < Nr; i++)
        UU[i][j] = b[i];
    }
    /* j=Nθにも値を入れておく (必要ある？) */
    for (i = 0; i <= Nr; i++)
      UU[i][Nt] = UU[i][0];
    
    /* 第n+1/2段から第n+1段 */
    /* θ方向を陰的に解く */
    for (i = 1; i < Nr; i++) {
      /* 差分方程式の右辺既知項を計算 (係数を事前に計算しておくと良い？) */
      ri = i * hr;
      rr = hr / (2 * ri);
      coef_a = 1 - 2 * lambda_r;
      coef_b = lambda_r * (1 + rr);
      coef_c = lambda_r * (1 - rr);
      for (j = 0; j < Nt; j++)
        b[j] = coef_a * UU[i][j] + coef_b * UU[i + 1][j] + coef_c * UU[i - 1][j];

      /*θ方向方程式を解く*/
      ptrisol0(Nt, Al[i], Ad[i], Au[i], Ab[i], Ar[i], b);
      /* U に計算結果を残す */
      for (j = 0; j < Nt; j++)
        U[i][j] = b[j];
      /* j=Nθにも値を入れておく (必要ある？) */
      U[i][Nt] = U[i][0];
    }
    /* 原点での値を計算 (陽解法) */
    s = 0;
    for (j = 0; j < Nt; j++)
      s += UU[1][j];
    U[0][0] = 4 * lambda_r * (s / Nt - UU[0][0]) + UU[0][0];
    for (j = 1; j <= Nt; j++)
      U[0][j] = U[0][0];
    /* 境界データの代入 */
    t = (n + 1) * tau;
    for (j = 0; j <= Nt; j++)
      U[Nr][j] = boundary_data(1.0, j * ht, t);

    if ((n + 1) % skip == 0) {
      /* 計算結果 (数値データ) の出力 */
      if (output)
        for (i = 0; i <= Nr; i++) {
          for (j = 0; j <= Nt; j++) {
            printf("%2d ", j);
            printf("%.3f,", U[i][j]);
          }
          printf("\n");
        }
      sprintf(label, "t=%g", t);
#ifdef USEGNUPLOT
      disk(Nr, Nt, U, label);
#endif
      /* 誤差を測る */
      if (zentai) {
        MaxError = 0.0;
        for (i = 0; i <= Nr; i++) {
          ri = i * hr;
          for (j = 0; j <= Nt; j++) {
            double e;
            e = fabs(exactu(ri, j * ht, t) - U[i][j]);
            if (e > MaxError)
              MaxError = e;
          }
        }
        printf("n=%d, norm=%e, t=%f, 誤差=%e\n",
               n + 1, maxnorm(Nr + 1, Nt + 1, U), t, MaxError);
      }
      else {
        double ex = exactu(0.0, 0.0, t);
        MaxError = fabs(U[0][0] - ex);
        printf("n=%d, norm=%g, t=%g, exactu=%g , 誤差=%g\n",
               n, maxnorm(Nr + 1, Nt + 1, U), t, ex, MaxError);
      }
    }
  }
  return 0;
}

#ifdef BESSEL
/* 同次 Dirichlet 境界条件で簡単な厳密解 (Bessel関数を利用) */
#define mu01    (2.404825557695773)

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

/* 境界値 */
double boundary_data(double r, double theta, double t)
{
  /* 同次 Dirichlet 境界条件 */
  return 0.0;
}
#else
/* 非同次 Dirichlet 境界条件で簡単な厳密解 … 定数解！ */
double u0(double r, double theta)
{
  return 1.0;
}

double exactu(double r, double theta, double t)
{
  return 1.0;
}

double boundary_data(double r, double theta, double t)
{
  return 1.0;
}
#endif

double maxnorm(int m, int n, matrix u)
{
  int i, j, ii, jj;
  double tmpmax, absu;
  ii = 0;
  jj = 0;
  tmpmax = fabs(u[0][0]);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      if ((absu = fabs(u[i][j])) > tmpmax) {
        tmpmax = absu;
        ii = i;
        jj = j;
      }
  printf("(i,j)=(%d,%d)", ii, jj);
  return tmpmax;
}
