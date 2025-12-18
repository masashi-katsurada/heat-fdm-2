/*
 * heat1d-e-glsc.c -- 1次元熱伝導方程式の初期値境界値問題を陽解法で解く。
 *   コンパイル:  ccmg heat1d-e-glsc.c
 *
 * オリジナルは fplot ライブラリィを利用した
 *  http://www.math.meiji.ac.jp/%7Emk/program/fdm/heat1d-e.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G_DOUBLE
#include <glsc.h>

int main()
{
  int i, n, nMax, N;
  double tau, h, lambda, Tmax;
  double *u, *newu;
  double f(double);
  double win_width, win_height, w_margin, h_margin;
  char message[100];

  /* N, λ を入力する */
  printf("区間の分割数 N = "); scanf("%d", &N);
  printf("λ (=τ/h^2) = "); scanf("%lf", &lambda);

  /* h, τ を計算する */
  h = 1.0 / N;
  tau = lambda * h * h;
  printf("τ=%g\n", tau);

  /* 最終時刻を入力する */
  printf("最終時刻 Tmax = "); scanf("%lf", &Tmax);

  /* ベクトル u, newu を用意する */
  u = malloc(sizeof(double) * (N+1));
  newu = malloc(sizeof(double) * (N+1));
  
  /* 初期値の代入 */
  for (i = 0; i <= N; i++)
    u[i] = f(i * h);

  /* ***************** グラフィックスの準備 ***************** */
  /* メタファイル名は "HEAT",
   * ウィンドウのサイズは、
      横 win_width  + 2 * w_margin, 縦 win_height + 2 * h_margin */
  win_width = 200.0; win_height = 200.0; w_margin = 10.0; h_margin = 10.0;
  g_init("HEAT", win_width  + 2 * w_margin, win_height + 2 * h_margin);
  /* 画面とメタファイルの両方に記録する */
  g_device(G_BOTH);
  /* 座標系の定義: [-0.1,1.1]×[-0.1,1.1] という閉領域を表示する */
  g_def_scale(0,
              -0.1, 1.1, -0.1, 1.1,
              w_margin, h_margin, win_width, win_height);
  /* 線を二種類用意する */
  g_def_line(0, G_BLACK, 0, G_LINE_SOLID);
  g_def_line(1, G_BLACK, 0, G_LINE_DOTS);
  /* 表示するための文字列の属性を定義する */
  g_def_text(0, G_BLACK, 3);
  /* 定義したものを選択する */
  g_sel_scale(0); g_sel_line(0); g_sel_text(0);

  /* タイトルと入力パラメーターを表示する */
  g_text(30.0, 30.0,
         "heat equation, homogeneous Dirichlet boundary condition");
  sprintf(message, "N=%d, lambda=%g, Tmax=%g", N, lambda, Tmax);
  g_text(30.0, 60.0, message);

  /* 座標軸を表示する */
  g_sel_line(1);
  g_move(-0.1, 0.0); g_plot(1.1, 0.0);
  g_move(0.0, -0.1); g_plot(0.0, 1.1);
  g_sel_line(0);

  /* t=0 の状態を表示する */
  g_move(0.0, u[0]);
  for (i = 1; i <= N; i++)
    g_plot(i * h, u[i]);

  /* ループを何回まわるか計算する (四捨五入) */
  nMax = rint(Tmax / tau);

  /* 時間に関するステップを進めるループ */
  for (n = 0; n < nMax; n++) {
    /* 差分方程式 (n -> n+1) */
    for (i = 1; i < N; i++)
      newu[i] = (1.0 - 2 * lambda) * u[i] + lambda * (u[i+1] + u[i-1]);
    /* 計算値を更新 */
    for (i = 1; i < N; i++)
      u[i] = newu[i];
    /* Dirichlet 境界条件 */
    u[0] = u[N] = 0.0;
    /* この時刻 (t=(n+1)τ) の状態を表示する */
    g_move(0.0, u[0]);
    for (i = 1; i <= N; i++)
      g_plot(i * h, u[i]);
  }
  
  printf("終りました。X の場合はウィンドウをクリックして下さい。\n");
  g_sleep(-1.0);
  /* ウィンドウを閉じる */
  g_term();
  return 0;
}

double f(double x)
{
  if (x <= 0.5)
    return x;
  else
    return 1.0 - x;
}
