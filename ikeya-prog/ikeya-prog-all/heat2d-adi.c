/*
 * heat2d-adi.c --- 2次元長方形領域における熱伝導方程式 (Dirichlet 条件) を
 *                  ADI 法で解く。
 *
 *   コンパイルするには
 *     gcc -c trid-lu1.c
 *   として trid-lu1.o を作っておいてから
 *     gcc -o heat2d-adi heat2d-adi.c trid-lu1.o -lmatrix -lglscd -lX11 -lm
 *   または
 *     ccmg heat2d-adi.c trid-lu1.c
 */

#include <stdio.h>
#include <math.h>
#include <matrix.h>
#define G_DOUBLE
#include <glsc.h>
#include "trid-lu1.h"

double u0(double, double), alpha(double, double);
void set_boundary_data(int, int, matrix, double, double, double, double);

int main()
{
    double a, b, c, d, c1, c2, c3, c4, c5, c6, c7, c8;
    int N_x, N_y, N, i, j, k, kmax, skip;
    double h_x, h_y, lambda_x, lambda_y, lambda, tau, Tmax, dt;
    double *alx, *adx, *aux, *aly, *auy, *ady, *ff, *gg;
    matrix U, newU;
    double t;

    /* 長方形領域の定義 */
    a = 0.0; b = 1.0; c = 0.0; d = 1.0;

    /* 分割数の決定 */
    printf("Nx, Ny: "); scanf("%d %d", &N_x, &N_y);

    N = (N_x + 1) * (N_y + 1);
    h_x = (b - a) / N_x;
    h_y = (d - c) / N_y;

    /* ベクトル、行列を記憶するための変数の準備 */
    if ((alx = (double *) malloc(N_x * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate alx[]\n");
	exit(1);
    }
    if ((adx = (double *) malloc(N_x * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate adx[]\n");
	exit(1);
    }
    if ((aux = (double *) malloc(N_x * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate aux[]\n");
	exit(1);
    }
    if ((aly = (double *) malloc(N_y * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate aly[]\n");
	exit(1);
    }
    if ((ady = (double *) malloc(N_y * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate ady[]\n");
	exit(1);
    }
    if ((auy = (double *) malloc(N_y * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate auy[]\n");
	exit(1);
    }
    if ((ff = (double *) malloc(N_x * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate ff[]\n");
	exit(1);
    }
    if ((gg = (double *) malloc(N_y * sizeof(double))) == NULL) {
	fprintf(stderr, "cannot allocate gg[]\n");
	exit(1);
    }
    if ((newU = new_matrix(N_x + 1, N_y + 1)) == NULL) {
	fprintf(stderr, "数列 newU を記憶する領域の確保に失敗\n");
	exit(1);
    }
    if ((U = new_matrix(N_x + 1, N_y + 1)) == NULL) {
	fprintf(stderr, "係数行列 U を記憶する領域の確保に失敗\n");
	exit(1);
    }
    /* 時間刻の決定 */
    printf("τ(陽解法の場合の安定性のための上限=%g): ",
	   0.5 / (1 / (2 * h_x * h_x) + 1 / (2 * h_y * h_y)));
    scanf("%lf", &tau);

    /* 最終時刻の決定 */
    printf("Tmax: "); scanf("%lf", &Tmax);

    /* パラメーターの計算 */
    lambda_x = tau / (2 * h_x * h_x);
    lambda_y = tau / (2 * h_y * h_y);
    lambda = lambda_x + lambda_y;

    /* グラフを描画する時間間隔の決定 */
    printf("Δt: "); scanf("%lf", &dt);
    skip = dt / tau;

    /* グラフを描くための準備 (GLSC の初期化その他) */
    g_init("Meta", 250.0, 160.0);
    g_device(G_BOTH);
    g_def_scale(0, 0.0, 1.0, 0.0, 1.0, 30.0, 70.0, 100.0, 72.0);
    g_def_scale(4, -1.0, 1.0, -1.0, 1.0, 30.0, 30.0, 100.0, 100.0);
    g_def_line(0, G_BLACK, 0, G_LINE_SOLID);
    g_sel_scale(0);

    /* 初期値 */
    for (i = 1; i < N_x; i++) {
	for (j = 1; j < N_y; j++)
	    U[i][j] = u0(a + i * h_x, c + j * h_y);
    }
    /* 境界条件 */
    set_boundary_data(N_x, N_y, U, a, b, c, d);
    set_boundary_data(N_x, N_y, newU, a, b, c, d);

    /* 初期条件のグラフを描く */
    g_cls();
    g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
	      150.0, 100.0, U, N_x + 1, N_y + 1,
	      1, G_SIDE_NONE, 2, 1);

    /* 係数行列の計算 (その一) */
    c1 = 1.0 + 2.0 * lambda_x;
    c2 = -lambda_x;
    c3 = 1.0 - 2.0 * lambda_y;
    c4 = lambda_y;
    for (i = 1; i < N_x; i++) {
	alx[i] = c2; adx[i] = c1; aux[i] = c2;
    }

    /* 係数行列の計算 (その二) */
    c5 = 1.0 + 2.0 * lambda_y;
    c6 = -lambda_y;
    c7 = 1.0 - 2.0 * lambda_x;
    c8 = lambda_x;
    for (j = 1; j < N_y; j++) {
	aly[j] = c6; ady[j] = c5; auy[j] = c6;
    }

    kmax = Tmax / tau + 0.5;
    printf("skip = %d \n", skip);

    /* LU分解 */
    trilu1(N_y - 1, aly, ady, auy);
    trilu1(N_x - 1, alx, adx, aux);

    /* 時間に関するループ */
    for (k = 1; k <= kmax; k++) {

	/*  第k段から第k＋1/2段までの連立1次方程式 */
	for (j = 1; j < N_y; j++) {
	    double y = c + j * h_y;
	    /* 右辺の準備 */
	    for (i = 1; i < N_x; i++) {
	      ff[i] = c3 * U[i][j] + c4 * (U[i][j + 1] + U[i][j - 1]);
	    }
	    /* 非同次境界条件の参入 */
	    ff[1]     += lambda_x * alpha(a, y);
	    ff[N_x-1] += lambda_x * alpha(b, y);
	    /* 連立1次方程式を解く */
	    trisol1(N_x - 1, alx, adx, aux, ff);
	    /* 解を記憶する */
	    for (i = 1; i < N_x; i++)
		newU[i][j] = ff[i];
	}
	/* 第k+1/2段から第k+1段までの連立1次方程式 */
	for (i = 1; i < N_x; i++) {
	    double x = a + i * h_x;
	    /* 右辺の準備 */
	    for (j = 1; j < N_y; j++) {
	      gg[j] = c7 * newU[i][j] + c8 * (newU[i + 1][j] + newU[i - 1][j]);
	    }
	    /* 非同次境界条件の参入 */
	    gg[1]     += lambda_y * alpha(x, c);
	    gg[N_y-1] += lambda_y * alpha(x, d);
	    /* 連立1次方程式を解く */
	    trisol1(N_y - 1, aly, ady, auy, gg);
	    /* 解を記憶する */
	    for (j = 1; j < N_y; j++)
		U[i][j] = gg[j];
	}
	/* Δt ごとにグラフを描く */
	if (k % skip == 0) {
	    g_cls();
	    g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
		      150.0, 100.0, U, N_x + 1, N_y + 1,
		      1, G_SIDE_NONE, 2, 1);
	}
	t = tau * k;
    }
    /* マウスでクリックされるのを待つ */
    g_sleep(-1.0);
    /* ウィンドウを消す */
    g_term();

    return 0;
}

double u0(double x, double y)
{
    /* ピラミッド型の関数 */
    if (y > 0.5)
	y = 1 - y;
    if (x > 0.5)
	x = 1 - x;
    if (y < x)
	return 5 * y;
    else
	return 5 * x;
}

double alpha(double x, double y)
{
    /* 同次 Dirichlet 境界条件 */
    return x * x - y * y;
    return 0.0;
}

void set_boundary_data(int N_x, int N_y, matrix u,
		       double a, double b, double c, double d)
{
  int i, j;
  double h_x, h_y, x, y;

  h_x = (b - a) / N_x;
  h_y = (d - c) / N_y;

  for (j = 0; j <= N_y; j++) {
    y = c + j * h_y;
    u[0][j]   = alpha(a, y);
    u[N_x][j] = alpha(b, y);
  }
  for (i = 0; i <= N_x; i++) {
    x = a + i * h_x;
    u[i][0]   = alpha(x, c);
    u[i][N_y] = alpha(x, d);
  }
}
