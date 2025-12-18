/*
 * heat1n-i-glsc.c -- 陰的スキームで熱方程式 (Neumann BC) を解く
 *
 *   菊地・山本, 微分方程式と計算機演習, 山海堂 (1991) 第11章
 *   p237 のプログラムを修正・拡張したもの
 *
 C     ** IMPLICIT FINITE DIFFERENCE METHOD **
 C     **         FOR HEAT EQUATION         **
 C     nfunc : 関数の番号（5種類の関数を用意している）
 C     theta : スキームのパラメーター
 C     N : 分割の数
 C     t : 時刻
 C     tau : 時間の刻み幅
 C     Tmax : 時刻の上限
 C     h : 空間の刻み幅
 C     AL,AD,AU : 係数行列
 C     u : DIMENSION FOR UNKNOWN FUNCTION
 *
 *  このプログラムは
 *    http://www.math.meiji.ac.jp/%7Emk/program/
 *  から入手可能
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G_DOUBLE
#include <glsc.h>

/* 円周率 (math.h にある M_PI の値を利用する
 * あるいは double pi; として、どこかで pi = 4.0 * atan(1.0); と代入
 * まあ Visual C++ とかでは M_PI はないという話もあります。
 */

const double pi = M_PI;

void trilu(int, double *, double *, double *);
void trisol(int, double *, double *, double *, double *);
double f(double, int);
void print100(int, double, double *);
void axis(double, double, double, double);
void print_parameters(double, double, double, double,
                      int, double, int,
                      double, double, double, double);

int main()
{
    int N, nfunc, i, n, nMax;
    double a, b, theta, h, Tmax, tau, c1, c2, c3, c4, lambda, t;
    double *AL, *AD, *AU, *ff, *u;
    double xleft, ybottom, xright, ytop;
    /* グラフを書き換える時間間隔 */
    double dt;
    /* 何ステップおきにグラフを書き換えるか */
    int skip;
    /* 以前描いたグラフを消すか？ (0=No, 1=Yes) */
    int erase_always = 0;
    /* 数値を表示するか? (0=No, 1=Yes) */
    int print_numerical_data = 0;
    /* */
    double win_width, win_height, w_margin, h_margin;

    /* 区間の左端、右端 */
    a = 0.0;
    b = 1.0;

    /* 入力 */
    printf("入力して下さい : nfunc(1..5)=");
    scanf("%d", &nfunc);
    printf("入力して下さい : θ=");
    scanf("%lf", &theta);
    printf("入力して下さい : N=");
    scanf("%d", &N);
    printf("入力して下さい : λ=");
    scanf("%lf", &lambda);

    /* ベクトルの準備 */
    if ((AL = malloc((N + 1) * sizeof(double))) == NULL) {
        fprintf(stderr, "cannot allocate AL[]\n");
        exit(1);
    }
    if ((AD = malloc((N + 1) * sizeof(double))) == NULL) {
        fprintf(stderr, "cannot allocate AD[]\n");
        exit(1);
    }
    if ((AU = malloc((N + 1) * sizeof(double))) == NULL) {
        fprintf(stderr, "cannot allocate AU[]\n");
        exit(1);
    }
    if ((ff = malloc((N + 1) * sizeof(double))) == NULL) {
        fprintf(stderr, "cannot allocate ff[]\n");
        exit(1);
    }
    if ((u = malloc((N + 1) * sizeof(double))) == NULL) {
        fprintf(stderr, "cannot allocate u[]\n");
        exit(1);
    }
    /* */
    h = (b - a) / N;
    tau = lambda * h * h;
    printf(" 時間の刻み幅τ= %g になりました。\n", tau);

    /* 入力 */
    printf("入力して下さい : 最終時刻 Tmax=");
    scanf("%lf", &Tmax);
    printf("入力して下さい : グラフ書き換え時間間隔 Δt=");
    scanf("%lf", &dt);
    skip = (int) rint(dt / tau);
    /* skip が 0 になっていないかチェック */
    if (skip == 0) {
        fprintf(stderr, "Δtが小さすぎるので、Δt=τ=%g にします。\n", tau);
        skip = 1;
        dt = tau;
    }

    /* 係数行列の準備 */
    c1 = - theta * lambda;
    c2 = 1.0 + 2.0 * theta * lambda;
    c3 = 1.0 - 2.0 * (1.0 - theta) * lambda;
    c4 = (1.0 - theta) * lambda;
    for (i = 1; i < N; i++) {
        AL[i] = c1;
        AD[i] = c2;
        AU[i] = c1;
    }
    AL[0] = 0;
    AD[0] = c2;
    AU[0] = 2 * c1;
    AL[N] = 2 * c1;
    AD[N] = c2;
    AU[N] = 0;

    /* 初期値 */
    for (i = 0; i <= N; i++)
        u[i] = f(a + i * h, nfunc);

    /* 出力（t=0） */
    /* 出力（ここでは画面に数値を表示すること）は FORTRAN と C で
       かなり異なるので、直接の翻訳は出来ない。ごたごたするので、
       print100 という関数にまとめた */
    t = 0.0;
    print100(N, t, u);

    /* グラフィックス画面の準備 -- まず画面の範囲を決めてから */
    xleft = a - (b - a) / 10;
    xright = b + (b - a) / 10;
    ybottom = - 0.1;
    ytop = 1.1;

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
  print_parameters(xleft, ybottom, xright, ytop,
                   nfunc, theta, N,
                   lambda, tau, Tmax, t);

  /* 座標軸を表示する */
  g_sel_line(1);
  g_move(-0.1, 0.0); g_plot(1.1, 0.0);
  g_move(0.0, -0.1); g_plot(0.0, 1.1);
  g_sel_line(0);

  /* t=0 の状態を表示する */
  g_move(0.0, u[0]);
  for (i = 1; i <= N; i++)
    g_plot(i * h, u[i]);

    /* 係数行列の LU 分解 */
    trilu(N + 1, AL, AD, AU);

    /* 繰り返し計算 */
    nMax = (int) rint(Tmax / tau);
    for (n = 1; n <= nMax; n++) {
        /* 連立1次方程式を作って解く */
        for (i = 1; i < N; i++)
            ff[i] = c3 * u[i] + c4 * (u[i - 1] + u[i + 1]);
        ff[0] = c3 * u[0] + 2 * c4 * u[1];
        ff[N] = c3 * u[N] + 2 * c4 * u[N - 1];
        trisol(N + 1, AL, AD, AU, ff);
        for (i = 0; i <= N; i++)
            u[i] = ff[i];

        /* 時刻 */
        t = n * tau;

        if (n % skip == 0) {
            /* 数値データの表示 */
            if (print_numerical_data)
                print100(N, t, u);

            /* 解のグラフを描く */
            if (erase_always) {
                /* これまで描いたものを消し、パラメーターを再表示する */
                g_cls();
                print_parameters(xleft, ybottom, xright, ytop,
                                 nfunc, theta, N, lambda, tau, Tmax, t);
            }
            g_move(0.0, u[0]);
            for (i = 1; i <= N; i++)
              g_plot(i * h, u[i]);
        }
    }

    printf("マウスでウィンドウをクリックして下さい。\n");
    g_sleep(-1.0);
    g_term();

    return 0;
}

/**********************************************************/

/* 三重対角行列の LU 分解 (pivoting なし) */
void trilu(int n, double *al, double *ad, double *au)
{
  int i, nm1 = n - 1;
  /* 前進消去 (forward elimination) */
  for (i = 0; i < nm1; i++) {
    al[i + 1] /= ad[i];
    ad[i + 1] -= au[i] * al[i + 1];
  }
}

/* LU 分解済みの三重対角行列を係数に持つ三項方程式を解く */
void trisol(int n, double *al, double *ad, double *au, double *b)
{
  int i, nm1 = n - 1;
  /* 前進消去 (forward elimination) */
  for (i = 0; i < nm1; i++) b[i + 1] -= b[i] * al[i + 1];
  /* 後退代入 (backward substitution) */
  b[nm1] /= ad[nm1];
  for (i = n - 2; i >= 0; i--) b[i] = (b[i] - au[i] * b[i + 1]) / ad[i];
}

/**********************************************************/

/* 初期条件を与える関数。複数登録して番号 nfunc で選択する。 */

double f(double x, int nfunc)
{
    /* f(x)=max(x,1-x) */
    if (nfunc == 1) {
        if (x <= 0.5)
            return x;
        else
            return 1.0 - x;
    }
    /* f(x)=1 */
    else if (nfunc == 2)
        return 1.0;
    else if (nfunc == 3)
        return sin(pi * x);
    /* f(x)= 変な形をしたもの(by M.Sakaue) */
    else if (nfunc == 4) {
        if (x <= 0.1)
            return 5 * x;
        else if (x <= 0.3)
            return -2 * x + 0.7;
        else if (x <= 0.5)
            return 4.5 * x - 1.25;
        else if (x <= 0.7)
            return -x + 1.5;
        else if (x <= 0.9)
            return x + 0.1;
        else
            return -10 * x + 10.0;
    }
    /* やはり非対称な形をしたもの(by M.Sakaue) */
    else if (nfunc == 5) {
        if (x <= 0.2)
            return -x + 0.8;
        else if (x <= 0.3)
            return 4 * x - 0.2;
        else if (x <= 0.4)
            return -4 * x + 2.2;
        else if (x <= 0.7)
            return -x + 1.0;
        else if (x <= 0.8)
            return 4 * x - 2.6;
        else
            return -x + 1.5;
    } else
        return 0;
}

/**********************************************************/

/* 配列 u[] の内容(u[0],u[1],...,u[n])を一行 5 個ずつ表示する。 */

void print100(int N, double t, double *u)
{
    int i;
    printf("T= %12.4e\n", t);
    for (i = 0; i < 5; i++)
        printf("  I     u(i)    ");
    printf("\n");
    for (i = 0; i <= N; i++) {
        printf("%3d%12.4e ", i, u[i]);
        if (i % 5 == 4)
            printf("\n");
    }
    printf("\n");
}

/**********************************************************/

/* ウィンドウに座標軸を表示する */
void axis(double xleft, double ybottom, double xright, double ytop)
{
    g_sel_line(1);
    g_move(xleft, 0.0); g_plot(xright, 0.0);
    g_move(0.0, ybottom); g_plot(0.0, ytop);
    g_sel_line(0);
}

/* パラメーターの値等をウィンドウに表示する */
void print_parameters(double xleft, double ybottom, double xright, double ytop,
                      int nfunc, double theta, int N,
                      double lambda, double tau, double Tmax, double t)
{
    char message[80];

    axis(xleft, ybottom, xright, ytop);

    g_text(30.0, 30.0, "heat equation, u(0)=u(1)=0");
    sprintf(message, "nfunc=%d, theta=%g, N=%d, lambda=%g, tau=%g, Tmax=%g",
            nfunc, theta, N, lambda, tau, Tmax);
    g_text(30.0, 60.0, message);
    if (t != 0.0) {
        sprintf(message, "t = %g", t);
        g_text(30.0, 90.0, message);
    }

}
