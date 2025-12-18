/*
 * heat2n-i.c --- 2次元熱方程式 (Neumann 境界条件) を陰解法で解く
 *   コンパイルするには
 *     gcc -c symbandlu.c
 *     gcc -o heat2n-i heat2n-i.c symbandlu.o -lmatrix -lglscd -lX11 -lm
 *   安直には
 *     ccmg heat2n-i.c symbandlu.c
 *
 *   FreeBSD では、事前に
 *      setenv LD_LIBRARY_PATH /usr/local/lib:/usr/X11R6/lib:/usr/lib
 *   あるいはコンパイル・オプションとして、-l何とか の前に
 *      -L/usr/local/lib -L/usr/X11R6/lib
 *   が必要になるかも
 *
 *   このプログラムについては 1998 年度卒研の学生だった深石君に感謝します。
 */

#include <stdio.h>
#include <math.h>
#include <matrix.h>
#define G_DOUBLE
#include <glsc.h>
#include "symbandlu.h"

#define phi(i,j) (j)*m+(i)

int main()
{
    double a, b, c, d;
    int N_x, N_y, m, N, i, j, p, q, L, n, nMax;
    matrix Uk, A;
    double *B, *vector_U, cond;
    int *iwork, skip;
    double h_x, h_y, lambda_x, lambda_y, lambda, lambda_limit, tau, Tmax, dt;
    double f(double, double), Phi(double, double, double);
    double theta, gamma, alpha, beta, gamma_p, alpha_p, beta_p;
    double x, y, t, t_p;

    /* 問題を考える区間 [a,b]×[c,d] */
    a = 0.0; b = 1.0; c = 0.0; d = 1.0;

    /* 区間の分割数 */
    printf("Nx, Ny: "); scanf("%d %d", &N_x, &N_y);

    m = N_x + 1;
    N = (N_x + 1) * (N_y + 1);
    /* 空間の刻み幅 */
    h_x = (b - a) / N_x;
    h_y = (d - c) / N_y;

    /* 行列、ベクトルを記憶する変数のメモリー割り当て */
    if ((Uk = new_matrix(N_x + 1, N_y + 1)) == NULL) {
        fprintf(stderr, "数列 U^k を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((A = new_matrix(N, m + 1)) == NULL) {
        fprintf(stderr, "係数行列 A を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((B = (double *) malloc(sizeof(double) * N)) == NULL) {
        fprintf(stderr, "B を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((vector_U = (double *) malloc(sizeof(double) * N)) == NULL) {
        fprintf(stderr, "vector_U を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((iwork = (int *) malloc(sizeof(int) * N)) == NULL) {
        fprintf(stderr, "iwork を記憶する領域の確保に失敗\n");
        exit(1);
    }

    /* θ法の重みの決定 */
    printf("θ (0≦θ≦1): "); scanf("%lf", &theta);

    if (theta == 1.0) {
        printf("τ: "); scanf("%lf", &tau);
    } else {
        printf("τ(≦%g≡最大値ノルムに関する安定性条件を満たすτの上限): ",
               0.5 / (1 - theta) / (1 / (h_x * h_x) + 1 / (h_y * h_y)));
        scanf("%lf", &tau);
    }

    lambda_x = tau / (h_x * h_x);
    lambda_y = tau / (h_y * h_y);
    lambda = lambda_x + lambda_y;

    /* 最大値ノルムに関する安定性を満たすλの上限 */
    lambda_limit = 1.0 / (2.0 * (1.0 - theta));

    if (lambda > lambda_limit)
        printf("注意: λ=%g>1/2(1-θ) となっています。\n", lambda);
    else
        printf("λ=%g\n", lambda);

    /* 初期値の設定 */
    for (i = 0; i <= N_x; i++) {
      x = a + i * h_x;
      for (j = 0; j <= N_y; j++)
        Uk[i][j] = f(x, c + j * h_y);
    }

    /* 連立1次方程式に現れる係数 
       γU_{ij}^{n+1}
         +α(U_{i+1,j}^{n+1}+U_{i-1,j}^{n+1})
         +β(U_{i,j+1}^{n+1}+U_{i,j-1}^{n+1})
      =γ'U_{ij}^n
         +α'(U_{i+1,j}^n+U_{i-1,j}^n)
         +β'(U_{i,j+1}^n+U_{i,j-1}^n)
     と書いたときのα,β,γ,α',β',γ' */
    gamma = 1.0 + 2.0 * theta * lambda;
    alpha = - theta * lambda_x;
    beta  = - theta * lambda_y;

    gamma_p = 1.0 - 2.0 * (1.0 - theta) * lambda;
    alpha_p = (1.0 - theta) * lambda_x;
    beta_p  = (1.0 - theta) * lambda_y;

    /* 係数行列の作成 */
    /* まず 0 クリア */
    for (p = 0; p < N; p++) {
        for (q = 0; q <= m; q++)
            A[p][q] = 0.0;
    }
    for (i = 0; i <= N_x; i++)
        for (j = 0; j <= N_y; j++) {
            L = phi(i, j);
	    /*
            if (j != 0)   A[L][L - m] = beta;
            if (i != 0)   A[L][L - 1] = alpha;
	    */
                          A[L][0] = gamma; /* A[L][L    ] = gamma; */
            if (i != N_x) A[L][1] = alpha; /* A[L][L + 1] = alpha; */
            if (j != N_y) A[L][m] = beta;  /* A[L][L + m] = beta;  */
        }

    /* 左または右 */
    for (j = 1; j < N_y; j++) {
        /* (0,j) */
        L = phi(0,j);
        A[L][0] /= 2.0; A[L][m] /= 2.0;
        /* A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0; */
        /* (N_x,j) */
        L = phi(N_x,j);
        A[L][0] /= 2.0; A[L][m] /= 2.0;
        /* A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0; */
    }

    /* 下または上 */
    for (i = 1; i < N_x; i++) {
        /* (i,0) */
        L = phi(i,0);
        A[L][0] /= 2.0; A[L][1] /= 2.0;
        /* A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0; */
        /* (i,N_y) */
        L = phi(i,N_y);
        A[L][0] /= 2.0; A[L][1] /= 2.0;
        /* A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0; */
    }
    /* 角の点 */
    /* 左下 */
    L = phi(0,0);
    A[L][0] /= 4.0; A[L][1] /= 2.0; A[L][m] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L+m] /= 2.0; */
    /* 左上 */
    L = phi(0,N_y);
    A[L][0] /= 4.0; A[L][1] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L-m] /= 2.0; */
    /* 右下 */
    L = phi(N_x,0);
    A[L][0] /= 4.0; A[L][m] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L+m] /= 2.0; */
    /* 右上 */
    L = phi(N_x,N_y);
    A[L][0] /= 4.0;
    /* A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L-m] /= 2.0; */

    /* 連立1次方程式の係数行列を表示する */
    if (N < 20) {
      printf("素朴に作った連立1次方程式の行列\n");
      for (p = 0; p < N; p++) {
        for (q = 0; q < N; q++)
          printf("%5.2f", A[p][q]);
        printf("\n");
      }
    }

    printf("備考: 1+2θλ=%5.2f, -θλx=%5.2f, -θλy=%5.2f\n",
           gamma, alpha, beta);

    printf("Tmax: "); scanf("%lf", &Tmax);
    printf("Δt: ");  scanf("%lf", &dt);
    skip = rint(dt / tau);
    if (skip == 0) skip = 1;

    nMax = rint(Tmax / tau);

    /* グラフィックス・ライブラリィ GLSC の呼び出し */
    g_init("Meta", 250.0, 160.0);
    g_device(G_BOTH);
    g_def_scale(0,
                0.0, 1.0, 0.0, 1.0,
                30.0, 70.0, 100.0, 72.0);
    g_def_scale(4,
                -1.0, 1.0, -1.0, 1.0,
                30.0, 30.0, 100.0, 100.0);
    g_def_line(0, G_BLACK, 0, G_LINE_SOLID);
    g_sel_scale(0);

    g_cls();
    g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
              150.0, 100.0, Uk, N_x + 1, N_y + 1,
              1, G_SIDE_NONE, 2, 1);

    /* 係数行列 LU 分解 */
    symbandlu(A, N, m + 1);
    if (cond + 1 == cond) {
        /* 条件数が大きければ、計算をあきらめる */
        printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
        return 0;
    }

    /* 時間に関するループ */
    for (n = 1; n <= nMax; n++) {
        t = n * tau; t_p = (n - 1) * tau;

        /* まず、素朴な連立1次方程式の右辺を用意する */
        /* 内部の格子点 */
        for (i = 1; i < N_x; i++)
            for (j = 1; j < N_y; j++) {
                L = phi(i,j);
                B[L] = gamma_p * Uk[i][j]
                    + alpha_p * (Uk[i + 1][j] + Uk[i - 1][j])
                    + beta_p * (Uk[i][j + 1] + Uk[i][j - 1]);
            }
        /* 以下、境界にある格子点での方程式を立てる。
         * 仮想格子点での値は境界条件を中心差分近似した方程式を
         * 用いて消去する
         *   ……右辺にも移項する量があるので、右辺について後で処理する */
        /* 下の辺、上の辺にある格子点 (角の点は含めない) */
        for (i = 1; i < N_x; i++) {
            x = a + i * h_x;
            /* (i, 0) */
            L = phi(i,0);
            B[L] = gamma_p * Uk[i][0]
                 + alpha_p * (Uk[i + 1][0] + Uk[i - 1][0])
                 + 2 * beta_p * Uk[i][1]
                 + 2 * h_y  * (beta * Phi(x, c, t) - beta_p * Phi(x, c, t_p));
            B[L] /= 2;
            /* (i, N_y) */
            L = phi(i,N_y);
            B[L] = gamma_p * Uk[i][N_y]
                 + alpha_p * (Uk[i + 1][N_y] + Uk[i - 1][N_y])
                 + 2 * beta_p * Uk[i][N_y - 1]
                 + 2 * h_y * (beta_p * Phi(x, d, t_p) - beta * Phi(x, d, t));
            B[L] /= 2;
        }
        /* 左の辺、右の辺にある格子点 (角の点は含めない) */
        for (j = 1; j < N_y; j++) {
            y = c + j * h_y;
            /* (0, j) */
            L = phi(0,j);
            B[L] = gamma_p * Uk[0][j]
                 + 2 * alpha_p * Uk[1][j]
                 + beta_p * (Uk[0][j + 1] + Uk[0][j - 1])
                 + 2 * h_x * (alpha * Phi(a, y, t) - alpha_p * Phi(a, y, t_p));
            B[L] /= 2;
            /* (N_x, j) */
            L = phi(N_x,j);
            B[L] = gamma_p * Uk[N_x][j]
                 + 2 * alpha_p * Uk[N_x - 1][j]
                 + beta_p * (Uk[N_x][j + 1] + Uk[N_x][j - 1])
                 + 2 * h_x * (alpha_p * Phi(b, y, t_p) - alpha * Phi(b, y, t));
            B[L] /= 2;
        }
        /* 左下 */
        L = phi(0,0);
        B[L] = gamma_p * Uk[0][0]
             + 2 * alpha_p * Uk[1][0] + 2 * beta_p * Uk[0][1]
             + 2 * h_x * (alpha * Phi(a,c,t) - alpha_p * Phi(a,c,t_p))
             + 2 * h_y * (beta * Phi(a,c,t) - beta_p * Phi(a,c,t_p));
            B[L] /= 4;
        /* 左上 */
        L = phi(0,N_y);
        B[L] = gamma_p * Uk[0][N_y]
             + 2 * alpha_p * Uk[1][N_y] + 2 * beta_p * Uk[0][N_y - 1]
             + 2 * h_x * (alpha * Phi(a,d,t) - alpha_p * Phi(a,d,t_p))
             + 2 * h_y * (beta_p * Phi(a,d,t_p) - beta * Phi(a,d,t));
            B[L] /= 4;
        /* 右下 */
        L = phi(N_x,0);
        B[L] = gamma_p * Uk[N_x][0]
             + 2 * alpha_p * Uk[N_x - 1][0] + 2 * beta_p * Uk[N_x][1]
             + 2 * h_x * (alpha_p * Phi(b,c,t_p) - alpha * Phi(b,c,t))
             + 2 * h_y * (beta * Phi(b,c,t) - beta_p * Phi(b,c,t_p));
	B[L] /= 4;
        /* 右上 */
        L = phi(N_x,N_y);
        B[L] = gamma_p * Uk[N_x][N_y]
             + 2 * alpha_p * Uk[N_x - 1][N_y] + 2 * beta_p * Uk[N_x][N_y - 1]
             + 2 * h_x * (alpha_p * Phi(b,d,t_p) - alpha * Phi(b,d,t))
             + 2 * h_y * (beta_p * Phi(b,d,t_p) - beta * Phi(b,d,t));
	B[L] /= 4;

        /* A vector_U = B を解く */
        symbandsolve(A, B, N, m + 1);
        /* */
        for (i = 0; i <= N_x; i++)
            for (j = 0; j <= N_y; j++)
                Uk[i][j] = B[phi(i,j)];

        /* データを数値で表示 */
        if (n % skip == 0) {
#ifdef  PRINT
            for (i = 0; i <= N_x; i++) {
                for (j = 0; j <= N_y; j++)
                    printf(" %5.1f", Uk[i][j]);
                printf("\n");
            }
#endif

            /* 鳥瞰図を描く */
            g_cls();
            g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
                      150.0, 100.0, Uk, N_x + 1, N_y + 1,
                      1, G_SIDE_NONE, 2, 1);
        }
    }
    /* マウスでクリックされるのを待つ */
    g_sleep(-1.0);
    /* ウィンドウを消す */
    g_term();

    return 0;
}

/* 初期値 */
double f(double x, double y)
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

/* Neumann 境界値 */
double Phi(double x, double y, double t)
{
    /* 同次 Neumann 境界条件 */
    return 0.0;
}
