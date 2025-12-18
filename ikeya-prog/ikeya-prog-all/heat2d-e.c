/*
 * heat2d-e.c -- solve the heat equation in a two dimensional square region
 * by explicit finite difference method. (version 2)
 *
 *	To compile this program on WS's in 6701,
 *		gcc -o heat2d-e heat2d-e.c -lmatrix -lglscd -lX11 -lm
 *      or try
 *              ccmg heat2d-e.c
 */

#include <stdio.h>
#include <math.h>

/* to use matrix, new_matrix() */
#include <matrix.h>

/* to use GLSC */
#define G_DOUBLE
#include <glsc.h>

int main()
{
    int             Nx, Ny, i, j, n, skip, nMax;
    double          hx, hy, lambdax, lambday, lambda, tau, Tmax, t, c, dt;
    matrix          u, newu;
    double          f(double, double);

    printf("Nx, Ny: "); scanf("%d %d", &Nx, &Ny);
    if ((u = new_matrix(Nx + 1, Ny + 1)) == NULL) {
        fprintf(stderr, "配列 u を確保できませんでした。");
        exit(1);
    }
    if ((newu = new_matrix(Nx + 1, Ny + 1)) == NULL) {
        fprintf(stderr, "配列 newu を確保できませんでした。");
        exit(1);
    }
    printf("Tmax: "); scanf("%lf", &Tmax);

    hx = 1.0 / Nx;
    hy = 1.0 / Ny;

    printf("τ(≦%g): ",
	 0.5 / (1 / (hx * hx) + 1 / (hy * hy)));
    scanf("%lf", &tau);

    lambdax = tau / (hx * hx);
    lambday = tau / (hy * hy);
    lambda = lambdax + lambday;
    c = 1.0 - 2.0 * lambda;
    printf("λ=%g になりました。\n", lambda);

    printf("Δt: ");
    scanf("%lf",  &dt);
    skip = rint(dt / tau);
    if (skip == 0) {
      printf("Δtが小さすぎるので、Δt=τとします。\n");
      skip = 1;
    }
    dt = skip * tau;

    g_init("Meta", 250.0, 160.0);
    g_device(G_BOTH);
    g_def_scale(0, 0.0, 1.0, 0.0, 1.0, 30.0, 70.0, 100.0, 72.0);
    g_def_scale(4, -1.0, 1.0, -1.0, 1.0, 30.0, 30.0, 100.0, 100.0);
    g_def_line(0, G_BLACK, 0, G_LINE_SOLID);
    g_sel_scale(0);

    for (i = 0; i <= Nx; i++)
        for (j = 0; j <= Ny; j++)
            u[i][j] = f(i * hx, j * hy);

    g_cls();
    g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
        150.0, 100.0, u, Nx + 1, Ny + 1,
        1, G_SIDE_NONE, 2, 1);

    nMax = rint(Tmax / tau);
    for (n = 1; n <= nMax; n++) {
        for (i = 1; i < Nx; i++)
            for (j = 1; j < Ny; j++)
                newu[i][j] = c * u[i][j]
                    + lambdax * (u[i + 1][j] + u[i - 1][j])
		    + lambday * (u[i][j + 1] + u[i][j - 1]);
        for (i = 0; i <= Nx; i++)
            newu[i][0] = newu[i][Ny] = 0.0;
        for (j = 0; j <= Ny; j++)
            newu[0][j] = newu[Nx][j] = 0.0;
        for (i = 0; i <= Nx; i++)
            for (j = 0; j <= Ny; j++)
                u[i][j] = newu[i][j];
        if (n % skip == 0) {
            g_cls();
            g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
                150.0, 100.0, u, Nx + 1, Ny + 1,
                1, G_SIDE_NONE, 2, 1);
        }
        t = tau * n;
    }

    /* マウスでクリックされるのを待つ */
    g_sleep(-1.0);
    /* ウィンドウを消す */
    g_term();
    return 0;
}

double f(double x, double y)
{
    return 16 * x * (1 - x) * y * (1 - y);
}
