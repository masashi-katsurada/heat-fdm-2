#include <stdio.h>
#include <math.h>
#include <matrix.h>
#include "call_gnuplot.h"

int main()
{
    double r, dr, dt, theta, t, pi, f(double, double);
    matrix u;
    char label[80];

    int i, j, k, nr, nt;
    pi = 4.0 * atan(1.0);

    nr = 20;
    nt = 80;
    u = new_matrix(nr + 1, nt + 1);

    dr = 1.0 / nr;
    dt = 2.0 * pi / nt;

    open_gnuplot();
    for (k = 0; k < 100; k++) {
	fprintf(stderr, "k=%d\n", k);
	t = k * 0.05;

	for (i = 0; i <= nr; i++) {
	    r = i * dr;
	    for (j = 0; j <= nt; j++) {
		theta = j * dt + k * 0.1;
		u[i][j] = exp(-t) * f(r * cos(theta), r * sin(theta));
	    }
	}
	sprintf(label, "t=%g", t);
	disk(nr, nt, u, label);
    }

    getchar();

    close_gnuplot();
    return 0;
}

double f(double x, double y)
{
    return x * x - y * y;
}
