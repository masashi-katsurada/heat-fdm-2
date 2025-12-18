/*
 * call_gnuplot.c --- gnuplot を用いて円板で定義された関数を可視化する
 *
 *  written by mk.
 *
 *  GNUPLOT manual (in Japanese)
 * http://www.meiji.ac.jp/servers/math/staffs/mk/on-computer/gnuplotj/index.html
 *
 */

#include <sys/stat.h>		/* S_IRUSR etc */
#include <stdio.h>
#include <unistd.h>		/* unlink() */
#include <math.h>		/* atan() etc */
#include <signal.h>
#include <matrix.h>
#include "call_gnuplot2.h"

int msleep(int);

#define PANIC(a) do { \
                       perror(a); \
                       if (temp_name) unlink(temp_name);\
                       exit(1);\
               } while(0)

void goodbye(int dummy)
{
    close_gnuplot();
    exit(0);
}

double pi;
char *temp_name = NULL;
FILE *gnuplot, *data;

int open_gnuplot()
{
    pi = 4.0 * atan(1.0);
    if ((temp_name = tmpnam((char *) 0)) == 0)
	PANIC("tmpnam failed");
    /*
     * stat.h
     *  S_IRUSR R for owner
     *  S_IWUSR W for owner
     */
    if (mkfifo(temp_name, S_IRUSR | S_IWUSR) != 0)
	PANIC("mkfifo failed");
    gnuplot = popen("gnuplot", "w");
    fprintf(gnuplot, "set data style lines\n");
    fprintf(gnuplot, "set parametric\n");
    fprintf(gnuplot, "set hidden3d\n");
    fprintf(gnuplot, "set contour\n");
    fprintf(gnuplot, "set cntrparam levels 10\n");
    fflush(gnuplot);
    signal(SIGINT, goodbye);
    return 0;
}

int disk(int nr, int nt, double R, matrix u, char *label)
{
    int i, j;
    double r, theta, x, y, z;
    double dr, dt;
    FILE *data;

    dr = R / nr;
    dt = 2.0 * pi / nt;

    fprintf(gnuplot, "set nolabel\n");
    fprintf(gnuplot, "set label \"%s\" at screen 0.4,0.85\n", label);
    fprintf(gnuplot, "set zrange [-1:3]\n");
    fprintf(gnuplot, "clear\n");
    fprintf(gnuplot, "splot \"%s\" with lines\n", temp_name);
    fflush(gnuplot);

    data = fopen(temp_name, "w");

    for (i = 0; i <= nr; i++) {
	r = i * dr;
	for (j = 0; j <= nt; j++) {
	    theta = (j) * dt;
	    x = r * cos(theta);
	    y = r * sin(theta);
	    z = u[i][j];
	    fprintf(data, "%f %f %f\n", x, y, z);
	}
	fprintf(data, "\n");
    }
    fclose(data);
    msleep(100);
    return 0;
}

void close_gnuplot()
{
    pclose(gnuplot);
    unlink(temp_name);
}

/*
 * msleep.c -- ミリ秒単位のタイマー
 */

#include <sys/types.h>
#include <sys/time.h>

int msleep(int ms)
{
    struct timeval timeout;
    timeout.tv_sec = ms / 1000;
    timeout.tv_usec = (ms % 1000) * 1000;
    if (select(0, (fd_set *) 0, (fd_set *) 0, (fd_set *) 0, &timeout) < 0) {
	perror("usleep");
	return -1;
    }
    return 0;
}
