/*
 * draw-Bessel-glsc.c --- Bessel 関数のグラフを描く。
 *   コンパイル: ccmg draw-Bessel-glsc.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G_DOUBLE
#include <glsc.h>

int main(int argc, char **argv)
{
  int n, n0, n1, N, i;
  double a, b, dx, margin, height, x;

  a = 0.0;
  b = 20.0;
  height = 1.2;
  margin = (b - a) / 10;

  g_init("BESSEL", 320.0, 220.0);
  g_device(G_BOTH);
  g_def_scale(0, a - margin, b + margin, - height, height,
	      10.0, 10.0, 300.0, 200.0);
  g_sel_scale(0);
  g_def_line(0, G_BLACK, 0, G_LINE_SOLID);

  if (argc == 3) {
    n0 = atoi(argv[1]);
    n1 = atoi(argv[2]);
  }
  else {
    printf("何次から何次までの Bessel 関数のグラフを書きますか?: ");
    scanf("%d %d", &n0, &n1);
  }

  N = 200;
  dx = (b - a) / N;
  for (n = n0; n <= n1; n++) {
    g_move(a, jn(n, a));
    for (i = 1; i <= N; i++) {
      x = a + i * dx;
      g_plot(x, jn(n, x));
    }
  }
  g_sleep(-1.0);
  g_term();
  
  return 0;
}
