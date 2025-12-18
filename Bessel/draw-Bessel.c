/*
 * draw-Bessel.c --- Bessel 関数のグラフを描く。
 *   コンパイル: ccx draw-Bessel.c
 */

#include <stdio.h>
#include <math.h>

main()
{
  int n, n0, n1, N, i;
  double a, b, dx, margin, height, x, y;

  a = 0.0;
  b = 20.0;
  height = 2.0;
  N = 200;

  dx = (b - a) / N;
  margin = (b - a) / 10;

  openpl();
  fspace2(a - margin, - height, b + margin, height);

  printf("何次から何次までの Bessel 関数のグラフを書きますか?: ");
  scanf("%d %d", &n0, &n1);

  for (n = n0; n <= n1; n++) {
    fmove(0.0, jn(n, 0.0));
    for (i = 1; i <= N; i++) {
      x = a + i * dx;
      y = jn(n, x);
      fcont(x, y);
    }
  }
  mkplot("Bessel.plot");
  closepl();
}
