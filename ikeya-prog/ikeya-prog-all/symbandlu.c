/*
 * symbandlu.c --- 菊地文雄『有限要素法概説』サイエンス社
 */

#include <matrix.h>

/* the Gauss elimination method for symmetric band matrices */
/* N は未知数の個数、 m は半バンド幅+1 */

void symbandlu(matrix at, int N, int m)
{
    int N1, i, j, k, mm;
    double aa;
    /* forward elimination; */
    N1 = N - 1;
    for (i = 0; i < N1; i++) {
        /* mm = min(i+m,N) */
        mm = i + m; if (mm > N) mm = N;
        for (j = i + 1; j < mm; j++) {
            aa = at[i][j-i] / at[i][0]; /* aa = a_{i,j} / a_{i,i} */
            for (k = j; k < mm; k++)
                at[j][k-j] -= aa * at[i][k-i]; /* a_{j,k} -= aa * a_{i,k} */
        }
    }
}

void symbandsolve(matrix at, vector f, int N, int m)
{
    int N1, i, j, k, mm;
    double aa;
    /* forward elimination; */
    N1 = N - 1;
    for (i = 0; i < N1; i++) {
        mm = i + m; if (mm > N) mm = N;
        for (j = i + 1; j < mm; j++) {
            aa = at[i][j-i] / at[i][0]; /* aa = a_{i,j} / a_{i,i} */
            f[j] -= aa * f[i];
        }
    }
    /* backward substitution */
    f[N-1] /= at[N-1][0]; /* f_{N-1} /= a_{N-1,N-1} */
    for (i = N - 2; i >= 0; i--) {
        mm = i + m;
        if (mm > N) mm = N;
         for (j = i + 1; j < mm; j++)
            f[i] -= at[i][j-i] * f[j]; /* f_i -= at_{i,j} * f_j */
         f[i] /= at[i][0]; /* f_i /= a_{i,i} */
    }
}
