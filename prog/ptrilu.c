/*
 * ptrilu.c --- 周期三重対角行列の LU 分解 (version 1.1)
 */

#include <matrix.h>

/*
   ptrilu(int n, vector al, vector ad, vector au, vector ab, vector ar)

   機能

   A =

   ad[0] au[0]                               ar[0]
   al[1] ad[1] au[1]                         ar[1]
   al[2] ad[2] au[2]                         ar[2]

                    al[n-3] ad[n-3] au[n-3] ar[n-3]
                            al[n-2] ad[n-2] au[n-2]
   ab[0] ab[1] ab[2] ...... ab[n-3] al[n-1] ad[n-1]

   という周期三重対角行列 A を LU 分解する。

   入力
   int n;                          行列 A の次数
   vector al, ad, au, ar, ab;      行列 A の成分 (上を参照)
   al, ad, au は n 次元ベクトル
   ar, ab は n-2 次元ベクトル
   出力
   vector al, ad, au, ar, ab;      行列 L, U の成分 (上を参照)

 */

void ptrilu0(int n, vector al, vector ad, vector au, vector ab, vector ar)
{
    int k, kp1, nm1 = n - 1, nm2 = n - 2;

    for (k = 0; k < nm2; k++) {
        kp1 = k + 1;

        al[kp1] /= ad[k];
        ad[kp1] -= al[kp1] * au[k];
        if (kp1 < nm2)
            ar[kp1] -= al[kp1] * ar[k];
        else
            au[kp1] -= al[kp1] * ar[k];

        ab[k] /= ad[k];
        if (kp1 < nm2)
            ab[kp1] -= ab[k] * au[k];
        else
            al[nm1] -= ab[k] * au[k];
        ad[nm1] -= ab[k] * ar[k];
    }
    al[nm1] /= ad[nm2];
    ad[nm1] -= al[nm1] * au[nm2];
    /* */
#ifdef NEVER
    ar[nm2] = au[nm2];
    ab[nm2] = al[nm1];
    ab[nm1] = ar[nm1] = ad[nm1];
#endif
}

void ptrisol0(int n,
	      vector al, vector ad, vector au, vector ab, vector ar,
	      vector b)
{
    int k, nm1 = n - 1;

    for (k = 0; k < nm1; k++) {
        b[k + 1] -= al[k + 1] * b[k];
        if (k + 1 != nm1)
            b[nm1] -= ab[k] * b[k];
    }
    b[nm1] /= ad[nm1];
    for (k = n - 2; k >= 0; k--) {
        b[k] -= au[k] * b[k + 1];
        if (k + 1 != nm1)
            b[k] -= ar[k] * b[nm1];
        b[k] /= ad[k];
    }
}

void ptrilu1(int n, vector al, vector ad, vector au, vector ab, vector ar)
{
    int k, kp1, nm1 = n - 1;

    for (k = 1; k < nm1; k++) {
        kp1 = k + 1;

        al[kp1] /= ad[k];
        ad[kp1] -= al[kp1] * au[k];
        if (kp1 < nm1)
            ar[kp1] -= al[kp1] * ar[k];
        else
            au[kp1] -= al[kp1] * ar[k];

        ab[k] /= ad[k];
        if (kp1 < nm1)
            ab[kp1] -= ab[k] * au[k];
        else
            al[n] -= ab[k] * au[k];
        ad[n] -= ab[k] * ar[k];
    }
    al[n] /= ad[nm1];
    ad[n] -= al[n] * au[nm1];
}

void ptrisol1(int n,
	      vector al, vector ad, vector au, vector ab, vector ar,
	      vector b)
{
    int k;

    for (k = 1; k < n; k++) {
        b[k + 1] -= al[k + 1] * b[k];
        if (k + 1 != n)
            b[n] -= ab[k] * b[k];
    }
    b[n] /= ad[n];
    for (k = n - 1; k >= 1; k--) {
        b[k] -= au[k] * b[k + 1];
        if (k + 1 != n)
            b[k] -= ar[k] * b[n];
        b[k] /= ad[k];
    }
}
