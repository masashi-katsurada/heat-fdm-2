/* trid-lu.c -- 3項方程式を Gauss の消去法で解く */

#include "trid-lu.h"

/* 3項方程式 (係数行列が三重対角である連立1次方程式のこと) Ax=b を解く
 *
 *   入力
 *     n: 未知数の個数
 *     al,ad,au: 連立1次方程式の係数行列
 *       (al: 対角線の下側 i.e. 下三角部分 (lower part)
 *        ad: 対角線       i.e. 対角部分 (diagonal part)
 *        au: 対角線の上側 i.e. 上三角部分 (upper part)
 *        つまり
 *
 *         ad[0] au[0]   0   .................. 0
 *         al[1] ad[1] au[1]   0  ............. 0
 *           0   al[2] ad[2] au[2]  0 ......... 0
 *                         ....................
 *                         al[n-2] ad[n-2] au[n-2]
 *                           0     al[n-1] ad[n-1]

 *        al[i] = A_{i,i-1}, ad[i] = A_{i,i}, au[i] = A_{i,i+1},
 *        al[0], au[n-1] は意味がない)
 *
 *     b: 連立1次方程式の右辺の既知ベクトル
 *        (添字は 0 から。i.e. b[0],b[1],...,b[n-1] にデータが入っている。)
 *   出力
 *     al,ad,au: 入力した係数行列を LU 分解したもの
 *     b: 連立1次方程式の解
 *   能書き
 *     一度 call すると係数行列を LU 分解したものが返されるので、
 *     以後は同じ係数行列に関する連立1次方程式を解くために、
 *     関数 trisol() が使える。
 *   注意
 *     Gauss の消去法を用いているが、ピボットの選択等はしていな
 *     いので、ピボットの選択をしていないので、係数行列が正定値である
 *     などの適切な条件がない場合は結果が保証できない。
 */

void trid(int n, double *al, double *ad, double *au, double *b)
{
  trilu(n,al,ad,au);
  trisol(n,al,ad,au,b);
}

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

/* LU 分解済みの三重対角行列を係数に持つ3項方程式を解く */
void trisol(int n, double *al, double *ad, double *au, double *b)
{
  int i, nm1 = n - 1;
  /* 前進消去 (forward elimination) */
  for (i = 0; i < nm1; i++) b[i + 1] -= b[i] * al[i + 1];
  /* 後退代入 (backward substitution) */
  b[nm1] /= ad[nm1];
  for (i = n - 2; i >= 0; i--) b[i] = (b[i] - au[i] * b[i + 1]) / ad[i];
}
