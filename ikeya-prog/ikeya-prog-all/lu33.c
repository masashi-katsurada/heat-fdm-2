/*
 * heat2d-i-en.c ---円盤状の2次元熱方程式を陰解法で解く
 *                *
 *   コンパイルするには
 *     gcc -c lu.c
 *     gcc -o heat2d-i-en heat2d-i-en.c lu.o -lmatrix -lglscd -lX11 -lm
 *   あるいは安直には
 *     ccmg heat2d-i-en.c lu.c   
 */
#include <stdio.h>
#include <math.h>
#include <matrix.h>
#include "lu.h"

#define psi(i,j) (j*mm+i)

double u0(double, double);
double exactu(double, double, double);
double maxnorm(int, int, matrix);
double pi;

int main()
{
  double ri ,ri2;
  int N_r, N_p, mm, mm2, NN, i, j, p, q, n, skip, nMax, L;
  matrix Uk, A;
  double *B, *vector_U, cond;
  int *iwork;
  double h_r, h_p, lambda_r, lambda_p, lambda, lambda_limit, tau, t, Tmax, dt;
  double theta;

  
   NN = 3;
  /* 行列、ベクトルを記憶する変数のメモリー割り当て */
    if ((Uk = new_matrix(NN, NN)) == NULL) {
        fprintf(stderr, "数列 U^k を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((A = new_matrix(NN, NN)) == NULL) {
        fprintf(stderr, "係数行列 A を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((B = (double *) malloc(sizeof(double) * NN)) == NULL) {
        fprintf(stderr, "B を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((vector_U = (double *) malloc(sizeof(double) * NN)) == NULL) {
        fprintf(stderr, "vector_U を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((iwork = (int *) malloc(sizeof(int) * NN)) == NULL) {
        fprintf(stderr, "iwork を記憶する領域の確保に失敗\n");
        exit(1);
    }
 

    /* 係数行列の値のセット */
    /* 0 クリア */
    for(p = 0; p < NN; p++){
      for(q = 0; q < NN; q++)
	A[p][q] = 0.0;
             /* Dirichle境界条件 原点のコピーであらかじめ単位行列を作る */
    }
    A[0][0] = 2; A[0][1] = 3; A[0][2] = -1;
    A[1][0] = 4; A[1][1] = 4; A[1][2] = -3;
    A[2][0] = -2; A[2][1] = 3; A[2][2] = -1;
    /* 連立方程式の係数行列の表示 */
    if (NN <= 20){
      printf(" 連立１次方程式の行列\n");
      for(p = 0; p < NN; p++){
	for(q = 0; q < NN; q++)
          printf(" %4.1f", A[p][q]);
	printf("\n");
      }
    }

 
   /* 係数行列 LU 分解 */
	decomp(NN, A,  &cond, iwork, B);
	if(cond + 1 == cond){
	  printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
	  return 0;
	}
  
	/* 連立方程式の係数行列の表示 */
    if (NN <= 20){
      printf(" 連立１次方程式の行列\n");
      for(p = 0; p < NN; p++){
	for(q = 0; q < NN; q++)
          printf(" %4.1f", A[p][q]);
	printf("\n");
      }
    }
     /* 連立1次方程式の右辺を用意する */
    B[0] = 5; B[1] = 3; B[2] = 1;

/* B[L]の表示 */
   printf(" B[L]の表示\n");
   for(i = 0; i < NN ; i++){
	L = i;	
      printf("B = %g\n",  B[L]);
   }
/* A vector_U = B を解く */
   solve(NN, A, B, iwork);

/* B[L]の表示 */
   printf(" B[L]の表示\n");
      for(i = 0; i < NN; i++)	
      printf("B = %g\n", B[i]);


}
 /* 厳密解を計算する関数 */

#define mu01    (2.404825557695771998)
#define mu11    (3.831705970207510692)

     /* 初期データ */
double u0(double r, double phi)
  {
    return j0(mu01 * r) + j1(mu11 * r) * (cos(phi) + sin(phi));
  }

    /* 厳密解 */
double exactu(double r, double phi, double t)
	  {
	    return exp(- mu01 * mu01 * t) * j0(mu01 * r) + exp(- mu11 * mu11 * t)*j1(mu11 * r) * (cos(phi) + sin(phi));
	  }

double maxnorm(int m, int n, matrix Uk)
      {
	int i, j, i0, j0;
	double tmpmax, absu;
	i0 = 0;
	j0 = 0;
	tmpmax = fabs(Uk[0][0]);
	for(i = 0; i < m; i++)
	  for(j = 0; j < n; j++)
	    if((absu = fabs(Uk[i][j])) > tmpmax) {
	      tmpmax = absu;
	      i0 = i;
	      j0 = j;
	    }
	printf("(i,j)=(%d,%d)", i0, j0);
	return tmpmax;
      }

