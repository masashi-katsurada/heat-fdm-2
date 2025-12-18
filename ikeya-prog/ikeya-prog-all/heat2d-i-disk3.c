/*
 * heat2d-i-en3.c ---円盤状の2次元熱方程式を陰解法で解く
 *         j方向番号付け 効率を考えたプログラム       
 *   コンパイルするには
 *     gcc -c lu.c
 *     gcc -o heat2d-i-en heat2d-i-en.c lu.o -lmatrix -lglscd -lX11 -lm
 *   あるいは安直には
 *     ccmg heat2d-i-en3.c lu.c   
 */
#include <stdio.h>
#include <math.h>
#include <matrix.h>
#include "bandlu.h"

#define psi(i,j) ((i-1)*mm + j+1)

double u0(double, double);
double exactu(double, double, double);
double maxnorm(int, int, matrix);
double pi;

int main()
{
  double ri ,ri2;
  int N_r, N_p, mm, mm2, NN, i, j, p, q, n, L;
  matrix Uk, A;
  double *B, *vector_U;
  double h_r, h_p, lambda_r, lambda_p, lambda, lambda_limit, tau;
  double theta;

  
  /* πの値*/
  pi = 4.0 * atan(1.0);

  /* 区間の分割数 */
   printf("Nr, Ntheta: "); scanf("%d %d", &N_r, &N_p);

  mm = N_p;
  NN = (N_r-1)*N_p + 1;

  /* 空間の刻み幅 */
  h_r = 1.0 / N_r;
  h_p = 2.0 * pi / N_p;


  /* 行列、ベクトルを記憶する変数のメモリー割り当て */
    if ((Uk = new_matrix(N_r+ 1, N_p+1)) == NULL) {
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

    /* θ法の重みの決定 */
    printf("θ (0≦θ≦1):  "); scanf("%lf", &theta);

    if (theta == 1.0) {
        printf("τ: "); scanf("%lf", &tau);
    } else {
      printf("τ(≦%g≡最大値ノルムに関する安定性条件を満たすτの上限): ", 
             0.5 * (h_r * h_r * h_p * h_p) / (1 + h_p * h_p));
        
       scanf("%lf", &tau);
    }
	
    /* λr, λθ */
    lambda_r = tau / (h_r * h_r);
    lambda_p = tau / (h_p * h_p);
    lambda = lambda_r + lambda_p;
    
/* 初期値の設定 */
    for (i = 0; i <= N_r; i++) {
      ri = i * h_r;
      for (j = 0; j <= N_p-1; j++)
        Uk[i][j] = u0(ri, j * h_p);
    }

    /* 係数行列の値のセット */
    /* 0 クリア */
    for(p = 0; p < NN; p++){
      for(q = 0; q < NN; q++)
	A[p][q] = 0.0;
    }
    
   /* 原点の係数 */ 
    A[0][0] = 1.0 + 4.0 * theta * lambda_r;
    
    for(j = 0; j <= N_p-1; j++){
      L =  psi(1,j);
    A[0][L] = - (4.0 * theta * lambda_r) / N_p;
      }

 /* 内点での係数 */
    for(i = 1; i < N_r; i++){
      ri = (i)*h_r;
      ri2 = ri*ri; 
      for(j = 0; j <= N_p-1; j++){
	L = psi(i,j);
	int L1 = L-1;
	int L2 = L+1;
	if(j == 0)
	  L1 = psi(i,N_p-1);

	if(j == N_p-1)
	  L2 = psi(i,0);

	if(i != N_r-1)
	   A[L][L+mm] = - theta* lambda_r * (1.0 + h_r/(2.0*ri));

	if(i == 1){
	  A[L][0] = - theta* lambda_r * (1.0 - h_r/(2.0*ri));
	}

	if(i != 1)
	  A[L][L-mm] = - theta* lambda_r * (1.0 - h_r/(2.0*ri));


	A[L][L] = 1.0 + 2.0*theta*lambda_r + (2.0*theta*lambda_p) / ri2;
	A[L][L1] = - theta * lambda_p / ri2;
	A[L][L2] = - theta * lambda_p / ri2; 

      }
    }

    /* 連立方程式の係数行列の表示 */
    if (NN <= 30){
      printf(" 連立１次方程式の行列\n");
      for(p = 0; p < NN; p++){
	for(q = 0; q < NN; q++)
          printf(" %4.1f", A[p][q]);
	printf("\n");
      }
    }

 
   /* 係数行列 LU 分解 */
	bandlu(A, NN, mm);
  
	/* 連立方程式の係数行列の表示 */
    if (NN <= 30){
      printf(" 連立１次方程式の行列\n");
      for(p = 0; p < NN; p++){
	for(q = 0; q < NN; q++)
          printf(" %4.1f", A[p][q]);
	printf("\n");
      }
    }
     /* 連立1次方程式の右辺を用意する */
 /* 内部の格子点 */ 
     for (i = 1; i < N_r; i++){
     ri = i * h_r; 
     ri2 = ri * ri;
     for(j = 0; j<= N_p-1; j++){
	 L = psi(i,j);
	 int jm1 = j-1;
	 int jm2 = j+1;
     if(jm1 == -1)
	   jm1 = N_p - 1;
     if(jm2 == N_p)
           jm2 = 0;
	 B[L] =  (1.0 - 2.0*(1.0 - theta)*lambda_r - (2.0*(1.0-theta)*lambda_p) / ri2)*Uk[i][j] + (1.0-theta)*lambda_r*((1.0 + h_r/2.0*ri)*Uk[i+1][j] + (1.0 - h_r/2.0*ri)*Uk[i-1][j]) +  ((1.0-theta)*lambda_p / ri2)*(Uk[i][jm2] + Uk[i][jm1]);
     }
     }

	 /* 原点 */
     B[0] = (1.0 - 4.0*(1.0 - theta)*lambda_r)*Uk[0][0];
   for(j = 0; j <= N_p-1; j++){
     B[0] += (4.0*(1.0 -theta)*lambda_r/N_p)*Uk[1][j];
   }


/* B[L]の表示 */
   printf(" B[L]の表示\n");
    for(i = 1; i <= N_r-1; i++)
      for(j = 0; j <= N_p-1; j++)	
      printf("B(%d, %d) = %g\n", i, j, B[psi(i,j)]);

/* A vector_U = B を解く */
   bandsolve(A, B, NN, mm);
 

  /* 値の更新 */
   for(i = 1; i< N_r; i++)
      for(j = 0; j <= N_p-1; j++)
       Uk[i][j] = B[psi(i,j)];
      
   for(j = 0; j <= N_p-1; j++)
     Uk[N_r][j] = 0.0;

   for(j = 1; j <=N_p - 1; j++)
     Uk[0][j] = Uk[0][0];

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
