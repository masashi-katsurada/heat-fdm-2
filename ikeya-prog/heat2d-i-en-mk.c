/*
 * heat2d-i-en-mk.c ---円盤状の2次元熱方程式を陰解法で解く
 *   heat2d-i-en.c の修正版
 *
 *    i方向番号付け のプログラム            
 *   コンパイルするには
 *     gcc -c lu.c
 *     gcc -c call_gnuplot.c
 *  
 *     ccmg heat2d-i-en-mk.c lu.c call_gnuplot.c  
 */

#include <stdio.h>
#include <math.h>
#include <matrix.h>
#include "lu.h"
#include "call_gnuplot.h"

#define  psi(i,j) ((j) * mm + (i))

double u0(double, double);
double exactu(double, double, double);
double maxnorm(int, int, matrix);
double pi;

int main()
{
  double ri ,ri2, phi_j;
  int N_r, N_p, mm, NN, i, j, p, q, n, skip, nMax, L;
  matrix Uk, A;
  double *B, *vector_U, cond;
  int *iwork;
  double h_r, h_p, lambda_r, lambda_p, lambda, tau, t, Tmax, dt, M, ex;
  double theta;
  char label[200];

  /* πの値*/
  pi = 4.0 * atan(1.0);

  /* 区間の分割数 */
  printf("Nr, Ntheta: "); scanf("%d %d", &N_r, &N_p);

  mm = N_r+1;
  NN = (N_r+1)*(N_p);

  /* 空間の刻み幅 */
  h_r = 1.0 / N_r;
  h_p = 2.0 * pi / N_p;

  /* 行列、ベクトルを記憶する変数のメモリー割り当て */
  if ((Uk = new_matrix(N_r+1, N_p+1)) == NULL) {
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

  /* θ法の重みの決定 */
  printf("θ (0≦θ≦1):  "); scanf("%lf", &theta);

  if (theta == 1.0) {
    printf("τ: "); scanf("%lf", &tau);
  }
  else {
    printf("τ(≦%g        ≡安定性条件): ", 
           (h_r * h_r * h_p * h_p) / (2.0 *(1 - theta) * (1 + h_p * h_p)));
    scanf("%lf", &tau);
  }
        
  /* λr, λθ */
  lambda_r = tau / (h_r * h_r);
  lambda_p = tau / (h_p * h_p);
  /* 結果を出力する時間間隔の決定 */
  printf("Tmax: "); scanf("%lf", &Tmax);
  printf("Δt(>=%g): ", tau);  scanf("%lf", &dt);
  if (dt < tau) {
    dt = tau;
  }
  skip = rint(dt / tau); 

  /* GNUPLOTの準備 */
  open_gnuplot(); 

  /* 初期値の設定 */
  for (i = 0; i <= N_r; i++) {
    ri = (i) * h_r;
    for (j = 0; j <= N_p; j++)
      Uk[i][j] = u0(ri, (j) * h_p);
  }
  /* 初期のグラフ */
  disk(N_r, N_p, Uk,"t=0");
  
  /* 係数行列の値のセット */
  /* 0 クリア */
  for (p = 0; p < NN; p++) {
    for (q = 0; q < NN; q++)
      A[p][q] = 0.0;
    A[p][p] = 1.0;
    /* Dirichle境界条件 原点のコピーのためあらかじめ単位行列を作る */
  }
  /* 原点の係数 */ 
  A[0][0] = 1.0 + 4.0*theta*lambda_r;
    
  for (j = 0; j <= N_p - 1; j++) {
    L =  psi(1,j);
    A[0][L] = - (4.0*theta*lambda_r) / N_p;
  }
 
  /* 原点のコピーの係数 */
  for (j = 1; j <= N_p - 1; j++) {
    L = psi(0,j);
    A[L][0] = - 1.0;
  }

  /* 内点での係数 */
  for (i = 1; i< N_r; i++) {
    ri = (i)*h_r;
    ri2 = ri*ri; 
    for (j = 0; j <= N_p - 1; j++) {
      int L1, L2;
      L = psi(i,j);
      L1 = L-mm;
      L2 = L+mm;
      if (j == 0)
        L1 = psi(i,N_p-1);
      /* j=0の時j-1が-1となるため-1の代わりにNθ-1を代入する */
      if(j == N_p-1)
        L2 = psi(i,0);
      /* j=N_p-1の時j+1がN_pとなるためN_pの代わりに0を代入する */
        
      A[L][L] = 1.0 + 2.0*theta*lambda_r + (2.0*theta*lambda_p) / ri2;
      A[L][L+1] = - theta* lambda_r * (1.0 + h_r/(2.0*ri));
      A[L][L-1] = - theta* lambda_r * (1.0 - h_r/(2.0*ri));
      A[L][L1] = - theta * lambda_p / ri2;
      A[L][L2] = - theta * lambda_p / ri2;
    }
  }

  /* 係数行列 LU 分解 */
  decomp(NN, A,  &cond, iwork, B);
  if (cond + 1 == cond){
    printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
    return 0;
  }

  /* 時間に関するループ */
  nMax = rint(Tmax / tau);
  for (n = 1; n <= nMax; n++) {
    /* 連立1次方程式の右辺を用意する */
    /* 内部の格子点 */ 
    for (i = 1; i < N_r; i++) {
      ri = (i)*h_r; 
      ri2 = ri*ri;
      for(j = 0; j<= N_p-1; j++){
        int jm1, jm2;
        L = psi(i,j);
        jm1 = j-1;
        jm2 = j+1;
        if (jm1 == -1)
          jm1 = N_p-1;
        if(jm2 == N_p)
          jm2 = 0;
        B[L] = (1.0 - 2.0*(1.0 - theta)*lambda_r
                    - (2.0*(1.0 - theta)*lambda_p)/ri2) * Uk[i][j]
          + (1.0 - theta) * lambda_r * ((1.0 + h_r/(2.0*ri)) *Uk[i+1][j]
          + (1.0 - h_r/(2.0*ri))*Uk[i-1][j])
          + ((1.0 - theta)*lambda_p / ri2)*(Uk[i][jm1] + Uk[i][jm2]);
      }
    }

    /* 原点 */
    B[0] = (1.0 - 4.0 * (1.0 - theta) * lambda_r) * Uk[0][0];
    for (j = 0; j <= N_p - 1; j++) {
      B[0] += (4.0*(1.0 - theta)*lambda_r / N_p)*Uk[1][j];
    }
   
    /* 原点のコピー */
    for (j = 1; j <= N_p-1; j++) {
      L = psi(0,j);
      B[L] = 0.0;
    }

    /* 境界条件 */
    for (j = 0; j <= N_p-1; j++) {
      L = psi(N_r,j);
      B[L] = 0.0;
    }

    /* A vector_U = B を解く */
    solve(NN, A, B, iwork);

    /* 更新 */
    for(i = 0; i<= N_r; i++){
      for(j = 0; j <= N_p-1; j++){
        L = psi(i,j);       
        Uk[i][j] = B[L];
      }
      Uk[i][N_p]= B[psi(i,0)];
    }
    t = n * tau;
 
    /* 誤差を測る */
    M = 0.0;
    for(i = 0; i <= N_r; i++){
      ri = (i)*h_r;
      for(j = 0; j <= N_p; j++){
        double e;
        phi_j = (j)*h_p;
        e = fabs(exactu(ri, phi_j, t) - Uk[i][j]);
        if(e > M)
          M = e;
      }
    }
    printf("n=%d, norm=%g, t=%g, 誤差=%g\n",
           n, maxnorm(N_r, N_p, Uk), t, M); 
    
    /* Δtの整数倍の時刻ではグラフをを描く */
    if (n % skip == 0){
      sprintf(label, "t=%g", t);
      disk(N_r, N_p, Uk, label);
    } 
  } 
  close_gnuplot();

  return 0;
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
  return exp(- mu01 * mu01 * t) * j0(mu01 * r)
    + exp(- mu11 * mu11 * t) * j1(mu11 * r) * (cos(phi) + sin(phi));
}

double maxnorm(int m, int n, matrix Uk)
{
  int i, j, i0, j0;
  double tmpmax, absu;
  i0 = 0;
  j0 = 0;
  tmpmax = fabs(Uk[0][0]);
  for (i = 0; i <= m; i++)
    for (j = 0; j <= n; j++)
      if ((absu = fabs(Uk[i][j])) > tmpmax) {
        tmpmax = absu;
        i0 = i;
        j0 = j;
      }
  printf("(i,j)=(%d,%d) ", i0, j0);
  return tmpmax;
}
