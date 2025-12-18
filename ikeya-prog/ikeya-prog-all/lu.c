/*
 * lu.c -- decomp() and solve()
 *
 *         They are translations of FORTRAN subroutines "DECOMP" and 
 *         "SOLVE" in
 *           G.E.Forsythe, M.A.Malcolm and C.B.Moler,
 *           Computer Methods for Mathmatical Computations,
 *           Prentice Hall (1977).
 *
 * decomp() -- decomposes a real matrix by Gaussian elimination
 *             and estimates the condition of the matrix.
 * solve() -- solution of linear system, A*x = b
 *            do not use if decomp() ahs detected singularity
 */

#include <math.h>
#include "lu.h"

void decomp(int n, double **A, double *cond, int *ipvt, double *work)
{
/*
 *
 *    decomposes a real matrix by Gaussian elimination
 *    and estimates the condition of the matrix.
 *
 *    use solve() to compute solutions to linear systems.
 *
 *    input..
 *
 *      n = order of the matrix.
 *
 *      A = matrix to be triangularized.
 *
 *    output..
 *
 *      A  contains an upper triangular matrix   U  and a permuted
 *        version of a lower triangular matrix  I-L so that
 *       (permutation matrix)*A = L*U
 *
 *      cond = an estimate of the condition of A.
 *        For the linear system  A*x = b, changes in  A  and  b
 *        may cause changes  cond  times as large in  x.
 *        If  cond+1.0 == cond , A is singular to working
 *        precision.  cond  is set to  1.0E+32  if exact
 *        singularity is detected.
 *
 *      ipvt = the pivot vector.
 *        ipvt[k] = the index of the k-th pivot row
 *        ipvt[n-1] = (-1)**(number of interchanges)
 *
 *    work space..  The vector  work  must be declared and included
 *              in the call.  Its input contents are ignored.
 *              Its output contents are usually unimportant.
 *
 *    The determinant of  A  can be obtained on output by
 *      det(A) = ipvt[n-1] * A[0][0] * A[1][1] * ... * A[n-1][n-1]
 *
 */
    double ek, t, anorm, ynorm, znorm;
    int nm1, i, j, k, kp1, m;

    ipvt[n - 1] = 1;
    if (n == 1) {       /* 1-by-1 */
        if (A[0][0] == 0.0) /* exact singularity */
            *cond = 1.0E+32;
        else
            *cond = 1.0;
        return;
    }
    nm1 = n - 1;
    /* compute 1-norm of A */
    anorm = 0.0;
    for (j = 0; j < n; j++) {
        t = 0.0;
        for (i = 0; i < n; i++) t += fabs(A[i][j]);
        if (t > anorm)
            anorm = t;
    }
    /* Gaussian elimination with partial pivoting */
    for (k = 0; k < nm1; k++) {
        kp1 = k + 1;
        /* find pivot */
        m = k;
        for (i = kp1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[m][k]))
                m = i;
        ipvt[k] = m;
	if (m != k)
	  ipvt[n - 1] = - ipvt[n - 1];
        t = A[m][k]; 
        A[m][k] = A[k][k]; 
        A[k][k] = t;
        /* skip step if pivot is zero */
        if (t == 0.0)
            continue;
        /* compute multipliers */
        for (i = kp1; i < n; i++) A[i][k] /= (- t);
        /* interchange and eliminate by columns */
        for (j = kp1; j < n; j++) {
            t = A[m][j]; 
            A[m][j] = A[k][j]; 
            A[k][j] = t;
            if (t != 0.0)
                for (i = kp1; i < n; i++) A[i][j] += A[i][k] * t;
        }
    }

/*
 *    cond = (1-norm of A)*(an estimate of 1-norm of A-inverse)
 *    estimate obtained by one step of inverse iteration for the
 *    small singular vector. This involves solving two systems
 *    of equations,(A-transpose)*y = e  and  A*z = y  where  e
 *    is a vector of +1 or -1 chosen to cause growth in y.
 *    estimate = (1-norm of z)/(1-norm of y)
 */

    /* solve (A-transpose)*y = e */
    for (k = 0; k < n; k++) {
        t = 0.0;
        if (k != 0) {
            for (i = 0; i < k; i++) t += A[i][k] * work[i];
        }
        ek = 1.0;
        if (t < 0.0)
            ek = -1.0;
        if (A[k][k] == 0.0) { /* exact singularity */
            *cond = 1.0E+32;
            return;
        }
        work[k] = -(ek + t) / A[k][k];
    }
    for (k = nm1 - 1; k >= 0; k--) {
        t = 0.0;
        kp1 = k + 1;
        for (i = kp1; i < n; i++) t += A[i][k] * work[k];
        work[k] = t;
        m = ipvt[k];
        if (m != k) {
            t = work[m]; 
            work[m] = work[k]; 
            work[k] = t;
        }
    }

    ynorm = 0.0;
    for (i = 0; i < n; i++) ynorm += fabs(work[i]);
    /* solve A*z = y */
    solve(n, A, work, ipvt);

    znorm = 0.0;
    for (i = 0; i < n; i++) znorm += fabs(work[i]);
    /* estimate condition */
    *cond = anorm * znorm / ynorm;
    if (*cond < 1.0)
        *cond = 1.0;
}

void solve(int n, double **A, double *b, int *ipvt)
{
/*
 *  solution of linear system, A*x = b
 *  do not use if decomp() has detected singularity
 *
 *  input..
 *
 *    n = order of matrix.
 *
 *    A = triangularized matrix obtained from decomp()
 *
 *    b = right hand side vector.
 *
 *    ipvt = pivot vector obtained from decomp()
 *
 *  output..
 *
 *    b = solution vector, x .
 */
    int kb, nm1, kp1, i, k, m;
    double t;
    /* forward elimination */
    if (n == 1) {
        b[0] /= A[0][0];
        return;
    }
    nm1 = n - 1;
    for (k = 0; k < nm1; k++) {
        kp1 = k + 1;
        m = ipvt[k];
        t = b[m]; 
        b[m] = b[k]; 
        b[k] = t;
        for (i = kp1; i < n; i++) b[i] += A[i][k] * t;
    }
    /* back substitution */
    for (k = nm1; k > 0; k--) {
        b[k] /= A[k][k];
        t = - b[k];
        for (i = 0; i < k; i++) b[i] += A[i][k] * t;
    }
    b[0] /= A[0][0];
}
