// c_SPSumInverse.c: This file contains the low-level c function

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

static double square(const double h) { return h*h; }
#define ABORT(exitcode,msg,...){ fprintf(stderr,"Error %d: " #msg "\n",exitcode,__VA_ARGS__); \
                                    return exitcode; }

/// <summary>Calculates the inverse sum of two single-pair matrices A+C</summary>
/// <param name="n">Dimension</param>
/// <param name="tol">Tolerance for zero tests</param>
/// <param name="oa">Generator a[1]...a[n] of A (array of size n+1)</param>
/// <param name="xb">Parameter b[0] = x and generator b[1]...b[n] of A (array of size n+1)</param>
/// <param name="oc">Generator c[1]...c[n] of C (array of size n+1)</param>
/// <param name="Q">Output upper-triangular inverse q(i,j) = Q[(2*n-i)*(i-1)/2+j-1] (array of size n(n+1)/2)</param>
/// <param name="Pi">Output products of beta coefficients (array of size n(n+1)/2-1)</param>
/// <param name="zv">Parameter v[0] = z,output continuants v[1]...v[n] (array of size n+1)</param>
/// <param name="obeta">Output beta coefficients (array of size n)</param>
/// <param name="verbose">Flag for warning output</param>
/// <returns>0 if successful,<0 if error,>0 if warning </returns>
#define q(i,j) Q[(2*n-(i))*(i-1)/2+j-1]
#define pi(i,j) Pi[(2*n+1-(i))*(i-2)/2+j-1]
#define a oa
#define b xb
#define c oc
#define v zv
#define beta obeta
int c_SPSumInverse(int n,double tol,double oa[],double xb[],double oc[],
    double Q[],double Pi[],double zv[],double obeta[],_Bool verbose)
{
    static int i, j, k, warning = 0;
    //Check arguments
    if (!(a && b && c && Q && Pi && v && beta))
        ABORT(-1, "An array pointer is NULL");
    if (n < 3 || !(isfinite(tol) && isfinite(b[0]) && isfinite(v[0])) || tol <= 0. || fabs(v[0]) < tol)
        ABORT(-2, "Invalid n = %d,tol = %e,x = b[0] = %e,or z = v[0] = %e", n, tol, b[0], v[0]);
    if (fabs(b[1] - b[0]) < tol)
        ABORT(-3, "Invalid x = b[0] = %e too close to b[1] = %e", b[0], b[1]);
    //Initialize
    a[0] = c[0] = beta[0] = 0.;
    v[1] = (v[0] * (a[1] * b[1] + c[1])) / square(b[1] - b[0]);
    if (fabs(v[1]) < tol)
        ABORT(-4, "Invalid x = b[0],z = v[0],a[1],b[1],or c[1] (v[1] = %e ~ 0)", v[1]);
    beta[1] = (a[1] * b[0] + c[1]) / square(b[1] - b[0]);
    //Compute continuants and betas
    for (i = 2; i <= n; i++) {
        if (fabs(b[i] - b[i - 1]) < tol || fabs(b[i - 1] - b[i - 2]) < tol)
            ABORT(-16, "Some consecutive b[] values are too close around index %d", i);
        v[i] = ((a[i] * b[i] - 2 * a[i - 1] * b[i] + a[i - 1] * b[i - 1] + c[i] - c[i - 1]) / square(b[i] - b[i - 1]) -
            (a[i - 1] * b[i - 1] - 2 * a[i - 1] * b[i - 2] + a[i - 2] * b[i - 2] - (c[i - 1] - c[i - 2])) / square(b[i - 1] - b[i - 2]))
            * v[i - 1] - square((a[i - 1] * b[i - 2] - a[i - 2] * b[i - 1] + c[i - 1] - c[i - 2]) /
                square(b[i - 1] - b[i - 2])) * v[i - 2];
        if (fabs(v[i]) < tol)
            ABORT(-17, "Low continuant error (v[%d] = %e ~0)", i, v[i]);
        beta[i] = (a[i] * b[i - 1] - a[i - 1] * b[i] - c[i - 1] + c[i]) / square(b[i] - b[i - 1]);
        if (fabs(v[i] - beta[i] * v[i - 1]) < tol) {
            if (verbose) fprintf(stderr, "Warning: continuant value v[%d] ~ beta[%d]v[%d] might cause \
    unreliable results.\n", i, i, i - 1);
            warning = 18;
        }
        pi(i, i) = beta[i];
        pi(i, i - 1) = 1.;
    }
    //Compute products of betas
    for (i = 2; i <= n; i++)
        for (k = i + 1; k <= n; k++)
            pi(i, k) = pi(i, k - 1) * beta[k];
    //Compute inverse
    for (j = 1; j <= n; j++) {
        q(j, j) = 0.;
        if (j < n - 1) {
            for (k = j + 2; k <= n; k++)
                q(j, j) += square(pi(j + 2, k - 1)) / (v[k] * v[k - 1]);
            q(j, j) *= square((v[j + 1] - beta[j + 1] * v[j]) / (b[j + 1] - b[j]) - (beta[j + 1] * (v[j] - beta[j] * v[j - 1]))
                / (b[j] - b[j - 1]));
        }
        q(j, j) += v[j - 1] / (square(b[j] - b[j - 1]) * v[j]);
        if (j < n)
            q(j, j) += square((v[j] - beta[j] * v[j - 1]) / (b[j] - b[j - 1]) + v[j] / (b[j + 1] - b[j])) / (v[j + 1] * v[j]);
        for (i = 1; i <= j - 1; i++) {
            q(i, j) = 0.;
            if (j < n - 1) {
                for (k = j + 2; k <= n; k++)
                    q(i, j) += (pi(i + 2, k - 1) * pi(j + 2, k - 1)) / (v[k] * v[k - 1]);
                q(i, j) *= ((v[i + 1] - beta[i + 1] * v[i]) / (b[i + 1] - b[i]) -
                    (beta[i + 1] * (v[i] - beta[i] * v[i - 1])) / (b[i] - b[i - 1])) *
                    ((v[j + 1] - beta[j + 1] * v[j]) / (b[j + 1] - b[j]) - (beta[j + 1] * (v[j] - beta[j] * v[j - 1]))
                        / (b[j] - b[j - 1]));
            }
            if (j < n)
                q(i, j) -= (((v[i + 1] - beta[i + 1] * v[i]) / (b[i + 1] - b[i]) - (beta[i + 1] * (v[i] - beta[i] * v[i - 1]))
                    / (b[i] - b[i - 1])) * ((v[j] - beta[j] * v[j - 1]) / (b[j] - b[j - 1]) + v[j] / (b[j + 1] - b[j]))
                    * pi(i + 2, j)) / (v[j + 1] * v[j]);
            if (i < j - 1)
                q(i, j) += (((v[i + 1] - beta[i + 1] * v[i]) / (b[i + 1] - b[i]) - (beta[i + 1] * (v[i] - beta[i] * v[i - 1]))
                    / (b[i] - b[i - 1])) * pi(i + 2, j - 1)) / ((b[j] - b[j - 1]) * v[j]);
            else // i = j - 1
                q(i, j) -= ((v[j - 1] - beta[j - 1] * v[j - 2]) / (b[j - 1] - b[j - 2]) + v[j - 1] / (b[j] - b[j - 1]))
                / ((b[j] - b[j - 1]) * v[j]);
        }
    }
    return warning;
}