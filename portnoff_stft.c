#include <complex.h> 

/* x is the destination to which the output will be written. x[k] must be
 * defined for 0 <= k < R.
 * s contains 2*Q pointers to arrays of complex double values of length N. If s_ is the
 * beginning of the space allocated to store the 2*Q pointers then &s_[Q-1]
 * should be passed to the function as it uses negative indexing.
 * f is a 1:R interpolating window of length 2*Q*R + 1. Pass a pointer to f_[QR]
 * where f_ is the beginning of the space holding the window.
 * n is pointer to an int keeping track of the current position in s.
 * N, R and Q are the values as described above.
 */
void portnoff_synth_stream(double *x, 
        complex double **s, double *f, int *n, int N, int R, int Q)
{
    int Lmin, Lmax, count;
    Lmin = -Q + 1;
    Lmax = Q;
    count = R;
    while (count--) {
        int r;
        for (r = Lmin; r <= Lmax; r++) {
            *x += f[-r*R]*creal(s[r][*n]);
        }
        *n = (*n+1)%N;
        x++;
        f++;
    }
}

/* x_n is an array of complex double numbers of length N. Typically this is initialized
 * to containing only zeroes.
 * x is the input signal where x[k] is defined for -P*N <= k < P*N
 * h is a pointer to the middle of a window. Typically h(k) = 0 for k = l*N,
 * l = (-P...-1,1...P), h(0) = 1 (see Portnoff (1976) for details).
 * n is a pointer to the sample index to which x refers. Typically this is
 * incremented by H after calling this function, where H is the hop-size, and is
 * reduced modulo N, i.e., *n = (*n + H)%N
 * N, P are as described above.
 */
void portnoff_analysis_stream(complex double *x_n,
        double *x, double *h, int *n, int N, int P)
{
    int m, l;
    for (m = 0; m < N; m++) {
        for (l = -P; l < P; l++) {
            x_n[(m+*n)%N] += x[l*N+m]*h[-l*N-m];
        }
    }
}

