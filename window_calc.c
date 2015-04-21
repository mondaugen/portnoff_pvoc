#include <math.h> 
#include <stdlib.h> 

/* Generate a window of length 2*N + 1 which is a sinc function whose zeros are
 * at x = lR for floor(-N/R) <= l <= floor(N/R)
 * x should point to the datum x_[N] as the algorithm uses negative indexing to
 * compute a symmetrical window.
 * R is as described above.
 */
void wc_sinc(double *x, double N, double R)
{
    double n;
    for (n = -N; n <= N; n += 1.) {
        x[(int)n] = (n == 0.) ? 1. : sin(M_PI*n/R) / (M_PI*n/R);
    }
}

/* Generate a Hann window of length 2*N + 1.
 * x should point to datum x_[N]
 */
void wc_hann(double *x, double N)
{
    double n;
    for (n = -N; n <= N; n += 1.) {
        x[(int)n] = 0.5*(cos(M_PI*n/N) + 1.);
    }
}

/* calculates a window of length 2*N+1 that is 0 every R samples. x is the address
 * of the centre of the space allocated to store the window */
void hann_windowed_sinc(double *x, int N, int R)
{
//    double w[2*N+1];
    double *w = (double*)malloc(sizeof(double)*(2*N+1));
    int n;
    wc_hann(&w[N],N);
    wc_sinc(x,N,R);
    for (n = -N; n <= N; n++) {
        x[n] *= w[n+N];
    }
    free(w);
}

/* calculates a window of length 2*N+1 that is 0 every R samples. x is the address
 * of the centre of the space allocated to store the window */
/*
static void hann_windowed_sinc(double *x, int N, int R)
{
    double w[2*N+1];
    int n;
    wc_hann(&w[N],N);
    wc_sinc(x,N,R);
    for (n = -N; n <= N; n++) {
        x[n] *= w[n+N];
    }
}
*/
