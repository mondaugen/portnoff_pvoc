/* Comparing two windowed sinc functions. Why would they differ? */
#include <stdio.h> 
#include "window_calc.h" 

/* calculates a window of length 2*N+1 that is 0 every R samples. x is the address
 * of the centre of the space allocated to store the window */
void win_calc(double *x, int N, int R)
{
    double w[2*N+1];
    int n;
    wc_hann(&w[N],N);
    wc_sinc(x,N,R);
    for (n = -N; n <= N; n++) {
        x[n] *= w[n+N];
    }
}

int main (int argc, char **argv) {
    if (argc !=3 ) {
        fprintf(stderr,"Give values for N and R as arguments\n");
        return(-1);
    }
    int N, R;
    int n;
    N = atoi(argv[1]);
    R = atoi(argv[2]);
    double x[2*N+1];
    double y[2*N+1];
    win_calc(&x[N],N,R);
    hann_windowed_sinc(&y[N],N,R);
    for (n = 0; n < 2*N+1; n++) {
        if (x[n] != y[n]) {
            fprintf(stderr,"x[%d] and y[%d] differ by %f\n",n,n,x[n]-y[n]);
        }
    }
    fwrite(x,sizeof(double),2*N+1,stdout);
    return(0);
}

