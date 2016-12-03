#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <complex.h>
#include <fftw3.h> 
#include "portnoff_stft.h"
#include "window_calc.h" 
#include <math.h> 
#include <signal.h> 

static volatile int Ha;
static int maxHa;
static int H;
static double a;

static int input(double *x, int H);
static int output(fftw_complex *X, int N);
static void sig_slowdown(int signum);
static void sig_speedup(int signum);

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr,
                "arguments are:\n"
                "N - The FFT size\n"
                "P - The window size factor. Window size is 2*P*N + 1\n"
                "H - The hop size.\n"
                "a - The speed factor: a < 1 slower, a > 1 faster\n"
                "After, the speed can be set anywhere between 0 and 1 using SIGUSR1 and SIGUSR2. SIGUSR1 slows by 5%% of maximum and SIGUSR2 speeds up by 5%%.\n");
        return(-1);
    }
    int N, P, n, m, k, done;
    fftw_complex *X_n, *X_n_H, *Y_n;
    double *x, *h;
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    H = atoi(argv[3]);
    a = atof(argv[4]);
    maxHa = (int)(((double)H)*a);
    Ha = maxHa;
    fftw_plan pf;
    n = 0; /* index of "current" analysis frame */
    m = -H; /* index of analysis frame H samples ago */
    X_n = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    X_n_H = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    Y_n = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    /* Allocate enough buffer space for two analyses separated by a hop, plus
     * the extra space needed to hop for timestretching (probably you could get
     * away with less). */
    x = (double*)malloc(sizeof(double)*(2*P*N+H+1+Ha));
    h = (double*)malloc(sizeof(double)*(2*P*N+1));
    hann_windowed_sinc(&h[P*N],P*N,N);
    memset(x,0,sizeof(double)*(2*P*N+H+1));
    /* Initialize Y_n with small values */
    for (k = 0; k < N; k++) {
        Y_n[k] = 0.00001;
    }
    pf = fftw_plan_dft_1d(N,X_n,X_n,FFTW_FORWARD,FFTW_ESTIMATE);
    done = 0;
    signal(SIGUSR1,sig_slowdown);
    signal(SIGUSR2,sig_speedup);
    while (!done) {
        memset(X_n,0,sizeof(fftw_complex)*N);
        memset(X_n_H,0,sizeof(fftw_complex)*N);
        /* Shift-in new values */
        for (k = 0; k < (2*P*N+1-Ha); k++) {
            x[k] = x[k+Ha];
        }
        /* if input fails, set values to 0 and set done to 1 */
        done = input(&x[k],Ha);
        if (done == -1) {
            fprintf(stderr,"Error reading file.\n");
        }
        /* Analyse most current block */
        portnoff_analysis_stream(X_n,x+P*N+H,h+P*N,&n,N,P);
        n = (n+H)%N;
        /* Analyse block H samples ago */
        portnoff_analysis_stream(X_n_H,x+P*N,h+P*N,&m,N,P);
        m = (m+H)%N;
        fftw_execute_dft(pf,X_n,X_n);
        fftw_execute_dft(pf,X_n_H,X_n_H);
        /* Calculate new Y_n */
        for (k = 0; k < N; k++) {
            /* avoid divide by 0 */
            X_n[k] += 10E-6;
            X_n_H[k] += 10E-6;
            Y_n[k] = X_n[k]*(Y_n[k]/X_n_H[k])/cabs(Y_n[k]/X_n_H[k]);
        }
        if (output(Y_n,N)) {
            fprintf(stderr,"Error writing file.\n");
        }
    }
    fftw_destroy_plan(pf);
    fftw_free(X_n);
    fftw_free(X_n_H);
    fftw_free(Y_n);
    free(x);
    free(h);
    return(0);
}
    
int input(double *x, int H)
{
    int k, result = 0;
    k = 0;
    while (k < H) {
        if (fread(x,sizeof(double),1,stdin) != 1) {
            if (feof(stdin)) {
                result = 1;
            } else {
                result = -1; /* error occured */
            }
            break;
        }
        x++;
        k++;
    }
    /* fill the rest with 0 */
    while (k < H) {
        *x = 0.;
        x++;
        k++;
    }
    return result;
}

int output(fftw_complex *X, int N)
{
    if (fwrite(X,sizeof(fftw_complex),N,stdout) != N) {
        return -1;
    }
    return 0;
}

void sig_slowdown(int signum)
{
    if ((a -= 0.01) < 0) {
        a = 0;
    }
    Ha = (int)(((double)H)*a); 
}

void sig_speedup(int signum)
{
    if ((a += 0.01)*H > maxHa) {
        a = ((double)maxHa)/H;
    }
    Ha = (int)(((double)H)*a); 
}
