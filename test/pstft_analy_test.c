#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <complex.h>
#include <fftw3.h> 
#include "portnoff_stft.h"
#include "window_calc.h" 
#include <math.h> 


static int input(double *x, int H);
static int output(fftw_complex *X, int N);
static void hann_windowed_sinc_(double *x, int N, int R);
static void write_hann_windowed_sinc(double *x, int N, char *name);

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr,
                "arguments are:\n"
                "N - The FFT size\n"
                "P - The window size factor. Window size is 2*P*N + 1\n"
                "H - The hop size.\n");
        return(-1);
    }
    int N, P, H, n, done;
    fftw_complex *X_n;
    double *x, *h, *h2;
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    H = atoi(argv[3]);
    fftw_plan pf;
    X_n = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    x = (double*)malloc(sizeof(double)*(2*P*N+1));
    h = (double*)malloc(sizeof(double)*(2*P*N+1));
    h2 = (double*)malloc(sizeof(double)*(2*P*N+1));
    hann_windowed_sinc_(&h[P*N],P*N,N);
    hann_windowed_sinc(&h2[P*N],P*N,N);
    write_hann_windowed_sinc(h,2*P*N+1,"WTF_better.raw");
    write_hann_windowed_sinc(h2,2*P*N+1,"WTF_worse.raw");
    memset(x,0,sizeof(double)*(2*P*N+1));
    pf = fftw_plan_dft_1d(N,X_n,X_n,FFTW_FORWARD,FFTW_ESTIMATE);
    done = 0;
    while (!done) {
        int k;
        memset(X_n,0,sizeof(fftw_complex)*N);
        /* Shift-in new values */
        for (k = 0; k < (2*P*N+1-H); k++) {
            x[k] = x[k+H];
        }
        /* if input fails, set values to 0 and set done to 1 */
        done = input(&x[k],H);
        if (done == -1) {
            fprintf(stderr,"Error reading file.\n");
        }
        portnoff_analysis_stream(X_n,x+P*N,h+P*N,&n,N,P);
        n = (n+H)%N;
        fftw_execute(pf);
        if (output(X_n,N)) {
            fprintf(stderr,"Error writing file.\n");
        }
    }
    fftw_destroy_plan(pf);
    fftw_free(X_n);
    free(x);
    free(h);
    free(h2);
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

/* calculates a window of length 2*N+1 that is 0 every R samples. x is the address
 * of the centre of the space allocated to store the window */
static void hann_windowed_sinc_(double *x, int N, int R)
{
    double w[2*N+1];
    /* The commented out configuration sounds worse when invoking:
     * sox -n -t f64 -c 1 -r 44100 - synth 10 sine 400-1000 | \
     * ./pstft_analy_test.elf 512 4 128 | \
     * ./pstft_synth_test.elf 512 128 4 | \
     * sox -t f64 -c 1 -r 44100 - -d
     * WTF? */
//    double *w;
//    w = (double*)malloc(sizeof(double)*(2*N+1));
    int n;
    wc_hann(&w[N],N);
    wc_sinc(x,N,R);
    for (n = -N; n <= N; n++) {
        x[n] *= w[n+N];
    }
//    free(w);
}

static void write_hann_windowed_sinc(double *x, int N, char *name)
{
    FILE *f = fopen(name,"w");
    fwrite(x,sizeof(double),N,f);
    fclose(f);
}
