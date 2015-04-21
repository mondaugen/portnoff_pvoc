#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <complex.h>
#include <fftw3.h> 
#include "portnoff_stft.h"
#include "window_calc.h" 
#include <math.h> 

static int input(fftw_complex *s, int N);
static int output(double *x, int R);

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr,"invoke with arguments N R Q where:\n"
                "N is the FFT-size of the spectra you are reading in\n"
                "R is the hop-size or interpolation factor\n"
                "Q is a measure of how many spectra to use."
                "The algorithm will use 2*Q spectra.");
        return(-1);
    }
    int N, R, Q, n, k, done;
    double *x, *f;
    fftw_complex *s_queue[2*Q];
    fftw_plan pb;
    N = atoi(argv[1]);
    R = atoi(argv[2]);
    Q = atoi(argv[3]);
    n = 0;
    x = (double*)malloc(sizeof(double)*R);
    f = (double*)malloc(sizeof(double)*(2*Q*R+1));
    /* Initialize window */
    hann_windowed_sinc(&f[Q*R],Q*R,R);
    /* Initialize queue of spectra to spectra containing 0s */
    for (k = 0; k < 2*Q; k++) {
        s_queue[k] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
        memset(s_queue[k],0,sizeof(fftw_complex)*N);
    }
    pb = fftw_plan_dft_1d(N,s_queue[2*Q-1],s_queue[2*Q-1],
            FFTW_BACKWARD,FFTW_ESTIMATE);
    done = 0;
    while (!done) {
        /* shift the spectrum queue */
        int l;
        fftw_complex *tmp;
        tmp = s_queue[0];
        for (l = 0; l < 2*Q-1; l++) {
            s_queue[l] = s_queue[l+1];
        }
        s_queue[l] = tmp;
        done = input(s_queue[l],N);
        if (done == -1) {
            fprintf(stderr,"Error reading file.\n");
        }
        /* make time-domain signal */
        fftw_execute_dft(pb,s_queue[l],s_queue[l]);
        /* scale values by 1/N */
        int m;
        for (m = 0; m < N; m++) {
            s_queue[l][m] /= (double)N;
        }
        /* set output buffer to 0 */
        memset(x,0,sizeof(double)*R);
        portnoff_synth_stream(x,&s_queue[Q-1],&f[Q*R],&n,N,R,Q);
        if (output(x,R)) {
            fprintf(stderr,"Error writing file.\n");
        }
    }
    fftw_destroy_plan(pb);
    for (k = 0; k < 2*Q; k++) {
        fftw_free(s_queue[k]);
    }
    free(x);
    free(f);
    return(0);
}

int input(fftw_complex *s, int N)
{
    int k, result = 0;
    k = 0;
    while (k < N) {
        if (fread(s,sizeof(fftw_complex),1,stdin) != 1) {
            if (feof(stdin)) {
                result = 1;
            } else {
                result = -1; /* error occured */
            }
            break;
        }
        s++;
        k++;
    }
    /* fill the rest with 0 */
    while (k < N) {
        *s = 0.;
        s++;
        k++;
    }
    return result;
}

int output(double *x, int R)
{
    if (fwrite(x,sizeof(double),R,stdout) != R) {
        return -1;
    }
    return 0;
}
