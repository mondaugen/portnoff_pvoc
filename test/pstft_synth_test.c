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
    int N, R, Q, n, k, rem, r_reads, exit_code = 0;
    double *x, *f;
    fftw_plan pb;
    N = atoi(argv[1]);
    R = atoi(argv[2]);
    Q = atoi(argv[3]);
    fftw_complex *s_queue[2*Q];
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
    /* r_reads is -1 until the end of the input file is reached, then at the end
     * it is set to the number of remaining reads */
    r_reads = -1; 
    rem = 0;
    while (r_reads) {
        /* shift the spectrum queue */
        int l, k;
        fftw_complex *tmp;
        tmp = s_queue[0];
        for (l = 0; l < 2*Q-1; l++) {
            s_queue[l] = s_queue[l+1];
        }
        s_queue[l] = tmp;
        if (r_reads < 0) {
            rem = input(s_queue[l],N);
            if (rem == -1) {
                fprintf(stderr,"Error reading file.\n");
                exit_code = -1;
                goto cleanup;
            }
            if (rem > 0) {
                /* We have reached the end of the input file. This is how many
                 * additional reads we want to make before the true end */
                r_reads = 2*Q;
            }
        } else {
            r_reads--;
            rem = N;
        }
        /* Read in the missing samples as just 0. */
        for (k = N-rem; k < N; k++) {
            s_queue[l][k] = 0.;
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
            exit_code = -1;
            goto cleanup;
        }
    }
cleanup:
    fftw_destroy_plan(pb);
    for (k = 0; k < 2*Q; k++) {
        fftw_free(s_queue[k]);
    }
    free(x);
    free(f);
    return(exit_code);
}

int input(fftw_complex *s, int N)
{
    int k;
    k = 0;
    while (k < N) {
        if (fread(s,sizeof(fftw_complex),1,stdin) != 1) {
            if (feof(stdin)) {
                return N - k;
            } else {
                return -1; /* error occured */
            }
        }
        s++;
        k++;
    }
    return 0;
}

int output(double *x, int R)
{
    if (fwrite(x,sizeof(double),R,stdout) != R) {
        return -1;
    }
    return 0;
}
