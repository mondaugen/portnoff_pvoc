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
    int N, P, H, H_cur, H_first, n = 0, r_hops, exit_code = 0, L_win;
    fftw_complex *X_n;
    double *x, *h, *h2;
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    H = atoi(argv[3]);
    L_win = 2*P*N+1;
    fftw_plan pf;
    X_n = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    x = (double*)malloc(sizeof(double)*(L_win));
    h = (double*)malloc(sizeof(double)*(L_win));
    hann_windowed_sinc(&h[P*N],P*N,N);
    memset(x,0,sizeof(double)*(L_win));
    pf = fftw_plan_dft_1d(N,X_n,X_n,FFTW_FORWARD,FFTW_ESTIMATE);
    r_hops = -1; /* the number of remaining hops, always -1 until the end */
    /* The first time we read in values, we want to read in enough so that the
     * centre of the buffer is aligned with the first sample of the input sound
     * file */
    H_first = P*N + 1;
    H_cur = H_first;
    while (r_hops) {
        int k, rem;
        memset(X_n,0,sizeof(fftw_complex)*N);
        /* Shift-in new values */
        for (k = 0; k < (L_win-H_cur); k++) {
            x[k] = x[k+H_cur];
        }
        /* If r_hops < 0, we read from stdin until the end of the file. If the
         * end of the file is reached, we set r_hops to the number of remaining
         * hops until the middle of the window has passed the end of the input
         * file */
        if (r_hops < 0) {
            /* If input returns a positive value it is the remaining number of
             * items to read in in order to obtain a full hop. */
            rem = input(&x[k],H_cur);
            if (rem == -1) {
                fprintf(stderr,"Error reading file.\n");
                exit_code = -1;
                goto cleanup;
            }
            if (rem > 0) {
                /* We have reached the end of the input file. This is how many
                 * additional hops we will do in order to get just past the last
                 * sample of the input file. */
                r_hops = (P*N+1)/H;
            }
        } else {
            r_hops--;
            rem = H;
        }
        /* Read in the missing samples as just 0. */
        for (k = L_win-rem; k < L_win; k++) {
            x[k] = 0.;
        }
        portnoff_analysis_stream(X_n,x+P*N,h+P*N,&n,N,P);
        n = (n+H)%N;
        fftw_execute(pf);
        if (output(X_n,N)) {
            fprintf(stderr,"Error writing file.\n");
            exit_code = -1;
            goto cleanup;
        }
        H_cur = H; /* After the first iteration, the hop size is always H */
    }
cleanup:
    fftw_destroy_plan(pf);
    fftw_free(X_n);
    free(x);
    free(h);
    return(exit_code);
}
    
/* When there's an error, return -1, when there are samples remaining but the
 * end of the file was reached, return a positive number indicating the number
 * of samples. Otherwise, return 0. */
int input(double *x, int H)
{
    int k, result = 0;
    k = 0;
    while (k < H) {
        if (fread(x,sizeof(double),1,stdin) != 1) {
            if (feof(stdin)) {
                return H - k;
            } else {
                return -1; /* error occured */
            }
        }
        x++;
        k++;
    }
    return 0; 
}

int output(fftw_complex *X, int N)
{
    if (fwrite(X,sizeof(fftw_complex),N,stdout) != N) {
        return -1;
    }
    return 0;
}
