#ifndef PORTNOFF_STFT_H
#define PORTNOFF_STFT_H 
#include <complex.h> 
void portnoff_analysis_stream(complex *x_n,
        double *x, double *h, int *n, int N, int P);
void portnoff_synth_stream(double *x, 
        complex **s, double *f, int *n, int N, int R, int Q);
#endif /* PORTNOFF_STFT_H */
