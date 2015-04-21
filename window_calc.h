#ifndef WINDOW_CALC_H
#define WINDOW_CALC_H 
void wc_sinc(double *x, double N, double R);
void wc_hann(double *x, double N);
void hann_windowed_sinc(double *x, int N, int R);
#endif /* WINDOW_CALC_H */
