#include <stdio.h> 
#include "window_calc.h" 

int main (int argc, char **argv) {
    if (argc !=2) {
        fprintf(stderr,"Give value for N as argument\n");
        return(-1);
    }
    int N;
    int n;
    N = atoi(argv[1]);
    double x[2*N+1];
    wc_hann(&x[N],N);
    fwrite(x,sizeof(double),2*N+1,stdout);
    return(0);
}
