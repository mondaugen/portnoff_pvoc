#include <stdio.h> 
#include "window_calc.h" 

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
    wc_sinc(&x[N],N,R);
    fwrite(x,sizeof(double),2*N+1,stdout);
    return(0);
}
