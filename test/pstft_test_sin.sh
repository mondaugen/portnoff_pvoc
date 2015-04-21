#!/bin/bash

E_BADARGS=65
ARGS=5

if [ $# -ne "$ARGS" ]
then
    echo 'Usage: '`basename $0`' FFT-size P-size hop-size length freq'
    echo 'all sizes in samples, freq is dependent on FFT-size (0 - FFT-size)'
    exit $E_BADARGS
fi

# BUG octave doesn't persist though
octave -qf --eval \
    'fwrite(stdout,sin('${5}'*(((1:'${4}')-1)/'${1}')*2*pi),"double")' | \
    ./pstft_analy_test $1 $2 $3 | \
    octave -qf --persist --eval \
    'X=fread(stdin,Inf,"double");
     Xr=X(1:2:end);
     Xi=X(2:2:end);
     Xr=reshape(Xr,['${1}' length(Xr)/'${1}']);
     Xi=reshape(Xi,['${1}' length(Xi)/'${1}']);
     Xm=sqrt(Xr.^2+Xi.^2);
     surf(Xm);'&
