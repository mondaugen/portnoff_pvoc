if [[ ! -e bin ]]
then
    mkdir bin
fi
gcc -I/usr/local/include pstft_analy_test.c ../*.c -I.. \
    -o bin/pstft_analy_test -lfftw3 -lm
    
