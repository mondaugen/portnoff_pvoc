if [[ ! -e bin ]]
then
    mkdir bin
fi
gcc -I/usr/local/include -L/usr/local/lib pstft_synth_test.c ../*.c -I.. \
    -o bin/pstft_synth_test -lfftw3 -lm -g -O0
