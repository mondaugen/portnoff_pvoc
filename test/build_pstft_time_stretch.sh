if [[ ! -e bin ]]
then
    mkdir bin
fi
gcc -I/usr/local/include -L/usr/local/lib pstft_time_stretch.c ../*.c -I.. \
    -o bin/pstft_time_stretch -lfftw3 -lm
