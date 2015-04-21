if [[ ! -e bin ]]
then
    mkdir bin
fi
gcc pstft_time_stretch.c ../*.c -I.. -o bin/pstft_time_stretch -lfftw3 -lm
