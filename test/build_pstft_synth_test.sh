if [[ ! -e bin ]]
then
    mkdir bin
fi
gcc pstft_synth_test.c ../*.c -I.. -o bin/pstft_synth_test -lfftw3 -lm -ggdb3
