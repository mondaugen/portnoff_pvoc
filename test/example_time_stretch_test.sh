#!/bin/bash
if [[ -z $1 ]]
then
    echo "give speed as argument"
    exit 65
fi
sox -n -t f64 -c 1 -r 44100 - synth 1 sine 400 | \
    bin/pstft_time_stretch 512 4 128 $1 | \
    bin/pstft_synth_test 512 128 4 | sox -t f64 -c 1 -r 44100 - -d
