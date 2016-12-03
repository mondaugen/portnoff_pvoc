#!/bin/bash
if [[ -z $1 ]]
then
    echo "give file as 1st argument"
    exit 65
fi
if [[ -z $2 ]]
then
    echo "give speed as 2nd argument"
    exit 65
fi
if [[ -z $3 ]]
then
    echo "give original speed as 3rd argument"
    exit 65
fi
mkfifo /tmp/fifo.$$
sox $1 -t f64 -c 1 -r 44100 - speed $3 | \
    bin/pstft_time_stretch 2048 4 512 $2 > /tmp/fifo.$$&
echo "$!" > /tmp/tsproc;
cat /tmp/fifo.$$ | bin/pstft_synth_test 2048 512 4 | \
    sox -t f64 -c 1 -r 44100 - -d

# Adjust speed with
# kill -30 `cat /tmp/tsproc` # to slow down and
# kill -31 `cat /tmp/tsproc` # to speed up
