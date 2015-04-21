sox -n -t f64 -c 1 -r 44100 - synth 1 sine 400 | bin/pstft_analy_test 512 4 128 | bin/pstft_synth_test.elf 512 128 4 | sox -t f64 -c 1 -r 44100 - -d
