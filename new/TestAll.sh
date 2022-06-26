#!/bin/bash
src="data/english.txt"
dest="data/dna.txt"
bus=("4" "2" "1" "0")
modes=("1" "2")
multiples=("2" "4")
rm cram_test
make cram_test
for mode in ${modes[@]}; do
    for multiple in ${multiples[@]}; do
        for bu in ${bus[@]}; do
            ./cram_test "$src" "$dest" 1 "${bu}" "${mode}" "${multiple}"
        done
    done
done
