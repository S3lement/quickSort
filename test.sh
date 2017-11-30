#!/usr/bin/env bash
for ((i=1;i<=10;i++)); do
    mpirun -np 8 ./a.out $i 1000 | sort -c -g
done

