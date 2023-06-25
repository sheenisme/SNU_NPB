#!/bin/bash

array=( 0 10 20 40 60 80 85 90 95 96 97 98 99 100)
for r in ${array[@]}
do
    make CLASS=A
    echo "ppcg --target c -R ${r} --no-reschedule ssor.c"
    ppcg --target c -R ${r} --no-reschedule ssor.c
    make CLASS=A
    taskset -c 0 ../bin/lu.A.x
    taskset -c 0 ../bin/lu.A.x 
    taskset -c 0 ../bin/lu.A.x 
    taskset -c 0 ../bin/lu.A.x 
    taskset -c 0 ../bin/lu.A.x 
    echo "${r}  over"
    echo ""
done

# ppcg --target c -R 99 --no-reschedule ssor.c