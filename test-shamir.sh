#!/bin/bash

usage () {
    echo "usage: $0 [number of parties] [network file]"
}

if [ $# -lt 2 ]; then
    usage
    exit 0
fi

n=$1
nf=$2
prog=shamir_test

log_dir=test_log

if ! [ -d $log_dir ]; then
    mkdir $log_dir
fi

if [ $n -lt 4 ]; then
    echo "n must be at least 4"
fi

make $prog

if [ $? -ne 0 ]; then
    echo "make failed"
    exit 0
fi

echo "running with prog=$prog, n=$n, network file=$nf"

for i in $(seq 0 $((n - 1))); do
    ./${prog}.exe $i $nf > "${log_dir}/party_${i}_of_${n}.txt" &
    pid=$!
done

wait $pid

echo "Done!"
