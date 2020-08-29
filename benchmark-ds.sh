#!/bin/bash

usage () {
    echo "usage: $0 [program] [first ID] [last ID] [numparties] [network file]"
}

if [ $# -lt 3 ]; then
    usage
    exit 0
fi

prog=$1
first_id=$2
last_id=$3
nparties=$4
nf=$5

echo "running with prog=$prog, network file=$nf"

if ! [ -d logs ]; then
    mkdir logs
fi



# for n in 1000 5000 10000 50000 100000 500000; do
for n in 1260 12600 63000 126000 630000; do

    for i in $(seq "${first_id}" 1 "${last_id}"); do
	echo "running $i ..."
    	./$prog $i $nf $n 2>> "logs/party_${i}_of_${nparties}_${n}_${prog}.txt" &
    	pid=$!
    done

    wait $pid

    sleep 2
    
done



echo "Done!"
