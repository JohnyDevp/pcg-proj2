#!/bin/bash

#
# @file      runFinalBenchmark.sh
#
# @author    David Bayer \n
#            Faculty of Information Technology \n
#            Brno University of Technology \n
#            ibayer@fit.vutbr.cz
#
# @brief     PCG Assignment 1
#
# @version   2024
#
# @date      04 October   2023, 09:00 (created) \n
#

if [[ $# < 1 ]]; then
    echo "Usage: <path to nbody binary>"
    exit -1
fi

BENCHMARK_DIR=Benchmark
INPUT_DIR=$BENCHMARK_DIR/Inputs
OUTPUT_DIR=$BENCHMARK_DIR/Outputs

NBODY_BINARY=$1
NBODY_BINARY_NAME=$(basename $NBODY_BINARY)
GEN_BINARY=$(dirname $NBODY_BINARY)/gen

DT=0.01f
STEPS=100
THREADS_PER_BLOCK=512
WRITE_INTENSITY=10
RED_THREADS=4096
RED_THREADS_PER_BLOCK=128

mkdir -p $BENCHMARK_DIR
mkdir -p $INPUT_DIR
mkdir -p $OUTPUT_DIR

printf "   N\t  Time\n"

for i in {1..8}; do
    N=$((2**i * 512))

    INPUT=$INPUT_DIR/Input$N.h5
    OUTPUT=$OUTPUT_DIR/${NBODY_BINARY_NAME}Output${N}.h5

    if [ ! -f $INPUT ]; then
    $GEN_BINARY $N $INPUT >> /dev/null
    fi
    
    ncu --set full --launch-count 1 -o prof-$N -f "$NBODY_BINARY" $N $DT $STEPS $WRITE_INTENSITY "$INPUT" "$OUTPUT"
done