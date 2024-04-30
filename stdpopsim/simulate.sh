#!/usr/bin/env bash

JOB_QUEUE_FILE="__job-queue.sh"
JOB_LOG_FILE="simulate.log"
NUM_SAMPLES=(5 10 20 40 80 160 320 640 12800 25600 51200 102400 204800 409600)
PARALLEL="$HOME/bin/parallel"

NUM_CORES=$(grep "^cpu\\scores" /proc/cpuinfo | uniq |  awk '{print $4}')
EXPECTED_MAX_MEM="40G"

# Empty the output file
echo -n "" > "$JOB_QUEUE_FILE"

# Generate job queue
for n in ${NUM_SAMPLES[@]}; do
    # for chr in {1..22}; do
    for chr in 20; do
        echo "[ -f simulated_${n}k_chr${chr}.trees ] || stdpopsim HomSap --genetic-map HapMapII_GRCh38 --chromosome $chr --output simulated_${n}k_chr${chr}.trees ${n}000" >> "$JOB_QUEUE_FILE";
    done
done

# Run the jobs
cat "$JOB_QUEUE_FILE" \
    | sort -n -k 11 \
    | sort -n -k 14 --stable \
    | "$PARALLEL" --jobs="$NUM_CORES" --joblog="$JOB_LOG_FILE" --memsuspend="$EXPECTED_MAX_MEM"

