#~/usr/bin/env bash

DATA_DIR="../sfkit/data"
MEASUREMENTS_DIR="."
SFKIT_BENCH="./build/release/benchmarks/sfkit-bench"
GIT_REV=$(git rev-parse --short HEAD)
MACHINE=$(hostname)
WARMUP_ITERS=1
ITERS=10

ls_files() {
    find "$DATA_DIR" -iname "*_chr*.trees" -print \
        | grep -E "(sgdp|1kg|unified)_*"
}

stats_cmd() {
    TREES_FILE="$1"
    STEM="$(basename "$TREES_FILE" .trees)"
    FOREST_FILE="${DATA_DIR}/${STEM}.forest"
    BPFOREST_FILE="${DATA_DIR}/${STEM}.bpforest"
    STATS_OUT="$MEASUREMENTS_DIR/${STEM}.stats.csv"

    echo "$SFKIT_BENCH stats --forest-file=$FOREST_FILE --trees-file=$TREES_FILE --bp-forest-file=$BPFOREST_FILE --revision=$GIT_REV --machine=$MACHINE > $STATS_OUT"
}

all_cmds() {
    for file in $(ls_files); do
        stats_cmd "$file"
    done
}

all_cmds

