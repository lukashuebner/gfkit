#!/usr/bin/env bash

DATASETS_CSV='./scaling-datasets.csv'
MLR="$(which mlr)"
PARALLEL="$HOME/bin/parallel"

"$MLR" --csv \
    put '
        $filename = "scaling-" . $name . ".trees";
        $cmd_msp_ancestry = "msp ancestry --random-seed " . $seed . " --length " . $sequence_length . " --ploidy 2 --population-size " . $population_size . " --recombination-rate " . $recombination_rate . " " . $num_individuals;
        $cmd_msp_mutations = "msp mutations --model jc69 --random-seed " . $seed . " --output " . $filename . " " . $mutation_rate;
        $cmd = "[[ -f " . $filename . " ]] || bash -c \"" . $cmd_msp_ancestry . " | " . $cmd_msp_mutations . "\"";
    ' \
    "$DATASETS_CSV" \
    | "$MLR" --icsv --opprint cut -f cmd \
    | sed "1d" \
    | "$PARALLEL"
