#!/usr/bin/env bash

DATASETS_CSV='data/scaling-datasets.csv'

mlr --csv \
    put '
        $filename = "scaling-" . $name . ".trees";
        $cmd_msp_ancestry = "msp ancestry --random-seed " . $seed . " --length " . $sequence_length . " --ploidy 2 --population-size " . $population_size . " --recombination-rate " . $recombination_rate . " " . $num_individuals;
        $cmd_msp_mutations = "msp mutations --model jc69 --random-seed " . $seed . " --output " . $filename . " " . $mutation_rate;
        $cmd = "[[ -f " . $filename . " ]] || bash -c \"" . $cmd_msp_ancestry . " | " . $cmd_msp_mutations . "\"";
    ' \
    "$DATASETS_CSV" \
    | mlr --icsv --opprint --headerless-csv-output cut -f cmd \
    | parallel --memsuspend 50G
