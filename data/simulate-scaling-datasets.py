#!/usr/bin/env python3

import csv
import msprime
import concurrent.futures

CONFIGS_FILE: str = "scaling-datasets.csv"
OUTPUT_DIR: str = "."

# --- Single population model ---
def single_population_model(pop_size: int, seq_len: int, mut_rate: float, num_individuals: int, recombination_rate: float, seed: int, filename: str):
    # These model parameters are taken from: https://tskit.dev/tutorials/getting_started.html
    sweep_model = msprime.SweepGenicSelection(
        position=seq_len / 2,
        start_frequency=0.0001,
        end_frequency=0.9999,
        s=0.25,
        dt=1e-6
    )

    ts = msprime.sim_ancestry(
        num_individuals,
        model=[sweep_model, msprime.StandardCoalescent()],
        population_size=pop_size,
        sequence_length=seq_len,
        recombination_rate=recombination_rate,
        random_seed=seed
    )

    ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)
    ts.dump(filename)

with concurrent.futures.ProcessPoolExecutor(100) as executor:
    with open(CONFIGS_FILE) as configs_file:
        configs_reader = csv.reader(configs_file)
        first_row = True
        for dataset in configs_reader:
            if first_row:
                first_row = False
                continue
            name = dataset[0]
            pop_size = float(dataset[1])
            seq_len = float(dataset[2])
            mut_rate = float(dataset[3])
            num_indiv = float(dataset[4])
            recomb_rate = float(dataset[5])
            seed = int(dataset[6])
            file = f'scaling-{name}.trees'
            executor.submit(
                single_population_model,
                pop_size=pop_size,
                seq_len=seq_len,
                mut_rate=mut_rate,
                num_individuals=num_indiv,
                recombination_rate=recomb_rate,
                seed=seed,
                filename=file
            )

