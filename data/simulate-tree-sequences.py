#!/usr/bin/env python3

import msprime
import concurrent.futures 

RANDOM_NAMES_FILE: str = "random-names.txt"
OUTPUT_DIR: str = "."
BASE_SEED: int = 1337
DEFAULT_SEQUENCE_LENGTH: int = int(1e+8)
DEFAULT_MUTATION_RATE: float = 1e-8
DEFAULT_NUM_INDIVIDUALS: int = 100
DEFAULT_POPULATION_SIZE: int = DEFAULT_NUM_INDIVIDUALS * 100
DEFAULT_RECOMBINATION_RATE: float = 1e-8
DEFAULT_SPLIT_TIME: int =  1000
DEFAULT_MIGRATION_PROPORTION: float = 0.5

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

# --- Model with two populations ---
def two_population_model(pop1_size: int, pop2_size: int, ancestral_size: int, split_time: int, num_samples_per_population: int, recombination_rate: float, sequence_length: int, mut_rate: float, seed: int, filename: str):
    # Building demographic scenario
    demography = msprime.Demography()
    demography.add_population(name="ANCESTRAL",
                              initial_size=ancestral_size)
    demography.add_population(name="Pop1",
                              initial_size=pop1_size,
                              default_sampling_time=0)
    demography.add_population(name="Pop2",
                              initial_size=pop2_size,
                              default_sampling_time=0)
    demography.add_population_split(time=split_time,
                                    derived=["Pop1", "Pop2"],
                                    ancestral="ANCESTRAL")

    # ALWAYS RUN "sort_events" starting the simulation!
    demography.sort_events()
    ts = msprime.sim_ancestry(
        demography=demography,
        # Specify number of samples from each population
        samples={"Pop1": num_samples_per_population,
                 "Pop2": num_samples_per_population},
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=seed
    )
    ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)
    ts.dump(filename)
    return ts

# --- Model with four populations ---
def four_population_model(pop1_size: int, pop2_size: int, pop3_size: int, pop4_size: int,
                    pop34_size: int, pop234_size: int, pop1234_size: int,
                    split_time_34: int, split_time_234: int, split_time_1234: int,
                    migration_proportion: int, recombination_rate: float,
                    sequence_length: int, mut_rate: float, seed: int, filename: str):
    # Building demographic scenario
    # In msprime you need to treat internal branches as ancestral populations.
    # So we need to specify
    demography = msprime.Demography()
    demography.add_population(name="anc_1234",
                              initial_size=pop1234_size)
    demography.add_population(name="anc_234",
                              initial_size=pop234_size)
    demography.add_population(name="anc_34",
                              initial_size=pop34_size)
    demography.add_population(name="Pop1",
                              initial_size=pop1_size,
                              default_sampling_time=0)
    demography.add_population(name="Pop2",
                              initial_size=pop2_size,
                              default_sampling_time=0)
    demography.add_population(name="Pop3",
                              initial_size=pop3_size,
                              default_sampling_time=0)
    demography.add_population(name="Pop4",
                              initial_size=pop4_size,
                              default_sampling_time=0)

    demography.add_population_split(time=split_time_34,
                                    derived=["Pop3", "Pop4"],
                                    ancestral="anc34")
    demography.add_population_split(time=split_time_234,
                                    derived=["Pop2", "anc_34"],
                                    ancestral="anc_234")
    demography.add_population_split(time=split_time_1234,
                                    derived=["Pop1", "anc_234"],
                                    ancestral="anc_1234")

    # Command to add migration. Note that migrations is backwards in time
    # that means that source-destination is decided as you move backwards
    # in terms of traversing a tree.
    demography.add_migration_rate_change(
        time=0,
        source="Pop4",
        dest="Pop1",
        rate=migration_proportion
    )
    # ALWAYS RUN "sort_events" starting the simulation!
    demography.sort_events()
    ts = msprime.sim_ancestry(
        demography=demography,
        samples={"Pop1": 5, "Pop2": 5,
                 "Pop3": 5, "Pop4": 5},  # Specify number of samples from each population. There are smarter ways to sample; see manual on sampling.
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=seed
    )
    ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)
    ts.dump(filename)

with concurrent.futures.ProcessPoolExecutor() as executor:
    with open(RANDOM_NAMES_FILE) as file:
        names = [line.rstrip() for line in file]
    name = iter(names)

    # Single population model
    with open("single-population-datasets.csv", "w") as datasets_csv:
        datasets_csv.write("name,population_size,sequence_length,mutation_rate,num_individuals,recombination_rate,seed\n")
        for population_size in [1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9]:
            file = next(name)
            datasets_csv.write(f"{file},{population_size},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
            executor.submit(
                single_population_model,
                pop_size=int(population_size),
                seq_len=DEFAULT_SEQUENCE_LENGTH,
                mut_rate=DEFAULT_MUTATION_RATE,
                num_individuals=DEFAULT_NUM_INDIVIDUALS,
                recombination_rate=DEFAULT_RECOMBINATION_RATE,
                seed=BASE_SEED,
                filename = file
            )

        for sequence_length in [1e+6, 1e+7, 1e+8, 1e+9, 1e+10]:
            file = next(name)
            datasets_csv.write(f"{file},{DEFAULT_POPULATION_SIZE},{sequence_length},{DEFAULT_MUTATION_RATE},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
            executor.submit(
                single_population_model,
                pop_size=DEFAULT_POPULATION_SIZE,
                seq_len=int(sequence_length),
                mut_rate=DEFAULT_MUTATION_RATE,
                num_individuals=DEFAULT_NUM_INDIVIDUALS,
                recombination_rate=DEFAULT_RECOMBINATION_RATE,
                seed=BASE_SEED,
                filename = file
            )

        for mutation_rate in [1e-8, 1e-7, 1e-6]:
            file = next(name)
            datasets_csv.write(f"{file},{DEFAULT_POPULATION_SIZE},{DEFAULT_SEQUENCE_LENGTH},{mutation_rate},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
            executor.submit(
                single_population_model,
                pop_size=DEFAULT_POPULATION_SIZE,
                seq_len=DEFAULT_SEQUENCE_LENGTH,
                mut_rate=mutation_rate,
                num_individuals=DEFAULT_NUM_INDIVIDUALS,
                recombination_rate=DEFAULT_RECOMBINATION_RATE,
                seed=BASE_SEED,
                filename = file
            )

        for num_individuals in [10, 100, 1000, 10000]:
            file = next(name)
            datasets_csv.write(f"{file},{DEFAULT_POPULATION_SIZE},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{num_individuals},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
            executor.submit(
                single_population_model,
                pop_size=DEFAULT_POPULATION_SIZE,
                seq_len=DEFAULT_SEQUENCE_LENGTH,
                mut_rate=DEFAULT_MUTATION_RATE,
                num_individuals=num_individuals,
                recombination_rate=DEFAULT_RECOMBINATION_RATE,
                seed=BASE_SEED,
                filename = file
            )

        for recombination_rate in [1e-8, 1e-7, 1e-6]:
            file = next(name)
            datasets_csv.write(f"{file},{DEFAULT_POPULATION_SIZE},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{DEFAULT_NUM_INDIVIDUALS},{recombination_rate},{BASE_SEED}\n")
            executor.submit(
                single_population_model,
                pop_size=DEFAULT_POPULATION_SIZE,
                seq_len=DEFAULT_SEQUENCE_LENGTH,
                mut_rate=DEFAULT_MUTATION_RATE,
                num_individuals=DEFAULT_NUM_INDIVIDUALS,
                recombination_rate=recombination_rate,
                seed=BASE_SEED,
                filename = file
            )

    # Two population model
    with open("two-population-datasets.csv", "w") as datasets_csv:
        datasets_csv.write("name,pop1_size,pop2_size,ancestral_size,split_time,num_samples_per_population,sequence_length,mutation_rate,recombination_rate,seed\n")
        for population_size in [1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9]:
            for pop1_fraction in [0.1, 0.2, 0.3, 0.4, 0.5]:
                pop2_fraction: float = 1 - pop1_fraction
                pop1_size = int(pop1_fraction * population_size)
                pop2_size = int(pop2_fraction * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{DEFAULT_SPLIT_TIME},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=int(0.5 * population_size),
                    split_time=DEFAULT_SPLIT_TIME,
                    num_samples_per_population=DEFAULT_NUM_INDIVIDUALS,
                    recombination_rate=DEFAULT_RECOMBINATION_RATE,
                    sequence_length=DEFAULT_SEQUENCE_LENGTH,
                    mut_rate=DEFAULT_MUTATION_RATE,
                    seed=BASE_SEED,
                    filename = file
                )

            for split_time in [100, 1000, 10000, 100000]:
                pop1_size = int(0.5 * population_size)
                pop2_size = int(0.5 * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{split_time},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=ancestral_size,
                    split_time=split_time,
                    num_samples_per_population=DEFAULT_NUM_INDIVIDUALS,
                    recombination_rate=DEFAULT_RECOMBINATION_RATE,
                    sequence_length=DEFAULT_SEQUENCE_LENGTH,
                    mut_rate=DEFAULT_MUTATION_RATE,
                    seed=BASE_SEED,
                    filename = file
                )
                
            for num_samples_per_population in [10, 100, 1000, 10000]:
                pop1_size = int(0.5 * population_size)
                pop2_size = int(0.5 * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{DEFAULT_SPLIT_TIME},{num_samples_per_population},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=ancestral_size,
                    split_time=DEFAULT_SPLIT_TIME,
                    num_samples_per_population=num_samples_per_population,
                    recombination_rate=DEFAULT_RECOMBINATION_RATE,
                    sequence_length=DEFAULT_SEQUENCE_LENGTH,
                    mut_rate=DEFAULT_MUTATION_RATE,
                    seed=BASE_SEED,
                    filename = file
                )

            for recombination_rate in [1e-8, 1e-7, 1e-6]:
                pop1_size = int(0.5 * population_size)
                pop2_size = int(0.5 * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{DEFAULT_SPLIT_TIME},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{recombination_rate},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=ancestral_size,
                    split_time=DEFAULT_SPLIT_TIME,
                    num_samples_per_population=DEFAULT_NUM_INDIVIDUALS,
                    recombination_rate=recombination_rate,
                    sequence_length=DEFAULT_SEQUENCE_LENGTH,
                    mut_rate=DEFAULT_MUTATION_RATE,
                    seed=BASE_SEED,
                    filename = file
                )

            for sequence_length in [1e+6, 1e+7, 1e+8, 1e+9, 1e+10]:
                pop1_size = int(0.5 * population_size)
                pop2_size = int(0.5 * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{DEFAULT_SPLIT_TIME},{DEFAULT_NUM_INDIVIDUALS},{sequence_length},{DEFAULT_MUTATION_RATE},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=ancestral_size,
                    split_time=DEFAULT_SPLIT_TIME,
                    num_samples_per_population=DEFAULT_NUM_INDIVIDUALS,
                    recombination_rate=DEFAULT_RECOMBINATION_RATE,
                    sequence_length=int(sequence_length),
                    mut_rate=DEFAULT_MUTATION_RATE,
                    seed=BASE_SEED,
                    filename = file
                )
            
            for mut_rate in [1e-8, 1e-7, 1e-6]:
                pop1_size = int(0.5 * population_size)
                pop2_size = int(0.5 * population_size)
                ancestral_size = int(0.5 * population_size)
                file = next(name)
                datasets_csv.write(f"{file},{pop1_size},{pop2_size},{ancestral_size},{DEFAULT_SPLIT_TIME},{DEFAULT_NUM_INDIVIDUALS},{DEFAULT_SEQUENCE_LENGTH},{mut_rate},{DEFAULT_RECOMBINATION_RATE},{BASE_SEED}\n")
                executor.submit(
                    two_population_model,
                    pop1_size=pop1_size,
                    pop2_size=pop2_size,
                    ancestral_size=ancestral_size,
                    split_time=DEFAULT_SPLIT_TIME,
                    num_samples_per_population=DEFAULT_NUM_INDIVIDUALS,
                    recombination_rate=DEFAULT_RECOMBINATION_RATE,
                    sequence_length=DEFAULT_SEQUENCE_LENGTH,
                    mut_rate=mut_rate,
                    seed=BASE_SEED,
                    filename = file
                )

    # Four population model
    with open("four-population-datasets.csv", "w") as datasets_csv:
        datasets_csv.write("dataset,pop1_size,pop2_size,pop3_size,pop4_size,pop34_size, pop234_size, pop1234_size,split_time_34,split_time_234,split_time_1234,migration_proportion,recombination_rate,sequence_length,mutation_rate,seed\n")
        for pop1_fraction in [0.1, 0.2, 0.3, 0.4, 0.5]:
            for pop2_fraction in [0.1, 0.2, 0.3, 0.4, 0.5]:
                for pop3_fraction in [0.1, 0.2, 0.3, 0.4, 0.5]:
                    if pop1_fraction >= 1.0: continue
                    else:
                        pop4_fraction: float = 1 - pop1_fraction - pop2_fraction - pop3_fraction
                        pop1_size = int(pop1_fraction * DEFAULT_POPULATION_SIZE)
                        pop2_size = int(pop2_fraction * DEFAULT_POPULATION_SIZE)
                        pop3_size = int(pop3_fraction * DEFAULT_POPULATION_SIZE)
                        pop4_size = int(pop4_fraction * DEFAULT_POPULATION_SIZE)
                        pop34_size = int((pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
                        pop234_size = int((pop2_fraction + pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
                        pop1234_size = int((pop1_fraction + pop2_fraction + pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
                        file = next(name)
                        datasets_csv.write(f"{file},{pop1_size},{pop2_size},{pop3_size},{pop4_size},{pop34_size},{pop234_size},{pop1234_size},{DEFAULT_SPLIT_TIME},{DEFAULT_SPLIT_TIME},{DEFAULT_SPLIT_TIME},{DEFAULT_MIGRATION_PROPORTION},{DEFAULT_RECOMBINATION_RATE},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{BASE_SEED}\n")
                        executor.submit(
                            four_population_model,
                            pop1_size = pop1_size,
                            pop2_size = pop2_size,
                            pop3_size = pop3_size,
                            pop4_size = pop4_size,
                            pop34_size = pop34_size,
                            pop234_size = pop234_size,
                            pop1234_size = pop1234_size,
                            split_time_34 = DEFAULT_SPLIT_TIME,
                            split_time_234 = DEFAULT_SPLIT_TIME,
                            split_time_1234 = DEFAULT_SPLIT_TIME,
                            migration_proportion = DEFAULT_MIGRATION_PROPORTION,
                            recombination_rate = DEFAULT_RECOMBINATION_RATE,
                            sequence_length = DEFAULT_SEQUENCE_LENGTH,
                            mut_rate = DEFAULT_MUTATION_RATE,
                            seed = BASE_SEED,
                            filename = file
                        )

        for migration_proportion in [0.1, 0.2, 0.3, 0.4, 0.5]:
            pop1_size = int(pop1_fraction * DEFAULT_POPULATION_SIZE)
            pop2_size = int(pop2_fraction * DEFAULT_POPULATION_SIZE)
            pop3_size = int(pop3_fraction * DEFAULT_POPULATION_SIZE)
            pop4_size = int(pop4_fraction * DEFAULT_POPULATION_SIZE)
            pop34_size = int((pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
            pop234_size = int((pop2_fraction + pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
            pop1234_size = int((pop1_fraction + pop2_fraction + pop3_fraction + pop4_fraction) * DEFAULT_POPULATION_SIZE)
            file = next(name)
            datasets_csv.write(f"{file},{pop1_size},{pop2_size},{pop3_size},{pop4_size},{pop34_size},{pop234_size},{pop1234_size},{DEFAULT_SPLIT_TIME},{DEFAULT_SPLIT_TIME},{DEFAULT_SPLIT_TIME},{migration_proportion},{DEFAULT_RECOMBINATION_RATE},{DEFAULT_SEQUENCE_LENGTH},{DEFAULT_MUTATION_RATE},{BASE_SEED}\n")
            executor.submit(
                four_population_model,
                pop1_size = pop1_size,
                pop2_size = pop2_size,
                pop3_size = pop3_size,
                pop4_size = pop4_size,
                pop34_size = pop34_size,
                pop234_size = pop234_size,
                pop1234_size = pop1234_size,
                split_time_34 = DEFAULT_SPLIT_TIME,
                split_time_234 = DEFAULT_SPLIT_TIME,
                split_time_1234 = DEFAULT_SPLIT_TIME,
                migration_proportion = migration_proportion,
                recombination_rate = DEFAULT_RECOMBINATION_RATE,
                sequence_length = DEFAULT_SEQUENCE_LENGTH,
                mut_rate = DEFAULT_MUTATION_RATE,
                seed = BASE_SEED,
                filename = file
            )
