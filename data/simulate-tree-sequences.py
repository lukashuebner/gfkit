import msprime

OUTPUT_DIR: str = "."
BASE_SEED: int = 1337
DEFAULT_SEQUENCE_LENGTH: int = int(1e+8)
DEFAULT_MUTATION_RATE: float = 1e-8
DEFAULT_NUM_INDIVIDUALS: int = 100
DEFAULT_POPULATION_SIZE: int = DEFAULT_NUM_INDIVIDUALS * 100
DEFAULT_RECOMBINATION_RATE: float = 1e-8

# --- Single population model ---


def single_population_model(pop_size: int, seq_len: int, mut_rate: float, num_individuals: int, recombination_rate: float, seed: int):
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
    ts.dump(f"{OUTPUT_DIR}/single-population-model-popsize{pop_size}-seqlen{seq_len}-mutrate{mut_rate}-numindv{num_individuals}.trees"
            )

for population_size in [1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9]:
    single_population_model(
        pop_size=int(population_size),
        seq_len=DEFAULT_SEQUENCE_LENGTH,
        mut_rate=DEFAULT_MUTATION_RATE,
        num_individuals=DEFAULT_NUM_INDIVIDUALS,
        recombination_rate=DEFAULT_RECOMBINATION_RATE,
        seed=BASE_SEED
    )

for sequence_length in [1e+6, 1e+7, 1e+8, 1e+9, 1e+10]:
    single_population_model(
        pop_size=DEFAULT_POPULATION_SIZE,
        seq_len=int(sequence_length),
        mut_rate=DEFAULT_MUTATION_RATE,
        num_individuals=DEFAULT_NUM_INDIVIDUALS,
        recombination_rate=DEFAULT_RECOMBINATION_RATE,
        seed=BASE_SEED
    )

for mutation_rate in [1e-8, 1e-7, 1e-6]:
    single_population_model(
        pop_size=DEFAULT_POPULATION_SIZE,
        seq_len=DEFAULT_SEQUENCE_LENGTH,
        mut_rate=mutation_rate,
        num_individuals=DEFAULT_NUM_INDIVIDUALS,
        recombination_rate=DEFAULT_RECOMBINATION_RATE,
        seed=BASE_SEED
    )

for num_individuals in [10, 100, 1000, 10000]:
    single_population_model(
        pop_size=DEFAULT_POPULATION_SIZE,
        seq_len=DEFAULT_SEQUENCE_LENGTH,
        mut_rate=DEFAULT_MUTATION_RATE,
        num_individuals=num_individuals,
        recombination_rate=DEFAULT_RECOMBINATION_RATE,
        seed=BASE_SEED
    )

for recombination_rate in [1e-8, 1e-7, 1e-6]:
    single_population_model(
        pop_size=DEFAULT_POPULATION_SIZE,
        seq_len=DEFAULT_SEQUENCE_LENGTH,
        mut_rate=DEFAULT_MUTATION_RATE,
        num_individuals=DEFAULT_NUM_INDIVIDUALS,
        recombination_rate=recombination_rate,
        seed=BASE_SEED
    )

# TODO Add sequence simulation
# --- Model with two populations ---


def two_population_model(pop1_size: int, pop2_size: int, ancestral_size: int, split_time: int, num_samples_per_population: int, recombination_rate: float, sequence_length: int, seed: int, filename: str):
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
    ts.dump(filename)
    return ts

# --- Model with four populations ---


def simulation_loop(pop1_size: int, pop2_size: int, pop3_size: int, pop4_size: int,
                    pop34_size: int, pop234_size: int, pop1234_size: int,
                    split_time_34: int, split_time_234: int, split_time_1234: int,
                    migration_proportion: int, recombination_rate: float,
                    sequence_length: int, seed: int, filename: str):
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
    ts.dump(filename)

# Run a loop that saves treesequences in files
# reps = 10 ## Number of repetition
# for i in range(reps):
#     ## 1st line in for loop is modern sizes
#     ## 2nd line is ancestral sizes
#     ## 3rd line is split times
#     temp_ts = simulation_loop( 1000, 1000, 1000, 1000,
#                                100, 100, 100,
#                                1000, 2000, 3000,
#                                0.05, 50 )
#     temp_ts.dump( "path/to/save/location/name".join(i) )
