import msprime
import argparse
import sys

arg_parser = argparse.ArgumentParser(
    prog = "generate-tree-sequence",
    description = "Generates a tree sequence by simulation with msprime."
)
arg_parser.add_argument('population_size', type=int)
arg_parser.add_argument('sequence_length', type=int)
arg_parser.add_argument('num_individuals', type=int)
arg_parser.add_argument('seed', default=1337, type=int)
arg_parser.add_argument('--add-mutations', action='store_true')
arg_parser.add_argument('--mutation-rate', default=1e-8, type=float)
args = arg_parser.parse_args()

# These model parameters are taken from: https://tskit.dev/tutorials/getting_started.html
sweep_model = msprime.SweepGenicSelection(
    position = args.sequence_length / 2,
    start_frequency = 0.0001,
    end_frequency = 0.9999,
    s = 0.25,
    dt = 1e-6
)

# Simualte tree sequence
ts = msprime.sim_ancestry(
    args.num_individuals,
    model = [sweep_model, msprime.StandardCoalescent()],
    population_size = args.population_size,
    sequence_length = args.sequence_length,
    recombination_rate = 1e-8,
    random_seed = args.seed
)

# Optionally add finite-site mutations to the ts using the Jukes & Cantor model, creating SNPs
if (args.add_mutations):
    ts = msprime.sim_mutations(ts, rate=args.mutation_rate, random_seed=args.seed)

ts.dump(sys.stdout)
