#!/usr/bin/env python3

import numpy
import tskit
from time import perf_counter
from os import path
import argparse

def s_to_ns(seconds):
    return seconds * 10**9

argparser = argparse.ArgumentParser(
                    prog='benchmark-tskits-python-only-funcs',
                    description='Benchmarks tskits Tajima\'s D (which is Python only)')
argparser.add_argument('-n', '--iterations', type=int, help='The number of times to run each benchmark', default=1)
argparser.add_argument('-f', '--trees-file', type=str, help='The tree sequence file to benchmark on. If not specified.', required=True)
argparser.add_argument('-r', '--revision', type=str, help='Revision of the benchmarked program (unique id, e.g. git commit hash)', default='undefined')
argparser.add_argument('-m', '--machine', type=str, help='Identifier of this computer (e.g. hostname)', default='undefined')

args = argparser.parse_args()
num_iterations = args.iterations
ts_file = args.trees_file
revision = args.revision
machine_id = args.machine

print('section,variant,dataset,revision,machine_id,iteration,variable,value,unit')
dataset = path.basename(ts_file)

ts = tskit.load(ts_file)

for iteration in range(0, num_iterations):
    start = perf_counter()
    tajimas_d = ts.Tajimas_D()
    end = perf_counter()

    walltime = s_to_ns(end - start)
    print(f'tajimas_d,tskit,{dataset},{revision},{machine_id},{iteration},walltime,{walltime},ns')
    print(f'tajimas_d,tskit,{dataset},{revision},{machine_id},{iteration},tajimas_d,{tajimas_d},1')

for iteration in range(0, num_iterations):
    sample_set_1 = numpy.array(ts.samples()[1::2], dtype = numpy.int32)
    sample_set_2 = numpy.array(ts.samples()[0::2], dtype = numpy.int32)

    start = perf_counter()
    fst = ts.Fst((sample_set_1, sample_set_2))
    end = perf_counter()

    walltime = s_to_ns(end - start)
    print(f'fst,tskit,{dataset},{revision},{machine_id},{iteration},walltime,{walltime},ns')
    print(f'fst,tskit,{dataset},{revision},{machine_id},{iteration},fst,{fst},1')

for iteration in range(0, num_iterations):
    for nth in range(2, 10):
        samples = ts.samples()[0::nth]

        start = perf_counter()
        for tree in ts.trees():
            lca = tree.mrca(*samples)
        end = perf_counter()

        walltime = s_to_ns(end - start)
        print(f'lca_{nth}th,tskit,{dataset},{revision},{machine_id},{iteration},walltime,{walltime},ns')
