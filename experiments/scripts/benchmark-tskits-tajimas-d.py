#!/usr/bin/env python3

import tskit
from time import perf_counter
from os import path
import argparse

def s_to_ns(seconds):
    return seconds * 10**9

argparser = argparse.ArgumentParser(
                    prog='tskit-tajimas-d-benchmark',
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
    ts.Tajimas_D()
    end = perf_counter()

    walltime = s_to_ns(end - start)
    print(f'tajimas_d,tskit,{dataset},{revision},{machine_id},{iteration},walltime,{walltime},ns')
