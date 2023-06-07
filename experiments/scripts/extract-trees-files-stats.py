#!/usr/bin/env python3

import argparse
import tskit
import re
import os


class TSKitMetadata():
    def __init__(self, filename):
        self._filename = filename
        self._ts = tskit.load(filename)

    def filename(self):
        return self._filename

    def num_trees(self):
        return self._ts.num_trees

    def sequence_length(self):
        return self._ts.sequence_length

    def num_sample_nodes(self):
        return self._ts.num_samples

    def total_size(self):
        return self._ts.nbytes

    def num_edges(self):
        return self._ts.num_edges

    def num_individuals(self):
        return self._ts.num_individuals

    def num_migrations(self):
        return self._ts.num_migrations

    def num_mutations(self):
        return self._ts.num_mutations

    def num_nodes(self):
        return self._ts.num_nodes

    def num_populations(self):
        return self._ts.num_populations

    def num_provenances(self):
        return self._ts.num_provenances

    def num_sites(self):
        return self._ts.num_sites


class FilenameMetadata():
    def __init__(self, filename, organism, chromosome, collection):
        self._filename = filename
        self._organism = organism
        self._chromosome = chromosome
        self._collection = collection

    def filename(self):
        return self._filename

    def organism(self):
        return self._organism

    def chromosome(self):
        return self._chromosome

    def collection(self):
        return self._collection

    _filename = None
    _chromosome = None
    _organism = None
    _collection = None


class FilenameMetadataFactory():
    @classmethod
    def build(cls, filename):
        filename = os.path.basename(filename)
        metadata = FilenameMetadataFactory._parse_with_first_successful(filename,
                                                                        cls._parse_tgp_filename,
                                                                        cls._parse_1kg_filename,
                                                                        cls._parse_sgdp_filename,
                                                                        cls._parse_unified_filename,
                                                                        cls._parse_anderson_filename)
        if metadata is None:
            raise ValueError("Unknown filename format: " + filename)
        else:
            return metadata

    @staticmethod
    def _parse_with_first_successful(filename, *parsers):
        for parser in parsers:
            metadata = parser(filename)
            if metadata is not None:
                return metadata
        return None

    @staticmethod
    def _parse_tgp_filename(filename):
        match = re.match(r'tgp_chr(\d+)\.trees', filename)
        if match is None:
            return None
        else:
            return FilenameMetadata(filename=filename, organism='human', chromosome=match.group(
                1), collection='1000 Genomes Project (our inference)')

    @staticmethod
    def _parse_1kg_filename(filename):
        match = re.match(r'1kg_chr(\d+)\.trees', filename)
        if match is None:
            return None
        else:
            return FilenameMetadata(filename=filename, organism='human', chromosome=match.group(
                1), collection='1000 Genomes Project')

    @staticmethod
    def _parse_sgdp_filename(filename):
        match = re.match(r'sgdp_chr(\d+)\.trees', filename)
        if match is None:
            return None
        else:
            return FilenameMetadata(filename=filename, organism='human', chromosome=match.group(
                1), collection='Simons Genome Diversity Project')

    @staticmethod
    def _parse_unified_filename(filename):
        match = re.match(r'unified_chr(\d+)\.trees', filename)
        if match is None:
            return None
        else:
            return FilenameMetadata(filename=filename, organism='human', chromosome=match.group(
                1), collection='Unified (Wohns et al.)')

    @staticmethod
    def _parse_anderson_filename(filename):
        match = re.match(r'anderson_chr(\d+)\.trees', filename)
        if match is None:
            return None
        else:
            return FilenameMetadata(filename=filename, organism='human', chromosome=match.group(
                1), collection='Anderson-Trocmé (simulated)')


def print_stats(filename_metadata, tskit_metadata):
    print(filename_metadata.filename(), filename_metadata.collection(), filename_metadata.organism(), filename_metadata.chromosome(), tskit_metadata.num_trees(), tskit_metadata.sequence_length(), tskit_metadata.num_sample_nodes(), tskit_metadata.total_size(
    ), tskit_metadata.num_edges(), tskit_metadata.num_individuals(), tskit_metadata.num_migrations(), tskit_metadata.num_mutations(), tskit_metadata.num_nodes(), tskit_metadata.num_populations(), tskit_metadata.num_provenances(), tskit_metadata.num_sites(), sep=',')


def main():
    # Parse command line arguments
    arg_parser = argparse.ArgumentParser(
        prog='extract-stats-from-trees-files.py',
        description='Uses tskit to extract stats from .trees files and outputs them in csv format.')

    arg_parser.add_argument('filename', nargs='+')

    args = arg_parser.parse_args()

    # Print CSV header
    print('filename,collection,organism,chromosome,num_trees,sequence_length,num_sample_nodes,total_size,num_edges,num_individuals,num_migrations,num_mutations,num_nodes,num_populations,num_provenances,num_sites')

    # Iterate over all files and print the statistics for each file.
    for filename in args.filename:
        filename_metadata = FilenameMetadataFactory.build(filename)
        tskit_metadata = TSKitMetadata(filename)
        print_stats(filename_metadata, tskit_metadata)


if __name__ == "__main__":
    main()
