from itertools import filterfalse, chain
import fs.path
import csv
from rich.pretty import pprint

from tugboat.config import Config
from tugboat.log import Logger

class EmpiricalDataset:
    def __init__(self, collection: str, chromosome: str, basename: str, tsz_url: str, config: Config):
        self._collection = collection
        self._chromosome = chromosome
        self._basename = basename
        self._tsz_url = tsz_url
        self._config = config

    def collection(self):
        return self._collection

    def chromosome(self):
        return self._chromosome

    def basename(self):
        return self._basename

    def tsz_url(self):
        return self._tsz_url

    def tsz_file(self):
        return fs.path.combine(self._config.DATA_DIR, self.basename() + ".tsz")

    def trees_file(self):
        return fs.path.combine(self._config.DATA_DIR, self.basename() + ".trees")

    def forest_file(self):
        return fs.path.combine(self._config.DATA_DIR, self.basename() + ".forest")

    def ops_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".ops_bench.csv")

    def tajimasD_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".tajimasD_bench.csv")

    def conversion_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".conversion_bench.csv")

    _collection: str
    _chromosome: str
    _basename: str
    _tsz_url: str
    _config: Config

class SimulatedDataset:
    def __init__(self, collection: str, basename: str, config: Config):
        self._collection = collection
        self._basename = basename
        self._config = config

    def chromosome(self):
        return 'simulated'

    def tsz_url(self):
        return ''

    def collection(self):
        return self._collection

    def basename(self):
        return self._basename

    def trees_file(self):
        return fs.path.combine(self._config.DATA_DIR, self.basename() + ".trees")

    def forest_file(self):
        return fs.path.combine(self._config.DATA_DIR, self.basename() + ".forest")

    def ops_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".ops_bench.csv")

    def tajimasD_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".tajimasD_bench.csv")

    def conversion_bench_file(self):
        return fs.path.combine(self._config.MEASUREMENTS_DIR, self.basename() + ".conversion_bench.csv")

    _collection: str
    _basename: str
    _tsz_url: str
    _config: Config

class Datasets:
    @classmethod
    def from_csv(cls, filename_empirical, filename_scaling, config: Config, log: Logger):
        return cls(
            cls.empirical_from_csv(filename_empirical, config, log),
            cls.simulated_from_csv(filename_scaling, 'scaling', config, log)
        )

    @classmethod
    def empirical_from_csv(cls, filename: str, config: Config, log: Logger):
        try:
            with open(filename) as file:
                datasets = list()
                reader = csv.reader(file, delimiter=',', quotechar='"')
                row = reader.__next__()
                assert(row[0] == "collection")
                assert(row[1] == "chromosome")
                assert(row[2] == "basename")
                assert(row[3] == "tsz_url")
                for row in reader:
                    datasets.append(
                        EmpiricalDataset(
                            collection=row[0],
                            chromosome=row[1],
                            basename=row[2],
                            tsz_url=row[3],
                            config=config
                        )
                    )
            return datasets
        except OSError:
            log.critical(f"Cannot open datasets description: {config.DATASETS_CSV}")

    @classmethod
    def simulated_from_csv(cls, filename: str, collection: str, config: Config, log: Logger):
        try:
            with open(filename) as file:
                datasets = list()
                reader = csv.reader(file, delimiter=',', quotechar='"')
                row = reader.__next__()
                assert(row[0] == "name")
                for row in reader:
                    datasets.append(
                        SimulatedDataset(
                            collection=collection,
                            basename=f'simulated-{row[0]}',
                            config=config
                        )
                    )
            return datasets
        except OSError:
            log.critical(f"Cannot open datasets description: {config.DATASETS_CSV}")

    def __init__(self, *datasets):
        self._datasets = [ds for ds in chain(*datasets)]

    def all(self):
        return self._datasets

    def by_collection(self, collections):
        return filter(lambda ds: ds.collection() in collections, self._datasets)

    # TODO Add existing flag
    def trees_files(self):
        return [ds.trees_file() for ds in self._datasets]

    def forest_files(self):
        return [ds.forest_file() for ds in self._datasets]

    def ops_bench_files(self):
        return [ds.ops_bench_file() for ds in self._datasets]

    def conversion_bench_files(self):
        return [ds.conversion_bench_file() for ds in self._datasets]

    def tajimasD_bench_files(self):
        return [ds.tajimasD_bench_file() for ds in self._datasets]

    def tsz_files(self):
        return [ds.tsz_file() for ds in self._datasets]

    def tsz_urls(self):
        return [ds.tsz_url() for ds in self._datasets]
