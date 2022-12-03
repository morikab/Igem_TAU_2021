import statistics
import typing
from dataclasses import dataclass
from enum import Enum


class Organism(object):
    CAI_PROFILE_ATTRIBUTE_NAME = "cai_profile"
    TAI_PROFILE_ATTRIBUTE_NAME = "tai_profile"

    def __init__(self,
                 name: str,
                 cai_profile,
                 tai_profile,
                 cai_scores: typing.Dict,
                 tai_scores: typing.Dict,
                 is_optimized: bool,
                 optimization_priority: float):
        self.name = name
        self.cai_profile = cai_profile
        self.tai_profile = tai_profile
        self.cai_scores = cai_scores
        self._cai_scores_values = cai_scores.values()
        self.tai_scores = tai_scores
        self._tai_scores_values = tai_scores.values()
        self.is_optimized = is_optimized
        self.optimization_priority = optimization_priority

    @property
    def cai_avg(self):
        return statistics.mean(self._cai_scores_values)

    @property
    def tai_avg(self):
        return statistics.mean(self._tai_scores_values)

    @property
    def cai_std(self):
        return statistics.stdev(self._cai_scores_values)

    @property
    def tai_std(self):
        return statistics.stdev(self._tai_scores_values)


class OptimizationMethod(Enum):
    single_codon_global_ratio = "single_codon_global_ratio"
    single_codon_local_ratio = "single_codon_local_ratio"
    single_codon_global_diff = "single_codon_global_diff"
    single_codon_local_diff = "single_codon_local_diff"
    hill_climbing_average = "hill_climbing_average"
    hill_climbing_bulk_aa_average = "hill_climbing_bulk_aa_average"
    hill_climbing_weakest_link = "hill_climbing_weakest_link"

    @property
    def is_single_codon_optimization(self) -> bool:
        return self in (OptimizationMethod.single_codon_global_ratio,
                        OptimizationMethod.single_codon_local_ratio,
                        OptimizationMethod.single_codon_global_diff,
                        OptimizationMethod.single_codon_local_diff)


class OptimizationCubScore(Enum):
    codon_adaptation_index = "CAI"
    trna_adaptation_index = "tAI"
    max_codon_trna_adaptation_index = "max_CAI_tAI"

    @property
    def is_codon_adaptation_score(self) -> bool:
        return self in (OptimizationCubScore.codon_adaptation_index,
                        OptimizationCubScore.max_codon_trna_adaptation_index)

    @property
    def is_trna_adaptation_score(self) -> bool:
        return self in (OptimizationCubScore.trna_adaptation_index,
                        OptimizationCubScore.max_codon_trna_adaptation_index)


@dataclass
class UserInput:
    organisms: typing.List[Organism]
    sequence: str
    zip_directory: str
    tuning_parameter: float
    clusters_count: int
    optimization_method: OptimizationMethod = OptimizationMethod.hill_climbing_bulk_aa_average
    optimization_cub_score: OptimizationCubScore = OptimizationCubScore.max_codon_trna_adaptation_index
