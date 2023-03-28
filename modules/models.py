import statistics
import typing
from dataclasses import dataclass
from enum import Enum


class Organism(object):
    CAI_PROFILE_ATTRIBUTE_NAME = "cai_profile"
    TAI_PROFILE_ATTRIBUTE_NAME = "tai_profile"

    def __init__(self,
                 name: str,
                 is_optimized: bool,
                 optimization_priority: float,
                 cai_profile: typing.Dict = None,
                 tai_profile: typing.Dict = None,
                 cai_scores: typing.Dict = None,
                 tai_scores: typing.Dict = None):
        self.name = name
        self.cai_profile = cai_profile
        self.tai_profile = tai_profile
        self.cai_scores = cai_scores
        self._cai_scores_values = cai_scores.values() if cai_scores else None
        self.tai_scores = tai_scores
        self._tai_scores_values = tai_scores.values() if tai_scores else None
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

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "name": self.name,
            "wanted": self.is_optimized,
            "optimization_priority": self.optimization_priority,
            "cai_weights": self.cai_profile,
            "tai_weights": self.tai_profile,
        }


class OptimizationMethod(Enum):
    single_codon_ratio = "single_codon_ratio"
    single_codon_diff = "single_codon_diff"
    zscore_single_aa_average = "zscore_single_aa_average"
    zscore_bulk_aa_average = "zscore_bulk_aa_average"
    zscore_single_aa_weakest_link = "zscore_single_aa_weakest_link"
    zscore_bulk_aa_weakest_link = "zscore_bulk_aa_weakest_link"

    @property
    def is_single_codon_optimization(self) -> bool:
        return self in (OptimizationMethod.single_codon_ratio, OptimizationMethod.single_codon_diff)

    @property
    def is_zscore_single_aa_optimization(self) -> bool:
        return self in (OptimizationMethod.zscore_single_aa_average, OptimizationMethod.zscore_single_aa_weakest_link)

    @property
    def is_zscore_bulk_aa_optimization(self) -> bool:
        return self in (OptimizationMethod.zscore_bulk_aa_average, OptimizationMethod.zscore_bulk_aa_weakest_link)

    @property
    def is_zscore_average_score_optimization(self) -> bool:
        return self in (OptimizationMethod.zscore_single_aa_average, OptimizationMethod.zscore_bulk_aa_average)

    @property
    def is_zscore_weakest_link_score_optimization(self) -> bool:
        return self in (OptimizationMethod.zscore_single_aa_weakest_link,
                        OptimizationMethod.zscore_bulk_aa_weakest_link)


class OptimizationCubIndex(Enum):
    codon_adaptation_index = "CAI"
    trna_adaptation_index = "tAI"
    max_codon_trna_adaptation_index = "max_CAI_tAI"

    @property
    def is_codon_adaptation_score(self) -> bool:
        return self in (OptimizationCubIndex.codon_adaptation_index,
                        OptimizationCubIndex.max_codon_trna_adaptation_index)

    @property
    def is_trna_adaptation_score(self) -> bool:
        return self in (OptimizationCubIndex.trna_adaptation_index,
                        OptimizationCubIndex.max_codon_trna_adaptation_index)


@dataclass
class UserInput:
    organisms: typing.List[Organism]
    sequence: str
    zip_directory: str
    tuning_parameter: float
    clusters_count: int
    optimization_method: OptimizationMethod = OptimizationMethod.zscore_bulk_aa_average
    optimization_cub_score: OptimizationCubIndex = OptimizationCubIndex.max_codon_trna_adaptation_index

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "sequence": self.sequence,
            "tuning_parameter": self.tuning_parameter,
            "optimization_method": self.optimization_method.value,
            "optimization_cub_score": self.optimization_cub_score.value,
            "organisms": [organism.summary for organism in self.organisms],
        }
