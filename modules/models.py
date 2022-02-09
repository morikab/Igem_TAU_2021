import statistics
import typing
from dataclasses import dataclass
from enum import Enum
from modules.configuration import Configuration


class Organism(object):
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


class TranslationFunction(Enum):
    single_codon_global = "single_codon_global"
    single_codon_local = "single_codon_local"
    zscore_hill_climbing_average = "zscore_hill_climbing_average"
    zscore_hill_climbing_weakest_link = "zscore_hill_climbing_weakest_link"


@dataclass
class ModelPreferences:
    restriction_enzymes: bool
    translation: bool
    translation_function: TranslationFunction = TranslationFunction.zscore_hill_climbing_average

    @classmethod
    def init_from_dictionary(cls, model_preferences_dict: typing.Dict[str, typing.Any]):
        return ModelPreferences(
            restriction_enzymes=model_preferences_dict["RE"],
            translation=model_preferences_dict["translation"],
            translation_function=TranslationFunction[model_preferences_dict["translation_function"]],

        )

    @classmethod
    def init_from_config(cls):
        config = Configuration.get_config()
        model_preferences = config["MODEL_PREFERENCES"]
        return ModelPreferences(
            restriction_enzymes=model_preferences["RESTRICTION_ENZYMES"],
            translation=model_preferences["TRANSLATION"],
            translation_function=TranslationFunction[model_preferences["TRANSLATION_FUNCTION"]],
        )


@dataclass
class UserInput:
    organisms: typing.List[Organism]
    sequence: str
    tuning_parameter: float
