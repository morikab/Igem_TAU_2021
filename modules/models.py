import typing
from dataclasses import dataclass
from enum import Enum
from modules.configuration import Configuration


class TranslationFunction(Enum):
    single_codon_global = "single_codon_global"
    zscore_hill_climbing_average = "zscore_hill_climbing_average"
    zscore_hill_climbing_weakest_link = "zscore_hill_climbing_weakest_link"


@dataclass
class ModelPreferences:
    restriction_enzymes: bool
    translation: bool
    transcription: bool
    translation_function: TranslationFunction = TranslationFunction.zscore_hill_climbing_average

    @classmethod
    def init_from_dictionary(cls, model_preferences_dict: typing.Dict[str, typing.Any]):
        return ModelPreferences(
            restriction_enzymes=model_preferences_dict["RE"],
            translation=model_preferences_dict["translation"],
            transcription=model_preferences_dict["transcription"],
            translation_function=TranslationFunction[model_preferences_dict["translation_function"]],
        )

    @classmethod
    def init_from_config(cls):
        config = Configuration.get_config()
        model_preferences = config["MODEL_PREFERENCES"]
        return ModelPreferences(
            restriction_enzymes=model_preferences["RESTRICTION_ENZYMES"],
            translation=model_preferences["TRANSLATION"],
            transcription=model_preferences["TRANSCRIPTION"],
            translation_function=TranslationFunction[model_preferences["TRANSLATION_FUNCTION"]],
        )
