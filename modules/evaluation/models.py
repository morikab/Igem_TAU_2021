import typing
from dataclasses import dataclass


@dataclass
class EvaluationModuleResult:
    # TODO - Add additional stats (% optimized, % deoptimized, etc.)
    sequence: str
    average_distance_score: float
    weakest_link_score: float
    ratio_score: float

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "final_sequence": self.sequence,
            "average_distance_score": self.average_distance_score,
            "weakest_link_score": self.weakest_link_score,
            "ratio_score": self.ratio_score,
        }
