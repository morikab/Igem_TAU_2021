import typing
from dataclasses import dataclass


@dataclass
class EvaluationModuleResult:
    sequence: str
    mean_opt_index: float
    mean_deopt_index: float
    optimization_index: float
    weakest_score: float

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "mean_opt_index": self.mean_opt_index,
            "mean_deopt_index": self.mean_deopt_index,
            "optimization_index": self.optimization_index,
            "weakest_score": self.weakest_score,
        }
