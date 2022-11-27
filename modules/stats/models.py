from dataclasses import dataclass


@dataclass
class EvaluationModuleResult:
    sequence: str
    mean_opt_index: float
    mean_deopt_index: float
    optimization_index: float
    weakest_score: float
