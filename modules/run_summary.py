import json
import os
import typing


class RunSummary(object):
    _run_summary = {}

    @classmethod
    def add_to_run_summary(cls, key: str, value: typing.Any) -> None:
        if key in cls._run_summary:
            raise KeyError(F"Key {key} already exists in run summary")
        cls._run_summary[key] = value

    @classmethod
    def put_in_run_summary(cls, key: str, value: typing.Any) -> None:
        cls._run_summary[key] = value

    @classmethod
    def delete_from_run_summary(cls, key:str) -> None:
        cls._run_summary.pop(key)

    @classmethod
    def save_run_summary(cls, output_directory: str) -> None:
        output_path = os.path.join(output_directory, "run_summary.json")
        with open(output_path, "w") as output_file:
            json.dump(cls._run_summary, output_file)

    @classmethod
    def reset(cls) -> None:
        cls._run_summary = {}

    @classmethod
    def get(cls) -> typing.Dict[str, typing.Any]:
        return cls._run_summary
