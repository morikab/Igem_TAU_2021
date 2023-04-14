import json
import os


def generate_summary():
    base_directory = "results"
    filename = "run_summary.json"

    for root, dirs, files in os.walk(base_directory):
        for file in files:
            if file == filename:
                filepath = os.path.join(root, file)
                parse_summary_file(filepath)


def parse_summary_file(file_path: str) -> None:
    with open(file_path, "r") as summary_file:
        summary = json.load(summary_file)
        print(summary["evaluation"])

        # TODO - continue from here: * create a row for each amino-acid with the chosen codon


if __name__ == "__main__":
    generate_summary()
