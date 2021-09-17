import typing
from pathlib import Path


def _create_organism_entry_if_possible(organism_key: str,
                                       organism_file_name: str,
                                       key_prefix: str,
                                       is_optimized: bool,
                                       data: typing.Dict) -> typing.Optional[typing.Tuple[str, typing.Dict]]:
    if not organism_key.startswith(F"{key_prefix}_seq_file_#"):
        return None

    organism_index = organism_key.split("#")[1]
    organism_name_key = F"{key_prefix}_name_#{organism_index}"
    organism_name = data[organism_name_key] if data[organism_name_key] != "" else Path(organism_file_name).stem

    organism_entry = {
        "genome_path": organism_file_name,
        "optimized": is_optimized,
        "expression_csv": None,
    }

    return organism_name, organism_entry


def process_input_for_modules(data: typing.Dict) -> typing.Dict:
    user_input = {}
    uploaded_data = data["uploaded_data"]

    user_input["sequence"] = uploaded_data.get("protein_seq", None)
    user_input["selected_promoters"] = None
    if data["promoter_optimization_options"] == "yes":
        pass        # FIXME
    input_organisms = {}
    for organism_key, organism_file_name in uploaded_data.items():
        optimized_organism_entry = _create_organism_entry_if_possible(
            organism_key=organism_key,
            organism_file_name=organism_file_name,
            data=data,
            key_prefix="optimized",
            is_optimized=True,
        )
        if optimized_organism_entry is not None:
            input_organisms[optimized_organism_entry[0]] = optimized_organism_entry[1]

        deoptimized_organism_entry = _create_organism_entry_if_possible(
            organism_key=organism_key,
            organism_file_name=organism_file_name,
            data=data,
            key_prefix="deoptimized",
            is_optimized=False,
        )
        if deoptimized_organism_entry is not None:
            input_organisms[deoptimized_organism_entry[0]] = deoptimized_organism_entry[1]

    user_input["organisms"] = input_organisms

    return user_input
