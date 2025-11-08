from pathlib import Path

import yaml

from alignment.alignment_constants import PloidType, coalescent_id_key, initial_assignments_key, pedigree_ids_key, \
    ploid_types


def dict_has_duplicate_values(dictionary: dict):
    seen = set()
    for value in dictionary.values():
        if value in seen:
            return True
        seen.add(value)
    return False


def dict_is_identity(dictionary: dict) -> bool:
    return all(k == v for k, v in dictionary.items())


def save_dictionary_to_file(dictionary_filepath: str | Path, dictionary: dict):
    dictionary_file = open(dictionary_filepath, 'w')
    for key, value in dictionary.items():
        dictionary_file.write(f"{key}: {value}\n")
    dictionary_file.close()


def parse_dictionary_from_file(file_path: str):
    result = dict()
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            key = int(parts[0].strip())
            value = int(parts[1].strip())
            result[key] = value
    return result


def convert_ploid_id_to_individual(ploid_id: int):
    individual_id = ploid_id // 2
    ploid_type = ploid_types[ploid_id % 2]
    return f"{individual_id}{ploid_type}"


def read_mapping_from_file(filepath: str | Path):
    parsed_dict = {}
    with open(filepath, 'r') as file:
        for line in file:
            # Split the line on ':' to separate key and value
            key, value = line.strip().split(':')
            # Convert key to an integer and value to a list of integers
            key = int(key.strip())
            value = [int(x) for x in value.strip()[1:-1].split(',')]
            # Add to the dictionary
            parsed_dict[key] = value
    return parsed_dict


def convert_ploid_str_to_ploid_id(ploid_str: str) -> int:
    try:
        ploid_type = PloidType(ploid_str[-1])
    except ValueError:
        raise ValueError(
            f"Invalid ploid_type: {ploid_str['ploid_type']}. Must be one of {[e.value for e in PloidType]}"
        )
    try:
        unprocessed_pedigree_id = int(ploid_str[:-1])
    except ValueError:
        raise ValueError(f"Invalid pedigree id {ploid_str[:-1]}")
    if unprocessed_pedigree_id < 0:
        raise ValueError(f"Negative pedigree id {unprocessed_pedigree_id}")
    if ploid_type == PloidType.Paternal:
        return 2 * unprocessed_pedigree_id
    return 2 * unprocessed_pedigree_id + 1


def read_initial_assignments_from_driver_file(driver_filepath: str | Path) -> dict[int, list[int]]:
    with open(driver_filepath, "r") as f:
        data = yaml.safe_load(f)
    assignments = data.get(initial_assignments_key, [])
    result = dict()
    for item in assignments:
        coalescent_id = item[coalescent_id_key]
        pedigree_ids = item[pedigree_ids_key]
        parsed_pedigree_ids = [convert_ploid_str_to_ploid_id(x) for x in pedigree_ids]
        result[coalescent_id] = parsed_pedigree_ids
    return result
