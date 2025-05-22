from alignment.driver_file import PloidType


def dict_has_duplicate_values(dictionary: dict):
    seen = set()
    for value in dictionary.values():
        if value in seen:
            return True
        seen.add(value)
    return False


def dict_is_identity(dictionary: dict) -> bool:
    return all(k == v for k, v in dictionary.items())


def save_dictionary_to_file(dictionary_filepath: str, dictionary: dict):
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
    ploid_type = PloidType.Paternal.value if ploid_id % 2 == 0 else PloidType.Maternal.value
    return f"{individual_id}{ploid_type}"


def save_alignment_to_file(dictionary_file, alignment: dict):
    for key, value in alignment.items():
        converted_value = convert_ploid_id_to_individual(value)
        dictionary_file.write(f"{key}: {converted_value}\n")


def read_mapping_from_file(filepath):
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
