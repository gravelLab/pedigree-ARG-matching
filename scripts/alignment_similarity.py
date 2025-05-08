from collections import defaultdict
from itertools import combinations

from matplotlib import pyplot as plt


def calculate_percentage_of_correct_assignments(alignment: dict):
    correct_assignments = sum(1 for k, v in alignment.items() if k == v)
    return correct_assignments / len(alignment)


def get_alignments_key_similarity_probands_ignored(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    # Calculate the number of shared values. Assuming here that the keys in both the dictionaries are the same
    return sum(1 for key in first if key not in probands and first[key] == second[key])


def get_alignments_pair_similarity(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    probands = set(probands).intersection(first)
    # Calculate the number of shared values
    shared_values = get_alignments_key_similarity_probands_ignored(first, second, probands)
    total_size = len(first) - len(probands)
    similarity = shared_values / total_size
    return similarity


def get_alignments_similarity(alignments: [dict], probands):
    total_similarity = 0
    for dict1, dict2 in combinations(alignments, 2):
        total_similarity += get_alignments_pair_similarity(first=dict1, second=dict2, probands=probands)
    return total_similarity / (len(alignments) * (len(alignments) - 1) / 2)


def get_alignments_pair_ind_similarity(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    probands = set(probands).intersection(first)
    # Calculate the number of shared values
    shared_values = sum(1 for key in first if key not in probands and first[key] // 2 == second[key] // 2)
    total_size = len(first) - len(probands)
    similarity = shared_values / total_size
    return similarity


def get_alignments_ind_similarity(alignments: [dict], probands):
    total_similarity = 0
    for dict1, dict2 in combinations(alignments, 2):
        total_similarity += get_alignments_pair_ind_similarity(first=dict1, second=dict2, probands=probands)
    return total_similarity / (len(alignments) * (len(alignments) - 1) / 2)


def get_alignments_proband_distance_probands_ignored(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    # Calculate the number of shared values. Assuming here that the keys in both the dictionaries are the same
    distance = sum(1 for key in first if key not in probands and first[key] != second[key])
    assert distance == len(first) - len(probands) - sum(1 for key in first if key not in probands and
                                                        first[key] == second[key])
    return distance


def get_distance_histogram_to_identity(alignments: [dict]):
    # Assuming that all the alignments have the same keys
    coalescent_vertices = alignments[0].keys()
    vertices_length = len(coalescent_vertices)
    identity_solution = {x: x for x in coalescent_vertices}
    distance_histogram = defaultdict(int)
    for alignment in alignments:
        alignment_distance = (vertices_length -
                              sum(1 for key in identity_solution if identity_solution[key] == alignment[key]))
        # alignment_distance = get_alignments_pair_distance_probands_ignored(alignment, identity_solution, probands)
        distance_histogram[alignment_distance] += 1
    # Verifying that the identity alignment is always present
    return distance_histogram


def save_distance_histogram(histogram_filepath: str, dictionary: dict):
    x = list(dictionary.keys())
    y = list(dictionary.values())
    plt.bar(x, y)
    plt.xlabel('Distances')
    plt.ylabel('Frequency')
    plt.title('Distance histogram')
    plt.savefig(f"{histogram_filepath}.png")
    plt.close()
