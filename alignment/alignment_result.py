import math
import threading
from abc import ABC
from collections import defaultdict, Counter
from dataclasses import dataclass
from itertools import chain
from pathlib import Path

import numpy
from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.PloidPedigree import PloidPedigree

from alignment.configuration import section_separator, subsection_separator, save_edge_alignments, \
    PosteriorProbabilitiesCalculationMode, save_edge_to_path_map
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.utility.alignment_utility import convert_ploid_id_to_individual
from scripts.utility.basic_utility import round_down, float_not_greater, verify_and_cap_probability


def save_text_to_file(text: str, path: str):
    with open(path, "a") as f:
        f.write(text)


def get_vertex_alignment_estimated_length(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                                          alignment: dict):
    def get_coalescent_vertex_mapped_level_difference(ancestor: int, descendant: int) -> int:
        return pedigree.get_minimal_path_length(ancestor=alignment[ancestor],
                                                descendant=alignment[descendant]
                                                )

    total_length = 0
    for parent, children in coalescent_tree.succ.items():
        if parent not in alignment:
            # The vertex does not belong to this clade
            continue
        for child in children:
            level_distance_approximation = get_coalescent_vertex_mapped_level_difference(ancestor=parent,
                                                                                         descendant=child)
            total_length += level_distance_approximation
    return total_length


def get_vertex_alignment_phasing_accuracy(vertex_alignment: dict[int, int], proband_truth_alignment: dict[int, int]):
    proband_number = len(proband_truth_alignment)
    cumulative_phasing_accuracy = 0
    for proband, proband_pedigree_candidate in proband_truth_alignment.items():
        proband_alignment_candidate = vertex_alignment[proband]
        if proband_pedigree_candidate == proband_alignment_candidate:
            cumulative_phasing_accuracy += 1
    resulting_accuracy = cumulative_phasing_accuracy / proband_number
    return resulting_accuracy


@dataclass
class AlignmentResult(ABC):
    clade_root: int | None


@dataclass
class FailedClimbingAlignmentResult(AlignmentResult):
    failed_vertex: int
    child_to_pedigree_candidates: dict[int, list[int]]


@dataclass
class VertexAlignmentPosteriorProbabilities:
    # The posterior probability for a vertex given a fixed vertex alignment. This value can be used
    # to calculate the posterior inclusion probability for a vertex
    # Specifically, for we can take the sum vertex_alignment_posterior_probability * vertex_posterior_probability
    vertex_posterior_probabilities_for_vertex_alignment: dict[int, float] | None
    # The sum of simple edge alignment probabilities for the given vertex alignment. We use this value later
    # to calculate the posterior probability of a vertex alignment
    vertex_alignment_probability: float


@dataclass
class FullAlignmentResult(AlignmentResult):
    vertex_alignment: dict
    is_valid: bool
    example_edge_alignment: dict | None = None
    edge_alignments: list[dict] | None = None
    posterior_probabilities: VertexAlignmentPosteriorProbabilities | None = None
    result_filepath: str | Path | None = None
    edge_to_pedigree_paths_map: dict | None = None

    def clean(self):
        # Removes all the heavy data that must not be used by the client code
        self.edge_alignments = None
        self.example_edge_alignment = None
        self.edge_to_pedigree_paths_map = None

    # Syntax sugar allowing us to write if alignment_result: ...
    # This can be useful when we only want to check whether the result is correct
    def __bool__(self):
        return self.is_valid

    # Syntax sugar allowing us to avoid accessing the vertex_alignment member
    def __getitem__(self, key):
        return self.vertex_alignment[key]

    # Syntax sugar for printing the alignment
    def __str__(self):
        return str(self.vertex_alignment)

    def save_to_file(self, filepath: str | Path, tree: CoalescentTree):
        self.result_filepath = filepath
        text = self._prepare_text_for_writing(tree=tree)
        # Launch the IO operation on a separate thread
        threading.Thread(target=save_text_to_file, args=(text, filepath), daemon=True).start()

    @staticmethod
    def format_path(path_to_format):
        path_to_format = [convert_ploid_id_to_individual(x) for x in path_to_format]
        return ", ".join(path_to_format)

    @staticmethod
    def _save_edge_alignments(text_lines: list[str], alignment, edge_alignments, tree: CoalescentTree):
        if not save_edge_alignments:
            return
        count = len(edge_alignments)
        text_lines.append(subsection_separator)
        verb = "is" if count == 1 else "are"
        noun = "alignment" if count == 1 else "alignments"
        text_lines.append(f"There {verb} {count} edge {noun}\n")
        tree_levels = tree.get_levels()[1:]
        padding = 4 * ' '
        text_lines.append(subsection_separator)
        # Print the edge alignments
        for edge_alignment in edge_alignments:
            for level in tree_levels:
                for vertex in level:
                    # Skip the vertices from other clades
                    if vertex not in alignment:
                        continue
                    vertex_children = tree.get_children(vertex)
                    text_lines.append(f"{vertex}:\n")
                    for vertex_child in vertex_children:
                        edge_mapping = edge_alignment[(vertex_child, vertex)]
                        formatted_path = FullAlignmentResult.format_path(edge_mapping)
                        text_lines.append(f"{padding}({vertex_child}, {vertex}): [{formatted_path}]\n")
            text_lines.append(subsection_separator)

    def calculate_alignment_probabilities(self, calculation_mode: PosteriorProbabilitiesCalculationMode):
        if calculation_mode == PosteriorProbabilitiesCalculationMode.SKIP:
            return
        assert self.edge_alignments
        all_lengths = [ea.length for ea in self.edge_alignments]
        min_length = min(all_lengths)
        scaled_sum = sum(2.0 ** (min_length - length) for length in all_lengths)
        log_scaled_sum = math.log2(scaled_sum) - min_length
        sum_edge_alignment_probabilities = 2 ** log_scaled_sum
        vertex_to_probability = None
        if calculation_mode == PosteriorProbabilitiesCalculationMode.INDIVIDUAL_INCLUSION_PROBABILITY:
            vertex_to_probability = dict()
            vertex_to_alignment_lengths = defaultdict(list)
            for index, edge_alignment in enumerate(self.edge_alignments):
                edge_alignment_length = all_lengths[index]
                for pedigree_path in edge_alignment.values():
                    for pedigree_vertex in pedigree_path[1:-1]:
                        vertex_to_alignment_lengths[pedigree_vertex].append(edge_alignment_length)
            for vertex, alignment_lengths in vertex_to_alignment_lengths.items():
                min_v_length = min(alignment_lengths)
                scaled_vertex_sum = sum(2.0 ** (min_v_length - L) for L in alignment_lengths)
                log_v_sum = math.log2(scaled_vertex_sum) - min_v_length
                vertex_to_probability[vertex] = 2 ** (log_v_sum - log_scaled_sum)
            # The pedigree candidates for coalescent vertices appear in all edge alignments
            for pedigree_vertex in self.vertex_alignment.values():
                vertex_to_probability[pedigree_vertex] = 1.0
        self.posterior_probabilities = VertexAlignmentPosteriorProbabilities(
            vertex_posterior_probabilities_for_vertex_alignment=vertex_to_probability,
            vertex_alignment_probability=sum_edge_alignment_probabilities
        )

    def save_vertex_inclusion_probabilities(self, global_vertex_inclusion_probability: dict[int, float]):
        assert self.result_filepath
        with (open(self.result_filepath, "a") as file):
            file.write(subsection_separator)
            file.write("Pedigree vertex appearance probabilities:\n")
            file.write("The format is:\n")
            file.write("{ploid_id}: {vertex_inclusion_probability_given_current_vertex_alignment}"
                       " ({vertex_inclusion_probability})\n")

            for vertex, vertex_posterior_inclusion_probability in sorted(
                    self.posterior_probabilities.vertex_posterior_probabilities_for_vertex_alignment.items(),
                    key=lambda x: x[1],
                    reverse=True):
                assert vertex_posterior_inclusion_probability <= 1.0
                global_vertex_probability = global_vertex_inclusion_probability[vertex]
                converted_vertex = convert_ploid_id_to_individual(vertex)
                rounded_prob = round_down(vertex_posterior_inclusion_probability, 5)
                global_probability_rounded_prob = round_down(global_vertex_probability, 5)
                file.write(f"{converted_vertex}: {rounded_prob} ({global_probability_rounded_prob})\n")

    @staticmethod
    def _build_edge_to_path_map(edge_alignments) -> dict:
        edge_to_path_map = defaultdict(set)
        for edge_alignment in edge_alignments:
            for edge, path in edge_alignment.items():
                edge_to_path_map[edge].add(tuple(path))
        return edge_to_path_map

    def _save_edge_to_path_map(self, text_lines: list[str], edge_alignments):
        def print_edge_to_path_map():
            path_padding = 4 * ' '
            text_lines.append("These are all the possible paths for every edge:\n")
            for edge_to_print, paths_to_print in edge_to_path_map.items():
                paths_count = len(paths_to_print)
                noun = "option" if paths_count == 1 else "options"
                text_lines.append(f"{edge_to_print} ({paths_count} {noun}):\n")
                for path_to_print in paths_to_print:
                    formatted_path = self.format_path(path_to_print)
                    text_lines.append(f"{path_padding}[{formatted_path}]\n")

        if not save_edge_to_path_map:
            return
        edge_to_path_map = self._build_edge_to_path_map(edge_alignments=edge_alignments)
        count = len(edge_alignments)
        # Identify whether the alignments are simply the Cartesian product of all the assignments for the edges
        text_lines.append(subsection_separator)
        cartesian_product_size = math.prod(len(v) for v in edge_to_path_map.values())
        if count == cartesian_product_size:
            text_lines.append(f"The set of all the {count} edge alignments is just the Cartesian product.\n"
                              "In other words, you can pick any possible path for every edge\n")
        else:
            text_lines.append(f"There are {cartesian_product_size - count} rejected edge alignments our of all "
                              f"{cartesian_product_size} theoretically possible (the Cartesian product)\n")
        print_edge_to_path_map()
        text_lines.append(subsection_separator)
        text_lines.append(f"Below are the unavoidable vertices for every edge:\n")
        # Identify the unavoidable vertices for every edge
        for edge, paths in edge_to_path_map.items():
            if not paths:
                unavoidable_vertices = set()
            else:
                iterator = iter(paths)
                unavoidable_vertices = set(next(iterator))
                for path in iterator:
                    unavoidable_vertices.intersection_update(path)
            formatted_vertices = self.format_path(unavoidable_vertices)
            text_lines.append(f"{edge}: [{formatted_vertices}]\n")

    def _prepare_text_for_writing(self, tree: CoalescentTree):
        text_lines = []
        alignment = self.vertex_alignment
        for key, value in alignment.items():
            converted_value = convert_ploid_id_to_individual(value)
            text_lines.append(f"{key}: {converted_value}\n")
        if self.example_edge_alignment is not None:
            text_lines.append(section_separator)
            text_lines.append(f"Example edge alignment:\n")
            for edge, path in self.example_edge_alignment.items():
                text_lines.append(f"{edge}: {path}\n")
        elif self.edge_alignments is not None:
            text_lines.append(section_separator)
            edge_alignments = self.edge_alignments
            if not edge_alignments:
                text_lines.append("There are no edges in the clade. Therefore, there are 0 edge alignments\n")
                full_text = "".join(text_lines)
                return full_text
            # Save the edge to path map
            self._save_edge_to_path_map(text_lines=text_lines,
                                        edge_alignments=edge_alignments)
            # Save edge alignments
            self._save_edge_alignments(text_lines=text_lines, alignment=alignment,
                                       edge_alignments=edge_alignments, tree=tree)
        full_text = "".join(text_lines)
        return full_text


def verify_probabilities_values(inclusion_probability_dictionary: dict[int, float]):
    assert all(v <= 1 for v in inclusion_probability_dictionary.values()), \
        f"Values >1 found: {[f'{k}: {v}' for k, v in inclusion_probability_dictionary.items() if v > 1]}"


def verify_pedigree_candidates_inclusion_probabilities(vertex_alignment: dict[int, int],
                                                       inclusion_probability_dictionary: dict[int, float]):
    for pedigree_candidate in vertex_alignment.values():
        inclusion_probability = inclusion_probability_dictionary[pedigree_candidate]
        assert inclusion_probability == 1.0, (f"The inclusion probability "
                                              f"of a pedigree candidate {pedigree_candidate}"
                                              f"must be 1, got {inclusion_probability} instead")


def verify_vertex_inclusion_likelihoods_are_consistent_for_vertex_alignment(
        inclusion_probability_dictionary: dict[int, float],
        vertex_alignment: dict[int, int],
        tree: CoalescentTree,
        clade_root: int,
        pedigree: PloidPedigree,
        initial_mapping: dict[int, [int]]):
    verify_probabilities_values(inclusion_probability_dictionary)
    verify_pedigree_candidates_inclusion_probabilities(vertex_alignment, inclusion_probability_dictionary)
    reverse_alignment = {v: k for k, v in vertex_alignment.items()}
    carriers = set(chain.from_iterable(initial_mapping.values()))
    pedigree_candidates = set(vertex_alignment.values())
    clade_root_candidate = vertex_alignment[clade_root]
    # Verify that the posterior probabilities are correct
    for parent_ploid in pedigree:
        other_ploid = pedigree.get_other_ploid(parent_ploid)
        parents = {parent_ploid, other_ploid}
        if parent_ploid in carriers or other_ploid in carriers:
            continue
        children = pedigree.get_children(parent_ploid)
        # Non-carriers with no children or parents of the clade root's pedigree candidate
        # can never appear in an edge alignment
        if not children or clade_root_candidate in children:
            non_carriers_with_no_children_inclusion_probability = sum(
                inclusion_probability_dictionary.get(x, 0) for x in parents)
            assert non_carriers_with_no_children_inclusion_probability == 0.0
            continue
        expected_parents = set(chain.from_iterable([pedigree.get_parents(x) for x in children]))
        assert set(pedigree.get_children(other_ploid)) == set(children)
        assert expected_parents == parents, f"Expected {expected_parents} as parents, got {parents} instead"
        children_inclusion_probability = sum(inclusion_probability_dictionary.get(x, 0) for x in children)
        parent_candidates = pedigree_candidates.intersection(parents)
        other_parents = [x for x in parents if x not in parent_candidates]
        parents_inclusion_probability = 0
        for parent_candidate in parent_candidates:
            candidate_probability = inclusion_probability_dictionary.get(parent_candidate, 0)
            tree_vertex = reverse_alignment[parent_candidate]
            children_number = len(tree.get_children(tree_vertex))
            parents_inclusion_probability += children_number * candidate_probability
        parents_inclusion_probability += sum(inclusion_probability_dictionary.get(x, 0) for x in other_parents)
        assert numpy.isclose(parents_inclusion_probability, children_inclusion_probability)


def verify_initial_mapping_vertex_probabilities(inclusion_probability_dictionary: dict[int, float],
                                                pedigree: PloidPedigree,
                                                initial_mapping: dict[int, [int]]):
    all_pedigree_candidates = [v for lst in initial_mapping.values() for v in lst]
    # Count occurrences
    value_counts = Counter(all_pedigree_candidates)
    # Keep only values that appear exactly once
    unique_pedigree_candidates = {v for v, count in value_counts.items() if count == 1}
    for candidate_vertices in initial_mapping.values():
        probabilities_sum = sum(inclusion_probability_dictionary.get(x, 0) for x in candidate_vertices)
        if unique_pedigree_candidates.issuperset(candidate_vertices):
            assert numpy.isclose(probabilities_sum, 1.0), \
                (f"The sum of posterior probabilities of unique pedigree candidates {candidate_vertices} "
                 f"must be 1, got {probabilities_sum} instead")
        else:
            assert float_not_greater(1, probabilities_sum), \
                (f"The total inclusion probabilities sum for the "
                 f"coalescent vertex candidates must be at least 1, "
                 f"got {probabilities_sum} instead")
        candidate_parents = chain.from_iterable(pedigree.get_parents(x) for x in candidate_vertices)
        parents_probabilities_sum = sum(inclusion_probability_dictionary.get(x, 0) for x in candidate_parents)
        assert float_not_greater(1, parents_probabilities_sum), \
            ("The initial_mapping's parents' probabilities "
             "must be at least 1, got "
             f"{parents_probabilities_sum} instead")


def verify_global_vertex_inclusion_likelihoods_are_consistent(inclusion_probability_dictionary: dict[int, float],
                                                              pedigree: PloidPedigree,
                                                              tree_root: int,
                                                              initial_mapping: dict[int, [int]],
                                                              vertex_alignments: [FullAlignmentResult]
                                                              ):
    verify_probabilities_values(inclusion_probability_dictionary)
    verify_initial_mapping_vertex_probabilities(initial_mapping=initial_mapping,
                                                pedigree=pedigree,
                                                inclusion_probability_dictionary=inclusion_probability_dictionary)
    tree_candidates = {alignment_result.vertex_alignment[tree_root] for alignment_result in vertex_alignments}
    for child_ploid in pedigree:
        if child_ploid in tree_candidates:
            continue
        parents = pedigree.get_parents(child_ploid)
        if not parents:
            continue
        child_inclusion_probability = inclusion_probability_dictionary.get(child_ploid, 0)
        parents_inclusion_probability = sum(inclusion_probability_dictionary.get(x, 0) for x in parents)
        assert float_not_greater(child_inclusion_probability, parents_inclusion_probability)


class CladeAlignmentResults(ABC):
    pass


@dataclass
class FailedClimbingCladeAlignmentResults(CladeAlignmentResults):
    failed_climbing_alignment_info: FailedClimbingAlignmentResult


@dataclass
class CladeAlignmentPosteriorProbabilities:
    # Dictionary mapping the vertex alignment index to the posterior probability
    vertex_alignment_to_posterior_probability: dict[int, float]
    vertex_posterior_probabilities: dict[int, float] | None


class SuccessCladeAlignmentResults(CladeAlignmentResults):
    clade_root: int
    alignments: list[FullAlignmentResult]
    clade_alignment_posterior_probabilities: CladeAlignmentPosteriorProbabilities | None = None

    def __init__(self, clade_root: int, alignments=None):
        self.clade_root = clade_root
        if alignments is None:
            alignments = []
        self.alignments = alignments

    def calculate_posterior_probabilities(self):
        if self.clade_alignment_posterior_probabilities:
            return
        calculate_vertex_inclusion_probabilities = True
        resulting_global_vertex_inclusion_probability_dict = None
        # Verify that all the vertex alignments have posterior probabilities calculated
        for alignment in self.alignments:
            if not alignment.posterior_probabilities:
                return
            if not alignment.posterior_probabilities.vertex_posterior_probabilities_for_vertex_alignment:
                calculate_vertex_inclusion_probabilities = False
        vertex_alignment_probabilities_sum = sum(x.posterior_probabilities.vertex_alignment_probability
                                                 for x in self.alignments)
        vertex_alignment_posterior_probability = dict()
        for index, vertex_alignment in enumerate(self.alignments):
            alignment_posterior_probability = (vertex_alignment.posterior_probabilities.vertex_alignment_probability /
                                               vertex_alignment_probabilities_sum)
            alignment_posterior_probability = verify_and_cap_probability(alignment_posterior_probability)
            vertex_alignment_posterior_probability[index] = alignment_posterior_probability
        if not calculate_vertex_inclusion_probabilities:
            self.clade_alignment_posterior_probabilities = CladeAlignmentPosteriorProbabilities(
                vertex_alignment_to_posterior_probability=vertex_alignment_posterior_probability,
                vertex_posterior_probabilities=resulting_global_vertex_inclusion_probability_dict
            )
            return
        global_vertex_inclusion_probability_dict = defaultdict(float)
        for vertex_alignment_index, posterior_probability in vertex_alignment_posterior_probability.items():
            vertex_alignment = self.alignments[vertex_alignment_index]
            for vertex, vertex_inclusion_probability in vertex_alignment.posterior_probabilities.vertex_posterior_probabilities_for_vertex_alignment.items():
                global_vertex_inclusion_probability_dict[vertex] += vertex_inclusion_probability * posterior_probability
        resulting_global_vertex_inclusion_probability_dict = {
            vertex: verify_and_cap_probability(probability)
            for vertex, probability in global_vertex_inclusion_probability_dict.items()
        }
        self.clade_alignment_posterior_probabilities = CladeAlignmentPosteriorProbabilities(
            vertex_alignment_to_posterior_probability=vertex_alignment_posterior_probability,
            vertex_posterior_probabilities=resulting_global_vertex_inclusion_probability_dict
        )


@dataclass
class TreeAlignmentResults:
    clade_root_to_clade_results: dict[int, CladeAlignmentResults]

    def get_unique_clade_results(self) -> CladeAlignmentResults:
        if len(self.clade_root_to_clade_results) != 1:
            raise ValueError("There must be exactly one clade in the tree")
        return self.clade_root_to_clade_results[next(iter(self.clade_root_to_clade_results))]
