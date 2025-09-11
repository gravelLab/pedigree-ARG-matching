import math
import threading
from abc import ABC
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import List

from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.configuration import section_separator, subsection_separator
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.utility.alignment_utility import convert_ploid_id_to_individual


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


def get_edge_alignment_probability(edge_alignment: dict) -> float:
    total_length = 0
    for pedigree_path in edge_alignment.values():
        total_length += len(pedigree_path)
    return Decimal(2) ** (-total_length)


@dataclass
class AlignmentResult(ABC):
    clade_root: int | None


@dataclass
class FailedClimbingAlignmentResult(AlignmentResult):
    failed_vertex: int
    child_to_pedigree_candidates: dict[int, list[int]]


@dataclass
class FullAlignmentResult(AlignmentResult):
    vertex_alignment: dict
    is_valid: bool
    example_edge_alignment: dict | None = None
    edge_alignments: list[dict] | None = None

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
        text = self._prepare_text_for_writing(tree=tree)
        # Launch the IO operation on a separate thread
        threading.Thread(target=save_text_to_file, args=(text, filepath), daemon=True).start()

    @staticmethod
    def format_path(path_to_format):
        path_to_format = [convert_ploid_id_to_individual(x) for x in path_to_format]
        return ", ".join(path_to_format)

    @staticmethod
    def _save_edge_alignments(text_lines: list[str], alignment, edge_alignments, tree: CoalescentTree):
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

    @staticmethod
    def _save_vertex_probability(text_lines: list[str], vertex_alignment, edge_alignments, tree: CoalescentTree):
        assert edge_alignments
        edge_alignment_probabilities = []
        vertex_to_probability = defaultdict(Decimal)
        for edge_alignment in edge_alignments:
            edge_alignment_probability = get_edge_alignment_probability(edge_alignment)
            assert edge_alignment_probability > 0
            edge_alignment_probabilities.append(edge_alignment_probability)
            for pedigree_path in edge_alignment.values():
                for pedigree_vertex in pedigree_path[1:-1]:
                    vertex_to_probability[pedigree_vertex] += edge_alignment_probability
            for coalescent_vertex in tree:
                pedigree_vertex = vertex_alignment[coalescent_vertex]
                vertex_to_probability[pedigree_vertex] += edge_alignment_probability
        sum_edge_alignment_probabilities = sum(edge_alignment_probabilities)
        assert sum_edge_alignment_probabilities != 0
        vertex_to_probability = {key: value / sum_edge_alignment_probabilities
                                 for key, value in vertex_to_probability.items()}
        text_lines.append(subsection_separator)
        text_lines.append("Pedigree vertex appearance probabilities:\n")

        def round_down(x, decimals=5):
            factor = 10 ** decimals
            return math.floor(x * factor) / factor

        for vertex, edge_alignment_probability in sorted(vertex_to_probability.items(), key=lambda x: x[1], reverse=True):
            assert edge_alignment_probability <= 1.0
            converted_vertex = convert_ploid_id_to_individual(vertex)
            rounded_prob = round_down(edge_alignment_probability, 5)
            text_lines.append(f"{converted_vertex}: {rounded_prob}\n")

    @staticmethod
    def _build_edge_to_path_map(edge_alignments) -> dict:
        edge_to_path_map = defaultdict(set)
        for edge_alignment in edge_alignments:
            for edge, path in edge_alignment.items():
                edge_to_path_map[edge].add(tuple(path))
        return dict(edge_to_path_map)

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
            edge_alignments = [x for x in self.edge_alignments if x]
            count = len(edge_alignments)
            if not count:
                text_lines.append("There are no edges in the clade. Therefore, there are 0 edge alignments\n")
                full_text = "".join(text_lines)
                return full_text
            # Save the edge to path map
            self._save_edge_to_path_map(text_lines=text_lines,
                                        edge_alignments=edge_alignments)
            # Calculate the probability for every pedigree vertex
            self._save_vertex_probability(text_lines=text_lines, vertex_alignment=alignment,
                                          edge_alignments=edge_alignments, tree=tree)
            # Save the edge alignments
            self._save_edge_alignments(text_lines=text_lines, alignment=alignment,
                                       edge_alignments=edge_alignments, tree=tree)
        full_text = "".join(text_lines)
        return full_text


class CladeAlignmentResults(ABC):
    pass


@dataclass
class FailedClimbingCladeAlignmentResults(CladeAlignmentResults):
    failed_climbing_alignment_info: FailedClimbingAlignmentResult


class SuccessCladeAlignmentResults(CladeAlignmentResults):
    clade_root: int
    alignments: List[FullAlignmentResult]
    pedigree_vertex_to_edge_alignment_appearance_number: defaultdict[int, int]
    edge_alignment_total_number: int = 0

    def __init__(self, clade_root: int, alignments=None,
                 pedigree_vertex_to_edge_alignment_appearance_number=None,
                 edge_alignment_total_number=None):
        self.clade_root = clade_root
        if alignments is None:
            alignments = []
        self.alignments = alignments
        if pedigree_vertex_to_edge_alignment_appearance_number is None:
            pedigree_vertex_to_edge_alignment_appearance_number = defaultdict(int)
        self.pedigree_vertex_to_edge_alignment_appearance_number = pedigree_vertex_to_edge_alignment_appearance_number
        if edge_alignment_total_number is None:
            edge_alignment_total_number = 0
        self.edge_alignment_total_number = edge_alignment_total_number


@dataclass
class TreeAlignmentResults:
    clade_root_to_clade_results: dict[int, CladeAlignmentResults]

    def get_unique_clade_results(self) -> CladeAlignmentResults:
        if len(self.clade_root_to_clade_results) != 1:
            raise ValueError("There must be exactly one clade in the tree")
        return self.clade_root_to_clade_results[next(iter(self.clade_root_to_clade_results))]
