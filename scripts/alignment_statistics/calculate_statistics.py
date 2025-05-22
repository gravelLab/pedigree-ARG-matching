import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from alignment.configuration import calculate_distances_histogram
from alignment.driver_file import PloidType
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.coalescent_tree import CoalescentTree
from scripts.alignment_similarity import get_alignments_similarity, get_alignments_ind_similarity

section_separator = "############################################################\n"


def convert_ploid_id_to_individual(ploid_id: int):
    individual_id = ploid_id // 2
    ploid_type = PloidType.Paternal.value if ploid_id % 2 == 0 else PloidType.Maternal.value
    return f"{individual_id}{ploid_type}"


def save_alignment_to_file(dictionary_file, alignment: dict):
    for key, value in alignment.items():
        converted_value = convert_ploid_id_to_individual(value)
        dictionary_file.write(f"{key}: {converted_value}\n")


def get_alignment_likelihood(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                             alignment: dict):
    def get_coalescent_vertex_mapped_level_difference(ancestor: int, descendant: int) -> int:
        return pedigree.get_minimal_path_length(ancestor=alignment[ancestor],
                                                descendant=alignment[descendant]
                                                )

    total_length = 0
    for parent, children in coalescent_tree.children_map.items():
        if parent not in alignment:
            # The vertex does not belong to this clade
            continue
        for child in children:
            level_distance_approximation = get_coalescent_vertex_mapped_level_difference(ancestor=parent,
                                                                                         descendant=child)
            total_length += level_distance_approximation
    return total_length


def calculate_coalescent_vertex_pedigree_candidates(clade_alignment_result: [dict]):
    vertex_candidate_count_dict = defaultdict(lambda: defaultdict(int))
    for alignment in clade_alignment_result:
        for key, value in alignment.items():
            vertex_candidate_count_dict[key][value] += 1
    return vertex_candidate_count_dict


def save_coalescent_vertex_pedigree_candidates(statistics_file, clade_alignment_result: [dict],
                                               clade: Iterable[int]):
    vertex_candidate_count_dict = calculate_coalescent_vertex_pedigree_candidates(clade_alignment_result)
    alignments_number = len(clade_alignment_result)
    statistics_file.write("##############################\n")
    statistics_file.write("Appearances of each candidate:\n")
    for vertex in clade:
        statistics_file.write(f"{vertex}:\n")
        for pedigree_candidate, count in sorted(vertex_candidate_count_dict[vertex].items(),
                                                key=lambda x: x[1], reverse=True):
            percentage = count / alignments_number
            formatted_pedigree_vertex = convert_ploid_id_to_individual(pedigree_candidate)
            statistics_file.write(f"        {formatted_pedigree_vertex} ({percentage});\n")


@dataclass
class CladeAlignmentsMetadata:
    calculate_similarity: bool
    calculate_distances_histogram: bool
    calculate_alignments_likelihoods: bool
    clade_alignments: [dict[int: int]]


@dataclass
class CladeMetadata:
    results_filepath: Path
    clade_root: int
    clade: Iterable[int]
    probands: Iterable[int]
    number_of_coalescing_events: int
    # TODO: Rename the class
    pedigree: PotentialMrcaProcessedGraph
    coalescent_tree: CoalescentTree
    clade_alignments_metadata: CladeAlignmentsMetadata = None

    @staticmethod
    def get_clade_basic_metadata(coalescent_tree: CoalescentTree, clade_root: int,
                                 pedigree: PotentialMrcaProcessedGraph, results_filepath: str | Path,
                                 clade_alignments_metadata: CladeAlignmentsMetadata = None):
        # Calculating the clade for the given vertex and sorting the vertices by their levels
        clade = coalescent_tree.get_vertex_descendants(clade_root)
        # clade = coalescent_tree.get_connected_component_for_vertex(clade_root)
        clade = sorted(clade, key=lambda v: coalescent_tree.vertex_to_level_map[v], reverse=True)
        # Calculating the rest of the data
        coalescing_events = len([x for x in clade if x in coalescent_tree.children_map
                                 and len(coalescent_tree.children_map[x]) > 1])
        proband_vertices = [x for x in clade if coalescent_tree.vertex_to_level_map[x] == 0]
        os.makedirs(results_filepath, exist_ok=True)
        results_filepath = Path(results_filepath)
        return CladeMetadata(
            clade_root=clade_root,
            clade=clade,
            probands=proband_vertices,
            number_of_coalescing_events=coalescing_events,
            clade_alignments_metadata=clade_alignments_metadata,
            coalescent_tree=coalescent_tree,
            pedigree=pedigree,
            results_filepath=results_filepath
        )

    def save_alignments_similarities(self, statistics_file):
        alignments_similarity = get_alignments_similarity(self.clade_alignments_metadata.clade_alignments,
                                                          self.probands)
        alignments_ind_similarity = get_alignments_ind_similarity(self.clade_alignments_metadata.clade_alignments,
                                                                  self.probands)
        statistics_file.write(f"The alignments similarity: {alignments_similarity}\n")
        statistics_file.write(f"The alignments individual similarity: {alignments_ind_similarity}\n")
        statistics_file.write(section_separator)

    def save_distance_histograms(self, statistics_file):
        min_length = sys.maxsize
        min_length_alignments = []
        # TODO: Cache the distance between two vertices per pedigree
        for counter, valid_alignment in enumerate(self.clade_alignments_metadata.clade_alignments):
            alignment_length = get_alignment_likelihood(pedigree=self.pedigree,
                                                        coalescent_tree=self.coalescent_tree,
                                                        alignment=valid_alignment)
            if calculate_distances_histogram:
                if min_length == alignment_length:
                    min_length_alignments.append(counter)
                elif alignment_length < min_length:
                    min_length = alignment_length
                    min_length_alignments = [counter]
        statistics_file.write(f"The minimum alignment length: {min_length}\n")
        statistics_file.write(f"The corresponding alignments: {min_length_alignments}\n")
        statistics_file.write(section_separator)

    def save_basic_alignments_metadata(self, statistics_file):
        statistics_file.write(section_separator)
        alignments_number = len(self.clade_alignments_metadata.clade_alignments)
        statistics_file.write(f"The total number of alignments is: {alignments_number}\n")

    def save(self):
        self.save_alignments()
        self.save_metadata()

    def save_alignments(self):
        for counter, alignment in enumerate(self.clade_alignments_metadata.clade_alignments):
            alignment_filename = f"alignment_{counter}"
            alignment_path = self.results_filepath / alignment_filename
            with open(alignment_path, "w") as dictionary_file:
                if self.clade_alignments_metadata.calculate_alignments_likelihoods:
                    alignment_likelihood = get_alignment_likelihood(pedigree=self.pedigree,
                                                                    coalescent_tree=self.coalescent_tree,
                                                                    alignment=alignment)
                    dictionary_file.write(section_separator)
                    dictionary_file.write(f"Approximated alignment length: {alignment_likelihood}\n")
                    dictionary_file.write(section_separator)
                save_alignment_to_file(dictionary_file=dictionary_file, alignment=alignment)

    def save_metadata(self):
        metadata_filepath = self.results_filepath / f"_clade_{self.clade_root}.txt"
        with open(metadata_filepath, 'w') as statistics_file:
            self.save_tree_metadata(statistics_file)
            self.save_number_of_pedigree_candidates_per_clade_vertex(statistics_file)
            if not self.clade_alignments_metadata:
                return
            self.save_basic_alignments_metadata(statistics_file)
            if self.clade_alignments_metadata.calculate_similarity:
                self.save_alignments_similarities(statistics_file)
            else:
                statistics_file.write("The similarities among the trees were not calculated\n")
            if self.clade_alignments_metadata.calculate_distances_histogram:
                self.save_distance_histograms(statistics_file)
            else:
                statistics_file.write("The distance histograms were not calculated\n")

    def save_number_of_pedigree_candidates_per_clade_vertex(self, statistics_file):
        coalescent_vertex_pedigree_candidates_number = {x: len({alignment[x] for alignment in
                                                                self.clade_alignments_metadata.clade_alignments})
                                                        for x in self.clade}
        statistics_file.write(section_separator)
        statistics_file.write("The number of pedigree candidates for every vertex:\n")
        for coalescent_vertex, vertex_candidates_number in coalescent_vertex_pedigree_candidates_number.items():
            if coalescent_vertex in self.coalescent_tree.probands:
                statistics_file.write(f"{coalescent_vertex} (proband): {vertex_candidates_number}\n")
            else:
                statistics_file.write(f"{coalescent_vertex}: {vertex_candidates_number}\n")
        save_coalescent_vertex_pedigree_candidates(
            statistics_file=statistics_file,
            clade_alignment_result=self.clade_alignments_metadata.clade_alignments,
            clade=self.clade
        )

    def save_tree_metadata(self, statistics_file):
        statistics_file.write(f"The root of the clade: {self.clade_root}\n")
        statistics_file.write(f"There are {len(self.clade)} vertices in the clade\n")
        statistics_file.write(f"There are {len(self.probands)} probands in the clade\n")
        statistics_file.write(f"There are {self.number_of_coalescing_events} coalescing events in the clade\n")
        statistics_file.write(section_separator)
        statistics_file.write("Number of coalescing events grouped by the children number:\n")
        coalescing_number_to_events_number = dict()
        for vertex in self.clade:
            if vertex not in self.coalescent_tree.children_map:
                continue
            children_number = len(self.coalescent_tree.children_map[vertex])
            previous_counter = coalescing_number_to_events_number.get(children_number, 0)
            coalescing_number_to_events_number[children_number] = previous_counter + 1
        coalescing_numbers = sorted(coalescing_number_to_events_number.keys())
        for coalescing_number in coalescing_numbers:
            statistics_file.write(f"{coalescing_number}: {coalescing_number_to_events_number[coalescing_number]}\n")


def save_statistics_to_file(clade_alignments_metadata: CladeAlignmentsMetadata, coalescent_tree: CoalescentTree,
                            clade_root: int, pedigree: PotentialMrcaProcessedGraph,
                            results_filepath: str):
    # Process the clade
    clade_metadata = CladeMetadata.get_clade_basic_metadata(coalescent_tree=coalescent_tree,
                                                            clade_root=clade_root,
                                                            pedigree=pedigree,
                                                            clade_alignments_metadata=clade_alignments_metadata,
                                                            results_filepath=results_filepath)
    # Printing the results to the file
    clade_metadata.save_tree_metadata(results_filepath)
