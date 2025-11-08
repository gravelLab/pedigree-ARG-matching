import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.alignment_result import FullAlignmentResult, get_vertex_alignment_estimated_length, \
    FailedClimbingCladeAlignmentResults, SuccessCladeAlignmentResults
from alignment.configuration import calculate_distances_histogram, section_separator
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.alignment_statistics.alignment_similarity import get_alignments_similarity, get_alignments_ind_similarity
from scripts.utility.alignment_utility import convert_ploid_id_to_individual
from scripts.utility.basic_utility import prepend_to_file, round_down


def calculate_coalescent_vertex_pedigree_candidates(clade_alignment_result: [FullAlignmentResult]):
    vertex_candidate_count_dict = defaultdict(lambda: defaultdict(int))
    for alignment_result in clade_alignment_result:
        alignment = alignment_result.vertex_alignment
        for key, value in alignment.items():
            vertex_candidate_count_dict[key][value] += 1
    return vertex_candidate_count_dict


def save_coalescent_vertex_pedigree_candidates(statistics_file, clade_alignment_result: [FullAlignmentResult],
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


def get_vertex_alignments_normalized_probabilities(alignments: list[dict], tree: CoalescentTree,
                                                   pedigree: PotentialMrcaProcessedGraph):
    alignment_likelihoods = []
    for alignment in alignments:
        alignment_likelihood = get_vertex_alignment_estimated_length(
            coalescent_tree=tree,
            pedigree=pedigree,
            alignment=alignment
        )
        alignment_likelihood = 2 ** (-alignment_likelihood)
        alignment_likelihoods.append(alignment_likelihood)
    likelihood_sum = sum(x for x in alignment_likelihoods)
    normalized_likelihoods = [x / likelihood_sum for x in alignment_likelihoods]
    return normalized_likelihoods


@dataclass
class CladeAlignmentStatisticsMetadata:
    # This class stores the flags indicating what statistical data is to be calculated for the results
    calculate_similarity: bool
    calculate_distances_histogram: bool
    calculate_alignments_likelihoods: bool


class CladeAlignmentMetadata:
    results_filepath: Path
    coalescent_tree: CoalescentTree
    clade_root: int
    probands: list[int] | None
    clade: list[int] | None

    def __init__(self, results_filepath, coalescent_tree, clade_root):
        self.results_filepath = results_filepath
        self.coalescent_tree = coalescent_tree
        self.clade_root = clade_root
        self.clade = None
        self.probands = None

    def save_tree_metadata(self, statistics_file):
        statistics_file.write(section_separator)
        # Calculating the clade for the given vertex and sorting the vertices by their levels
        clade = list(self.coalescent_tree.get_descendants_for_vertex(vertex_id=self.clade_root, include_self=True))
        clade = sorted(clade, key=lambda v: self.coalescent_tree.get_vertex_level(v), reverse=True)
        self.clade = clade
        # Calculating the rest of the data
        number_of_coalescing_events = len([x for x in clade if len(self.coalescent_tree.get_children(x)) > 1])
        self.probands = [x for x in clade if self.coalescent_tree.get_vertex_level(x) == 0]
        statistics_file.write(f"The root of the clade: {self.clade_root}\n")
        vertex_count = len(clade)
        verb = "is" if vertex_count == 1 else "are"
        vertex_noun = "vertex" if vertex_count == 1 else "vertices"
        statistics_file.write(f"There {verb} {vertex_count} {vertex_noun} in the clade\n")
        probands_count = len(self.probands)
        verb = "is" if probands_count == 1 else "are"
        probands_noun = "proband" if probands_count == 1 else "probands"
        statistics_file.write(f"There {verb} {probands_count} {probands_noun} in the clade\n")
        statistics_file.write(f"There are {number_of_coalescing_events} coalescing events in the clade\n")
        statistics_file.write(section_separator)
        statistics_file.write("Number of coalescing events grouped by the children number:\n")
        coalescing_number_to_events_number = dict()
        for vertex in clade:
            tree_children = self.coalescent_tree.get_children(vertex)
            children_number = len(tree_children)
            if not children_number:
                continue
            previous_counter = coalescing_number_to_events_number.get(children_number, 0)
            coalescing_number_to_events_number[children_number] = previous_counter + 1
        coalescing_numbers = sorted(coalescing_number_to_events_number.keys())
        for coalescing_number in coalescing_numbers:
            statistics_file.write(f"{coalescing_number}: {coalescing_number_to_events_number[coalescing_number]}\n")

    def save(self):
        os.makedirs(self.results_filepath, exist_ok=True)
        clade_data_filepath = self.get_clade_metadata_path()
        with open(clade_data_filepath, "a") as statistics_file:
            self.save_tree_metadata(statistics_file)

    def get_clade_metadata_path(self):
        return self.results_filepath / f"_clade_{self.clade_root}.txt"


class ClimbingFailedCladeAlignmentMetadata(CladeAlignmentMetadata):
    failed_alignment_information: FailedClimbingCladeAlignmentResults

    def __init__(self, results_filepath: str | Path, coalescent_tree: CoalescentTree, clade_root: int,
                 clade_alignment_result: FailedClimbingCladeAlignmentResults):
        super().__init__(results_filepath, coalescent_tree, clade_root)
        self.failed_alignment_information = clade_alignment_result

    def save(self):
        super().save()
        alignment_path = self.get_clade_metadata_path()
        os.makedirs(self.results_filepath, exist_ok=True)
        climbing_information = self.failed_alignment_information.failed_climbing_alignment_info
        with open(alignment_path, "a") as statistics_file:
            statistics_file.write(f"Failed climbing at {climbing_information.failed_vertex}\n")
            statistics_file.write(f"There are {len(climbing_information.child_to_pedigree_candidates)} "
                                  f"children for the failed vertex\n")
            statistics_file.write("The pedigree candidates for every child vertex:\n")
            for vertex, pedigree_candidates in climbing_information.child_to_pedigree_candidates.items():
                formatted_pedigree_candidates = [convert_ploid_id_to_individual(x) for x in pedigree_candidates]
                statistics_file.write(f"{vertex}: {formatted_pedigree_candidates}\n")


class SuccessCladeAlignmentMetadata(CladeAlignmentMetadata):
    # This class stores all the information about a successful alignment for a clade with all the
    # additional information
    # TODO: The clade root information is duplicated here and in clade_alignments_metadata
    pedigree: PotentialMrcaProcessedGraph
    clade_alignment_result: SuccessCladeAlignmentResults
    clade_alignment_statistics_metadata: CladeAlignmentStatisticsMetadata = None

    def __init__(self, results_filepath: str | Path, coalescent_tree: CoalescentTree, clade_root: int,
                 clade_alignment_result: SuccessCladeAlignmentResults,
                 pedigree: PotentialMrcaProcessedGraph,
                 clade_alignments_metadata: CladeAlignmentStatisticsMetadata = None):
        results_filepath = Path(results_filepath)
        super().__init__(results_filepath, coalescent_tree, clade_root)
        self.pedigree = pedigree
        self.clade_alignment_result = clade_alignment_result
        self.clade_alignment_statistics_metadata = clade_alignments_metadata

    def save(self):
        super().save()
        if self.clade_alignment_statistics_metadata.calculate_alignments_likelihoods:
            self.clade_alignment_result.calculate_posterior_probabilities()
        self._save_alignment_individual_metadata()
        self._save_clade_metadata()

    def _save_alignments_similarities(self, statistics_file):
        alignments = [x.vertex_alignment for x in self.clade_alignment_result.alignments]
        alignments_similarity = get_alignments_similarity(alignments, self.probands)
        alignments_ind_similarity = get_alignments_ind_similarity(alignments, self.probands)
        statistics_file.write(f"The alignments similarity: {alignments_similarity}\n")
        statistics_file.write(f"The alignments individual similarity: {alignments_ind_similarity}\n")
        statistics_file.write(section_separator)

    def _save_distance_histograms(self, statistics_file):
        min_length = sys.maxsize
        min_length_alignments = []
        # TODO: Cache the distance between two vertices per pedigree
        for counter, valid_alignment in enumerate(self.clade_alignment_result.alignments):
            valid_alignment: FullAlignmentResult
            alignment_length = get_vertex_alignment_estimated_length(pedigree=self.pedigree,
                                                                     coalescent_tree=self.coalescent_tree,
                                                                     alignment=valid_alignment.vertex_alignment)
            if calculate_distances_histogram:
                if min_length == alignment_length:
                    min_length_alignments.append(counter)
                elif alignment_length < min_length:
                    min_length = alignment_length
                    min_length_alignments = [counter]
        statistics_file.write(f"The minimum alignment length: {min_length}\n")
        statistics_file.write(f"The corresponding alignments: {min_length_alignments}\n")
        statistics_file.write(section_separator)

    def _save_basic_alignments_metadata(self, statistics_file):
        statistics_file.write(section_separator)
        alignments_number = len(self.clade_alignment_result.alignments)
        statistics_file.write(f"The total number of alignments is: {alignments_number}\n")

    def _save_alignments_posterior_probabilities(self):
        if not self.clade_alignment_result.clade_alignment_posterior_probabilities:
            return
        probability_dictionary = self.clade_alignment_result.clade_alignment_posterior_probabilities.vertex_alignment_to_posterior_probability
        for index, alignment_likelihood in probability_dictionary.items():
            assert alignment_likelihood >= 0.0
            rounded_alignment_likelihood = round_down(alignment_likelihood)
            alignment_filename = f"alignment_{index}"
            alignment_path = self.results_filepath / alignment_filename
            prefix_string = f"Posterior alignment likelihood: {rounded_alignment_likelihood}\n{section_separator}"
            prepend_to_file(alignment_path, prefix_string)

    def _save_alignment_individual_metadata(self):
        self._save_alignments_posterior_probabilities()

    def _save_clade_metadata(self):
        metadata_filepath = self.get_clade_metadata_path()
        with open(metadata_filepath, 'w') as statistics_file:
            self.save_tree_metadata(statistics_file)
            self._save_basic_alignments_metadata(statistics_file)
            self._save_number_of_pedigree_candidates_per_clade_vertex(statistics_file)
            if not self.clade_alignment_statistics_metadata:
                return
            if self.clade_alignment_result.clade_alignment_posterior_probabilities:
                probabilities = self.clade_alignment_result.clade_alignment_posterior_probabilities
                if probabilities.vertex_alignment_to_posterior_probability:
                    self._save_vertex_alignments_posterior_probabilities(statistics_file)
                if probabilities.vertex_posterior_probabilities:
                    self._save_pedigree_vertex_inclusion_probabilities(statistics_file)
            else:
                statistics_file.write("The posterior probabilities were not calculated\n")
            if self.clade_alignment_statistics_metadata.calculate_similarity:
                self._save_alignments_similarities(statistics_file)
            else:
                statistics_file.write("The similarities among the trees were not calculated\n")
            if self.clade_alignment_statistics_metadata.calculate_distances_histogram:
                self._save_distance_histograms(statistics_file)
            else:
                statistics_file.write("The distance histograms were not calculated\n")

    def _save_number_of_pedigree_candidates_per_clade_vertex(self, statistics_file):
        coalescent_vertex_pedigree_candidates_number = \
            {x: len({alignment.vertex_alignment[x] for alignment in
                     self.clade_alignment_result.alignments})
             for x in self.clade
             }
        statistics_file.write(section_separator)
        statistics_file.write("The number of pedigree candidates for every vertex:\n")
        tree_probands = self.coalescent_tree.get_sink_vertices()
        for coalescent_vertex, vertex_candidates_number in coalescent_vertex_pedigree_candidates_number.items():
            if coalescent_vertex in tree_probands:
                statistics_file.write(f"{coalescent_vertex} (proband): {vertex_candidates_number}\n")
            else:
                statistics_file.write(f"{coalescent_vertex}: {vertex_candidates_number}\n")
        save_coalescent_vertex_pedigree_candidates(
            statistics_file=statistics_file,
            clade_alignment_result=self.clade_alignment_result.alignments,
            clade=self.clade
        )

    def _save_vertex_alignments_posterior_probabilities(self, statistics_file):
        statistics_file.write(section_separator)
        statistics_file.write("The vertex alignments posterior probabilities:\n")
        probability_dictionary = self.clade_alignment_result.clade_alignment_posterior_probabilities.vertex_alignment_to_posterior_probability
        for index, probability in sorted(probability_dictionary.items(), key=lambda x: x[1], reverse=True):
            statistics_file.write(f"alignment_{index}: {probability}\n")

    def _save_pedigree_vertex_inclusion_probabilities(self, statistics_file):
        statistics_file.write(section_separator)
        statistics_file.write("The pedigree vertex inclusion probabilities:\n")
        probability_dictionary = self.clade_alignment_result.clade_alignment_posterior_probabilities.vertex_posterior_probabilities
        for vertex, probability in sorted(probability_dictionary.items(), key=lambda x: x[1], reverse=True):
            converted_vertex_id = convert_ploid_id_to_individual(vertex)
            statistics_file.write(f"{converted_vertex_id}: {probability}\n")
        # Save the vertex inclusion probabilities for every alignment separately
        for alignment_result in self.clade_alignment_result.alignments:
            alignment_result.save_vertex_inclusion_probabilities(
                global_vertex_inclusion_probability=probability_dictionary
            )
