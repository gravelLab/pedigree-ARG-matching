from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from itertools import chain
from pathlib import Path

import networkx as nx
from lineagekit.core.CoalescentTree import CoalescentTree
from networkx import NetworkXNoCycle, HasACycle

from alignment.configuration import *
import yaml

from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.utility.alignment_utility import convert_ploid_str_to_ploid_id
from alignment.alignment_constants import *


class YAMLValidationError(Exception):
    """Custom exception for YAML validation errors."""
    pass


@dataclass
class GraphParsingRules:
    filepath: str | Path
    missing_parent_notation: str
    separation_symbol: str
    skip_first_line: bool
    verify_graph_has_no_cycles: bool

    @staticmethod
    def parse_graph_parsing_rules(yaml_dict: dict, verify_graph_has_no_cycles: bool) -> GraphParsingRules:
        path = yaml_dict[path_key]
        separation_symbol = default_separation_symbol
        missing_parent_notation = default_missing_parent_notation
        skip_first_line = default_skip_first_line
        if separation_symbol_key in yaml_dict:
            separation_symbol = yaml_dict[separation_symbol_key]
        if missing_parent_notation_key in yaml_dict:
            missing_parent_notation = yaml_dict[missing_parent_notation_key]
        if skip_first_line_key in yaml_dict:
            skip_first_line = yaml_dict[skip_first_line_key]
        if verify_graph_has_no_cycles_key in yaml_dict:
            verify_graph_has_no_cycles = yaml_dict[verify_graph_has_no_cycles_key]
        return GraphParsingRules(
            filepath=path,
            separation_symbol=separation_symbol,
            missing_parent_notation=missing_parent_notation,
            skip_first_line=skip_first_line,
            verify_graph_has_no_cycles=verify_graph_has_no_cycles
        )


@dataclass
class ParsedDriverFile:
    pedigree_parsing_rules: GraphParsingRules
    coalescent_tree_parsing_rules: GraphParsingRules
    initial_assignments: dict[int, [int]]
    output_path: str | Path
    driver_file_path: str | Path
    vertex_alignment_mode: AlignmentVertexMode
    edge_alignment_mode: AlignmentEdgeMode
    probability_calculation_mode: PosteriorProbabilitiesCalculationMode

    def _identify_path(self, path: str | Path):
        # Firstly, verify is the path is valid for the current working directory
        if os.path.exists(path):
            return Path(path)
        # Otherwise, try appending it to the driver's absolute path. This allows the user to specify
        # a relative path regarding the driver file's location
        driver_relative_path = Path(self.driver_file_path).parent / path
        if os.path.exists(driver_relative_path):
            return driver_relative_path
        return None

    def identify_pedigree_path(self):
        return self._identify_path(self.pedigree_parsing_rules.filepath)

    def identify_coalescent_tree_path(self):
        return self._identify_path(self.coalescent_tree_parsing_rules.filepath)

    @staticmethod
    def parse_driver_file_and_validate_initial_assignments(filepath: str | Path) -> ParsedDriverFile:
        """
        Validates the specified YAML file and parses it into a dictionary mapping
        coalescent_id to an InitialAssignment object.

        Args:
            filepath (str | Path): Path to the YAML file.
        Returns:
            Dictionary where keys are coalescent_ids and values are InitialAssignment objects.
        Raises:
             YAMLValidationError: if the file is invalid.
        """
        root_required_keys = {initial_assignments_key, coalescent_tree_key, pedigree_key, output_path_key}
        initial_assignments_required_keys = {"coalescent_id", "pedigree_ids"}
        try:
            # Load the YAML file
            with open(filepath, 'r') as f:
                data = yaml.safe_load(f)

            # Check if all required keys are present
            if not root_required_keys.issubset(data.keys()):
                raise YAMLValidationError(
                    f"Missing required keys. Expected: {root_required_keys}, Found: {data.keys()}")

            # Parse the parsing data for the pedigree and the tree
            pedigree_parsing_data = data[pedigree_key]
            coalescent_tree_parsing_data = data[coalescent_tree_key]
            parsing_data = [("pedigree", pedigree_parsing_data), ("coalescent_tree", coalescent_tree_parsing_data)]
            for data_name, graph_data in parsing_data:
                if path_key not in graph_data:
                    raise YAMLValidationError(f"The path for the {data_name} is not specified")
            pedigree_parsing_rules = GraphParsingRules.parse_graph_parsing_rules(pedigree_parsing_data,
                                                                                 verify_graph_has_no_cycles=False)
            coalescent_tree_parsing_rules = GraphParsingRules.parse_graph_parsing_rules(coalescent_tree_parsing_data,
                                                                                        verify_graph_has_no_cycles=True)
            output_path = data[output_path_key]
            if not output_path:
                raise YAMLValidationError("Output path is empty")
            # Verifying if the specified path is absolute. If not, the specified path is treated as a relative path
            # with regard to the driver's file location
            if not os.path.isabs(output_path):
                output_path = Path(filepath).parent / output_path
            # Validate "initial_assignments" format
            if not isinstance(data[initial_assignments_key], list):
                raise YAMLValidationError("initial_assignments should be a list.")
            initial_assignments_dictionary = {}
            for entry in data[initial_assignments_key]:
                if not isinstance(entry, dict):
                    raise YAMLValidationError(f"Each entry in initial_assignments must be a dictionary. Found: {entry}")
                if not initial_assignments_required_keys.issubset(entry.keys()):
                    raise YAMLValidationError(f"Missing keys in initial_assignments entry. Found: {entry.keys()}")

                # Validate coalescent_id uniqueness
                try:
                    coalescent_id = int(entry[coalescent_id_key])
                except ValueError:
                    raise YAMLValidationError(f"Coalescent id must be an integer: {entry[coalescent_id_key]}")

                # Check if coalescent_id already exists
                if coalescent_id in initial_assignments_dictionary:
                    raise YAMLValidationError(f"Duplicate coalescent_id found: {coalescent_id}")

                pedigree_unprocessed_ids = entry[pedigree_ids_key]
                pedigree_processed_ids = []
                for pedigree_id in pedigree_unprocessed_ids:
                    try:
                        processed_pedigree_id = convert_ploid_str_to_ploid_id(pedigree_id)
                    except ValueError as e:
                        raise YAMLValidationError(e)
                    pedigree_processed_ids.append(processed_pedigree_id)
                if not pedigree_processed_ids:
                    raise YAMLValidationError(f"No pedigree ids specified for {coalescent_id}")
                # Add to the dictionary
                initial_assignments_dictionary[coalescent_id] = pedigree_processed_ids
            vertex_alignment_mode = AlignmentVertexMode.ALL_ALIGNMENTS
            edge_alignment_mode = AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT
            probability_calculation_mode = PosteriorProbabilitiesCalculationMode.SKIP
            if alignment_vertex_mode_key in data:
                specified_alignment_mode = data[alignment_vertex_mode_key]
                if specified_alignment_mode not in alignment_vertex_mode_dict:
                    raise YAMLValidationError(f"Invalid alignment mode specified: \"{specified_alignment_mode}\"")
                vertex_alignment_mode = alignment_vertex_mode_dict[specified_alignment_mode]
            if alignment_edge_mode_key in data:
                specified_alignment_mode = data[alignment_edge_mode_key]
                if specified_alignment_mode not in alignment_edge_mode_dict:
                    raise YAMLValidationError(f"Invalid alignment mode specified: \"{specified_alignment_mode}\"")
                edge_alignment_mode = alignment_edge_mode_dict[specified_alignment_mode]
            if posterior_probability_calculation_mode_key in data:
                specified_calculation_mode = data[posterior_probability_calculation_mode_key]
                if specified_calculation_mode not in posterior_probability_calculation_mode_dict:
                    raise YAMLValidationError(f"Invalid probability calculation mode"
                                              f" specified: \"{specified_calculation_mode}\"")
                probability_calculation_mode = posterior_probability_calculation_mode_dict[specified_calculation_mode]
                if (probability_calculation_mode != PosteriorProbabilitiesCalculationMode.SKIP and
                        (edge_alignment_mode != AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS or
                         vertex_alignment_mode != AlignmentVertexMode.ALL_ALIGNMENTS)):
                    raise YAMLValidationError(f"Probabilities can be calculated only when all the vertex and edge "
                                              f"alignments are calculated.")
            if (vertex_alignment_mode == AlignmentVertexMode.EXAMPLE_PER_ROOT_ASSIGNMENT and
                    edge_alignment_mode == AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS):
                warnings.warn("It is unusual to calculate only example vertex alignments, but with the list"
                              " of all the edge alignments. Consider double-checking your driver file. The"
                              " inference will continue with the specified options")
            return ParsedDriverFile(pedigree_parsing_rules=pedigree_parsing_rules,
                                    coalescent_tree_parsing_rules=coalescent_tree_parsing_rules,
                                    initial_assignments=initial_assignments_dictionary,
                                    driver_file_path=filepath,
                                    output_path=output_path,
                                    vertex_alignment_mode=vertex_alignment_mode,
                                    edge_alignment_mode=edge_alignment_mode,
                                    probability_calculation_mode=probability_calculation_mode,
                                    )
        except yaml.YAMLError as e:
            raise YAMLValidationError(f"Error parsing YAML file: {e}")


@dataclass
class ProcessedDriverFile:
    pedigree: PotentialMrcaProcessedGraph
    coalescent_tree: CoalescentTree
    output_path: str | Path
    initial_assignments: dict[int, [int]]
    alignment_vertex_mode: AlignmentVertexMode
    alignment_edge_mode: AlignmentEdgeMode
    probability_calculation_mode: PosteriorProbabilitiesCalculationMode

    def preprocess_graphs_for_alignment(self):
        self.preprocess_pedigree()
        self.preprocess_coalescent_tree()

    def get_pedigree_probands_for_alignment(self):
        return list(chain.from_iterable(self.initial_assignments.values()))

    def preprocess_pedigree(self):
        pedigree_probands = self.get_pedigree_probands_for_alignment()
        self.pedigree.reduce_to_ascending_graph(probands=pedigree_probands)
        self.pedigree.get_levels()
        self.pedigree.initialize_potential_mrca_map()

    def preprocess_coalescent_tree(self):
        self.coalescent_tree.get_levels()
        self.coalescent_tree.remove_unary_nodes()

    @staticmethod
    def finish_driver_file_processing(parsed_driver_file: ParsedDriverFile) -> ProcessedDriverFile:
        pedigree_path = parsed_driver_file.identify_pedigree_path()
        if not pedigree_path:
            raise ValueError(f"The specified pedigree file cannot be found: "
                             f"{parsed_driver_file.pedigree_parsing_rules.filepath}")
        coalescent_tree_path = parsed_driver_file.identify_coalescent_tree_path()
        if not coalescent_tree_path:
            raise ValueError(f"The specified coalescent tree file cannot be found:"
                             f" {parsed_driver_file.coalescent_tree_parsing_rules.filepath}")
        # Parse the coalescent tree
        tree_missing_parent_notation = parsed_driver_file.coalescent_tree_parsing_rules.missing_parent_notation
        tree_separation_symbol = parsed_driver_file.coalescent_tree_parsing_rules.separation_symbol
        tree_skip_first_line = parsed_driver_file.coalescent_tree_parsing_rules.skip_first_line
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(
            filepath=coalescent_tree_path,
            missing_parent_notation=[tree_missing_parent_notation],
            separation_symbol=tree_separation_symbol,
            skip_first_line=tree_skip_first_line
        )
        if parsed_driver_file.coalescent_tree_parsing_rules.verify_graph_has_no_cycles:
            try:
                cycle = nx.find_cycle(coalescent_tree, orientation='original')
                raise HasACycle(f"The coalescent tree has a cycle: {cycle}")
            except NetworkXNoCycle:
                # No errors are found
                pass
        initial_assignments = parsed_driver_file.initial_assignments
        specified_probands = initial_assignments.keys()
        coalescent_tree.reduce_to_ascending_graph(probands=specified_probands)
        # Parse the pedigree
        pedigree_missing_parent_notation = parsed_driver_file.pedigree_parsing_rules.missing_parent_notation
        pedigree_separation_symbol = parsed_driver_file.pedigree_parsing_rules.separation_symbol
        pedigree_skip_first_line = parsed_driver_file.pedigree_parsing_rules.skip_first_line
        pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(
            filepath=pedigree_path,
            preprocess_graph=False,
            missing_parent_notation=[pedigree_missing_parent_notation],
            separation_symbol=pedigree_separation_symbol,
            skip_first_line=pedigree_skip_first_line
        )
        if parsed_driver_file.pedigree_parsing_rules.verify_graph_has_no_cycles:
            try:
                cycle = nx.find_cycle(pedigree, orientation='original')
                raise HasACycle(f"The pedigree has a cycle: {cycle}")
            except NetworkXNoCycle:
                # No errors are found
                pass
        # Calculate the leaf vertices in the coalescent tree for which mapping isn't specified
        coalescent_tree_probands = coalescent_tree.get_sink_vertices()
        for coalescent_vertex, pedigree_vertices in initial_assignments.items():
            if coalescent_vertex not in coalescent_tree_probands:
                raise ValueError(f"The specified coalescent vertex {coalescent_vertex} either does not exist or isn't"
                                 f" a leaf vertex")
            for pedigree_vertex in pedigree_vertices:
                if pedigree_vertex not in pedigree:
                    raise ValueError(f"The specified pedigree vertex {pedigree_vertex} does not exist")
        return ProcessedDriverFile(
            coalescent_tree=coalescent_tree,
            pedigree=pedigree,
            initial_assignments=initial_assignments,
            output_path=parsed_driver_file.output_path,
            alignment_vertex_mode=parsed_driver_file.vertex_alignment_mode,
            alignment_edge_mode=parsed_driver_file.edge_alignment_mode,
            probability_calculation_mode=parsed_driver_file.probability_calculation_mode
        )

    @staticmethod
    def process_driver_file(filepath: str | Path) -> ProcessedDriverFile:
        parsed_driver_file = ParsedDriverFile.parse_driver_file_and_validate_initial_assignments(filepath=filepath)
        return ProcessedDriverFile.finish_driver_file_processing(parsed_driver_file=parsed_driver_file)
