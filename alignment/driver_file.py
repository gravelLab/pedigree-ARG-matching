from __future__ import annotations

import os
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from alignment.configuration import *
import yaml
from alignment.potential_mrca_processed_graph import *
from graph.coalescent_tree import CoalescentTree


class YAMLValidationError(Exception):
    """Custom exception for YAML validation errors."""
    pass


class PloidType(Enum):
    Paternal = "P"
    Maternal = "M"


@dataclass
class GraphParsingRules:
    filepath: str | Path
    missing_parent_notation: str
    separation_symbol: str
    skip_first_line: bool

    @staticmethod
    def parse_graph_parsing_rules(yaml_dict: dict) -> GraphParsingRules:
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
        return GraphParsingRules(
            filepath=path,
            separation_symbol=separation_symbol,
            missing_parent_notation=missing_parent_notation,
            skip_first_line=skip_first_line
        )


@dataclass
class ParsedDriverFile:
    pedigree_parsing_rules: GraphParsingRules
    coalescent_tree_parsing_rules: GraphParsingRules
    initial_assignments: dict[int, [int]]
    output_path: str | Path
    driver_file_path: str | Path
    alignment_mode: MatchingMode

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
            pedigree_parsing_rules = GraphParsingRules.parse_graph_parsing_rules(pedigree_parsing_data)
            coalescent_tree_parsing_rules = GraphParsingRules.parse_graph_parsing_rules(coalescent_tree_parsing_data)
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
                    raise YAMLValidationError(f"Invalid coalescent id: {entry[coalescent_id_key]}")

                # Check if coalescent_id already exists
                if coalescent_id in initial_assignments_dictionary:
                    raise YAMLValidationError(f"Duplicate coalescent_id found: {coalescent_id}")

                pedigree_unprocessed_ids = entry[pedigree_ids_key]
                pedigree_processed_ids = []
                for pedigree_id in pedigree_unprocessed_ids:
                    if not pedigree_id:
                        raise YAMLValidationError(f"Pedigree id is empty")
                    try:
                        ploid_type = PloidType(pedigree_id[-1])
                    except ValueError:
                        raise YAMLValidationError(
                            f"Invalid ploid_type: {entry['ploid_type']}. Must be one of {[e.value for e in PloidType]}"
                        )
                    try:
                        unprocessed_pedigree_id = int(pedigree_id[:-1])
                    except ValueError:
                        raise YAMLValidationError(f"Invalid pedigree id {pedigree_id[:-1]}")
                    if unprocessed_pedigree_id < 0:
                        raise YAMLValidationError(f"Negative pedigree id {unprocessed_pedigree_id}")
                    if ploid_type == PloidType.Paternal:
                        processed_pedigree_id = 2 * unprocessed_pedigree_id
                    else:
                        processed_pedigree_id = 2 * unprocessed_pedigree_id + 1
                    pedigree_processed_ids.append(processed_pedigree_id)
                if not pedigree_processed_ids:
                    raise YAMLValidationError(f"No pedigree ids specified for {coalescent_id}")
                # Add to the dictionary
                initial_assignments_dictionary[coalescent_id] = pedigree_processed_ids
            alignment_mode = MatchingMode.ALL_ALIGNMENTS
            if alignment_mode_key in data:
                specified_alignment_mode = data[alignment_mode_key]
                if specified_alignment_mode not in alignment_mode_dict:
                    raise YAMLValidationError("Invalid alignment mode specified")
                alignment_mode = alignment_mode_dict[specified_alignment_mode]
            return ParsedDriverFile(pedigree_parsing_rules=pedigree_parsing_rules,
                                    coalescent_tree_parsing_rules=coalescent_tree_parsing_rules,
                                    initial_assignments=initial_assignments_dictionary,
                                    driver_file_path=filepath,
                                    output_path=output_path,
                                    alignment_mode=alignment_mode
                                    )
        except yaml.YAMLError as e:
            raise YAMLValidationError(f"Error parsing YAML file: {e}")


@dataclass
class ProcessedDriverFile:
    pedigree: PotentialMrcaProcessedGraph
    coalescent_tree: CoalescentTree
    output_path: str | Path
    initial_assignments: dict[int, [int]]
    alignment_mode: MatchingMode

    def preprocess_graphs_for_alignment(self):
        self.preprocess_pedigree()
        self.preprocess_coalescent_tree()

    def get_pedigree_probands_for_alignment(self):
        return list(chain.from_iterable(self.initial_assignments.values()))

    def preprocess_pedigree(self):
        pedigree_probands = self.get_pedigree_probands_for_alignment()
        self.pedigree.reduce_to_ascending_genealogy(probands=pedigree_probands, recalculate_levels=True)
        self.pedigree.initialize_vertex_to_level_map()
        self.pedigree.initialize_potential_mrca_map()

    def preprocess_coalescent_tree(self):
        self.coalescent_tree.initialize_vertex_to_level_map()
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
            initialize_levels=False,
            missing_parent_notation=[tree_missing_parent_notation],
            separation_symbol=tree_separation_symbol,
            skip_first_line=tree_skip_first_line
        )
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
        initial_assignments = parsed_driver_file.initial_assignments
        # Calculate the leaf vertices in the coalescent tree for which mapping isn't specified
        coalescent_tree_probands = coalescent_tree.get_probands()
        tree_leaf_vertices_missing_assignments = set(initial_assignments.keys()).difference(coalescent_tree_probands)
        if tree_leaf_vertices_missing_assignments:
            raise YAMLValidationError("Missing initial assignments for tree's leaf vertices "
                                      f"{tree_leaf_vertices_missing_assignments}")
        for coalescent_vertex, pedigree_vertices in initial_assignments.items():
            if coalescent_vertex not in coalescent_tree_probands:
                raise ValueError(f"The specified coalescent vertex {coalescent_vertex} either does not exist or isn't"
                                 f"a leaf vertex")
            for pedigree_vertex in pedigree_vertices:
                if pedigree_vertex not in pedigree:
                    raise ValueError(f"The specified pedigree vertex {pedigree_vertex} does not exist")
        return ProcessedDriverFile(
            coalescent_tree=coalescent_tree,
            pedigree=pedigree,
            initial_assignments=initial_assignments,
            output_path=parsed_driver_file.output_path,
            alignment_mode=parsed_driver_file.alignment_mode
        )

    @staticmethod
    def process_driver_file(filepath: str | Path) -> ProcessedDriverFile:
        parsed_driver_file = ParsedDriverFile.parse_driver_file_and_validate_initial_assignments(filepath=filepath)
        return ProcessedDriverFile.finish_driver_file_processing(parsed_driver_file=parsed_driver_file)


# YAML keys used in the driver file
initial_assignments_key = "initial_assignments"
coalescent_tree_key = "coalescent_tree"
pedigree_key = "pedigree"
coalescent_id_key = "coalescent_id"
path_key = "path"
pedigree_ids_key = "pedigree_ids"
separation_symbol_key = "separation_symbol"
missing_parent_notation_key = "missing_parent_notation"
skip_first_line_key = "skip_first_line"
output_path_key = "output_path"
alignment_mode_key = "alignment_mode"

alignment_mode_dict = {
    "default": MatchingMode.ALL_ALIGNMENTS,
    "example_per_root_assignment": MatchingMode.EXAMPLE_PER_ROOT_ASSIGNMENT
}
