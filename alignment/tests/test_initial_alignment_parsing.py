from pathlib import Path

import pytest

from alignment.graph_matcher import validate_and_parse_yaml, YAMLValidationError


@pytest.fixture
def correct_initial_alignment_directory_name() -> Path:
    return Path("initial_alignments")


@pytest.fixture
def correct_initial_alignment_filepath(correct_initial_alignment_directory_name):
    return correct_initial_alignment_directory_name / "example_initial_alignment.yaml"


@pytest.fixture
def invalid_alignment_negative_id_filepath(correct_initial_alignment_directory_name):
    return correct_initial_alignment_directory_name / "negative_pedigree_id.yaml"


@pytest.fixture
def invalid_alignment_invalid_pedigree_id_filepath(correct_initial_alignment_directory_name):
    return correct_initial_alignment_directory_name / "invalid_pedigree_id.yaml"


@pytest.fixture
def invalid_alignment_missing_pedigree_ids_filepath(correct_initial_alignment_directory_name):
    return correct_initial_alignment_directory_name / "missing_pedigree_ids.yaml"


@pytest.fixture
def invalid_alignment_duplicate_coalescent_ids_filepath(correct_initial_alignment_directory_name):
    return correct_initial_alignment_directory_name / "duplicate_coalescent_ids.yaml"


def test_correct_parsing(correct_initial_alignment_filepath):
    parsed_assignments = validate_and_parse_yaml(correct_initial_alignment_filepath)
    expected_assignments = {1: [22], 2: [23], 3: [24, 26, 27], 4: [29, 30, 39]}
    assert parsed_assignments.keys() == expected_assignments.keys(), \
        "The parsed alignment coalescent vertices are invalid"
    for coalescent_id, pedigree_ids in expected_assignments.items():
        parsed_pedigree_ids = parsed_assignments[coalescent_id]
        assert frozenset(pedigree_ids) == frozenset(parsed_pedigree_ids), "The parsed pedigree ids are invalid"


def test_invalid_pedigree_id(invalid_alignment_invalid_pedigree_id_filepath):
    with pytest.raises(YAMLValidationError):
        validate_and_parse_yaml(invalid_alignment_invalid_pedigree_id_filepath)


def test_negative_pedigree_id(invalid_alignment_negative_id_filepath):
    with pytest.raises(YAMLValidationError):
        validate_and_parse_yaml(invalid_alignment_negative_id_filepath)


def test_missing_pedigree_ids(invalid_alignment_missing_pedigree_ids_filepath):
    with pytest.raises(YAMLValidationError):
        validate_and_parse_yaml(invalid_alignment_missing_pedigree_ids_filepath)


def test_duplicate_coalescent_ids(invalid_alignment_duplicate_coalescent_ids_filepath):
    with pytest.raises(YAMLValidationError):
        validate_and_parse_yaml(invalid_alignment_duplicate_coalescent_ids_filepath)
