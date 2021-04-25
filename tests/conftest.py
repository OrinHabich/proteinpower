import pytest
import sys
sys.path.append("..")
from src.protein import Protein


class Helpers:
    """Helper methods for tests."""
    @staticmethod
    def assert_correct_order(p, protein_acids):
        assert len(p.acids) == len(protein_acids)
        for i, acid_type in enumerate(protein_acids):
            assert p.acids[i]['acid_type'], acid_type


@pytest.fixture
def helpers():
    return Helpers


@pytest.fixture
def acids_string():
    return 'HHPHHHPH'


@pytest.fixture
def protein():
    """
        P___H
         |  |
    H___H|  |H__H__P__H
    """

    folding_score_1 = [{'x': 0, 'y': 0, 'z': 0, 'acid_type': 'H'},
                       {'x': 1, 'y': 0, 'z': 0, 'acid_type': 'H'},
                       {'x': 1, 'y': 1, 'z': 0, 'acid_type': 'P'},
                       {'x': 2, 'y': 1, 'z': 0, 'acid_type': 'H'},
                       {'x': 2, 'y': 0, 'z': 0, 'acid_type': 'H'},
                       {'x': 3, 'y': 0, 'z': 0, 'acid_type': 'H'},
                       {'x': 4, 'y': 0, 'z': 0, 'acid_type': 'P'},
                       {'x': 5, 'y': 1, 'z': 0, 'acid_type': 'H'}]

    return Protein('HHPHHHPH', folding=folding_score_1)
