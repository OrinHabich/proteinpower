import sys
sys.path.append("..")
from src.protein import Protein


def test_no_folding(acids_sequence, helpers):
    protein = Protein(acids_sequence)
    helpers.assert_correct_order(protein, acids_sequence)


def test_fold(protein, acids_sequence, helpers):
    protein.fold(3)
    helpers.assert_correct_order(protein, acids_sequence)


def test_own_folding(acids_sequence):
    own_folding = [{'x': 1, 'y': 0, 'z': 0, 'acid_type': 'H'},
                   {'x': 1, 'y': 0, 'z': 0, 'acid_type': 'H'},
                   {'x': 2, 'y': 0, 'z': 0, 'acid_type': 'P'},
                   {'x': 3, 'y': 0, 'z': 0, 'acid_type': 'H'},
                   {'x': 4, 'y': 0, 'z': 0, 'acid_type': 'H'},
                   {'x': 5, 'y': 0, 'z': 0, 'acid_type': 'H'},
                   {'x': 6, 'y': 0, 'z': 0, 'acid_type': 'P'},
                   {'x': 6, 'y': 1, 'z': 0, 'acid_type': 'H'}]
    protein = Protein(acids_sequence, folding=own_folding)
    assert protein.acids, own_folding
