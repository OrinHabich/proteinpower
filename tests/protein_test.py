import sys
sys.path.append("..")
from src.protein import Protein


def test_fold(acids_sequence, helpers):
    protein = Protein(acids_sequence)
    end_point_before_folding = protein.acids[-1].copy()
    protein.fold(1)
    end_point_after_folding = protein.acids[-1]
    assert end_point_before_folding != end_point_after_folding
    helpers.assert_correct_order(protein, acids_sequence)
