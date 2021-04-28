import sys
sys.path.append("..")
from src.algorithms import Algorithms


def test_score(protein):
    assert Algorithms.score(protein) == -1


def test_find_bonds(protein):
    weak_bonds, _ = Algorithms.find_bonds(protein.acids)
    assert len(weak_bonds) == 1


def test_fold_100_times(protein, acids_sequence, helpers):
    Algorithms.fold_n_times(100, protein)
    helpers.assert_correct_order(protein, acids_sequence)


def test_random_folding(protein, acids_sequence, helpers):
    Algorithms.random_folding(protein)
    helpers.assert_correct_order(protein, acids_sequence)
    assert sum(protein.acids[0][c] == 0 for c in 'xyz') == 3
    # assert protein.acids[0]['x'] == 0
    # assert protein.acids[0]['y'] == 0
    # assert protein.acids[0]['z'] == 0
    assert protein._injective()

def test_cube_folding(protein, acids_sequence, helpers):
    Algorithms.cube_folding(protein)
    helpers.assert_correct_order(protein, acids_sequence)
    assert protein._injective()
