import sys
sys.path.append("..")
from src.algorithms import Algorithms


def test_score(protein):
    assert Algorithms.score(protein) == -1


def test_find_connections(protein):
    weak_bonds, _ = Algorithms.find_connections(protein.acids)
    assert len(weak_bonds) == 1


def test_fold_100_times(protein, acids_string, helpers):
    Algorithms.fold_n_times(100, protein)
    helpers.assert_correct_order(protein, acids_string)


def test_random_folding(protein, acids_string, helpers):
    Algorithms.random_folding(protein)
    helpers.assert_correct_order(protein, acids_string)
    assert protein.acids[0]['x'] == 0
    assert protein.acids[0]['y'] == 0
    assert protein.acids[0]['z'] == 0
    assert protein._injective()
