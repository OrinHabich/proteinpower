from protein import Protein
from algorithms import Algorithms
from protein_plotter import ProteinPlotter
from functools import wraps

# Without C's
acids_sequence1 = 'HHPHHHPH'
acids_sequence2 = 'HHPHHHPHPHHHPH'
acids_sequence3 = 'HPHPPHHPHPPHPHHPPHPH'  # Best: -11
acids_sequence4 = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP'
acids_sequence5 = 'HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH'  # Best: -23

# With C's
acids_sequence6 = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'  # Best: -31
acids_sequence7 = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
acids_sequence8 = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH'  # Best: -40
acids_sequence9 = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH'


def best_of_experiments(acids_sequence, n, m, folding=None):
    """Take a sequence of acids. Give it a certain folding. Try to fold the
    protein m times to improve this folding. Print the resulting score. Repeat
    this n times. Return the best result.
    """
    print('Start folding:', folding)
    best_score = 1
    for _ in range(n):
        protein = Protein(acids_sequence)
        if folding == 'cube_folding':
            if not Algorithms.cube_folding(protein, shift='', d3=True):
                print("failed to get a cube folding as start")
                continue
        elif folding == 'random_folding':
            if not Algorithms.random_folding(protein):
                print("failed to get a random folding as start")
                continue
        start_score = Algorithms.score(protein)
        if not Algorithms.fold_n_times(m, protein):
            print("failed to fold n times")
        end_score = Algorithms.score(protein)
        if end_score < best_score:
            best_score = end_score
            best_result = [acid.copy() for acid in protein.acids]
        print('Start with score:\t{}\t\tEnd with score:\t{}\t\tBest score:\t{}'.format(start_score, end_score, best_score))
    protein.acids = best_result
    return protein

protein = best_of_experiments(acids_sequence6, 50, 1000, 'random_folding')
ProteinPlotter.plot(protein)



# protein = Protein(acids_sequence8)
# Algorithms.cube_folding(protein, shift='', d3=True)
# print('score of cube_folding:', Algorithms.score(protein))
# ProteinPlotter.plot(protein)

