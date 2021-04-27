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
acids_sequence6 = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'  # Best: -28
acids_sequence7 = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
acids_sequence8 = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH'  # Best: -55
acids_sequence9 = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH'


def timer(fn):
    """Generic decorator to test performance."""
    import time

    @wraps(fn)
    def wrapper_timer(*args, **kwargs):
        start_time = time.time()
        result = fn(*args, *kwargs)
        end_time = time.time()
        execution_time = round(end_time - start_time, 3)
        print(f'{fn.__name__} execution time: {execution_time} seconds')
        return result
    return wrapper_timer

@timer
def best_of(protein, n, m, folding=None):
    """Take a protein. Give it a random folding. Try to fold the protein m
    times to improve this folding. Print the resulting score. Repeat this n
    times. Return the best result.
    """
    best_score = 0
    for _ in range(n):
        if folding == 'cube_folding':
            if not Algorithms.cube_folding(protein, shift='', d3=True):
                print("failed to get a cube folding as start")
                continue
        elif folding == 'random_folding':
            if not Algorithms.random_folding(protein):
                print("failed to get a random folding as start")
                continue
        if not Algorithms.fold_n_times(m, protein, max_tolerance=5):
            print("failed to fold n times")
        score = Algorithms.score(protein)
        if score < best_score:
            best_score = score
            best_result = [acid.copy() for acid in protein.acids]
        print('score: {}, best_score: {}'.format(score, best_score))
        protein.acids = best_result
    return


protein = Protein(acids_sequence8)
best_of(protein, 50, 1000, folding='random_folding')
print('score van best_folding', Algorithms.score(protein))
ProteinPlotter.plot(protein)

# # Cube folding, while shifting the sequence through the zigzag
# for n in range(3):
#     protein = Protein(acids_sequence8)
#     Algorithms.cube_folding(protein, shift=n*'1', d3=True)
#     print('score:', Algorithms.score(protein))
#     ProteinPlotter.plot(protein)

