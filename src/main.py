from protein import Protein
from algorithms import Algorithms
from protein_plotter import ProteinPlotter
from functools import wraps

# Without C's
protein_string1 = 'HHPHHHPH'
protein_string2 = 'HHPHHHPHPHHHPH'
protein_string3 = 'HPHPPHHPHPPHPHHPPHPH'  # Best: -11
protein_string4 = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP'
protein_string5 = 'HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH'  # Best: -23

# With C's
protein_string6 = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'  # Best: -28
protein_string7 = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
protein_string8 = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH'  # Best: -55
protein_string9 = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH'


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
def best_of(protein, n, m):
    """Take a protein. Give it a random folding. Try to fold the protein m
    times to improve this folding. Print the resulting score. Repeat this n
    times. Return the best result.
    """
    best_score = 0
    for _ in range(n):
        if not Algorithms.random_folding(protein):
            print("oke let's continue")
            continue
        Algorithms.fold_n_times(m, protein)
        score = Algorithms.score(protein)
        if score < best_score:
            best_score = score
            best_result = protein
            # print('score: {}, best_score: {}'.format(score, best_score))
    print('best_score: {}'.format(best_score))
    return best_result


protein = Protein(protein_string8)
best_folding = best_of(protein, 1000, 200)
ProteinPlotter.plot(best_folding)

