import itertools
import random
import math


class Algorithms():
    """Different algorithms to fold a protein. As well as a score-function."""

    @classmethod
    def score(cls, protein):
        """Put a number to a protein folding to indicate its quality. This
        enables comparison of protein foldings. The lower the number, the
        better.
        """
        weak_bonds, strong_bonds = cls.find_bonds(protein.acids)
        return len(weak_bonds) * -1 + len(strong_bonds) * -5

    @classmethod
    def find_bonds(cls, acids):
        """Find all bonds in a folding."""
        relevant_acids = [acid for acid in acids if acid['acid_type'] != 'P']
        weak_bonds = []
        strong_bonds = []
        for acid1, acid2 in itertools.combinations(relevant_acids, 2):
            if abs(acids.index(acid1) - acids.index(acid2)) > 1 and\
               cls._distance(acid1, acid2) == 1:
                if acid1['acid_type'] == 'H' or acid2['acid_type'] == 'H':
                    weak_bonds.append([acid1, acid2])
                else:
                    strong_bonds.append([acid1, acid2])
        return weak_bonds, strong_bonds

    @classmethod
    def fold_n_times(cls, n, protein, max_tolerance=0):
        """Fold the given protein n times. Stop if folding is not possible.
        Note that this is basically a hillclimber.
        """
        highscore = cls.score(protein) or 1
        for i in range(n):
            indices_possible = list(range(1, len(protein.acids)-1))
            random.shuffle(indices_possible)
            success = False
            while not success:
                tolerance = 0
                previous_folding = [acid.copy() for acid in protein.acids]
                while tolerance <= max_tolerance:
                    index = indices_possible[0]
                    if protein.fold(index):
                        score = cls.score(protein)
                        if score > highscore:
                            protein.acids = previous_folding
                        else:
                            highscore = score
                        success = True
                        break
                    else:
                        tolerance += 1
                        del indices_possible[0]
                        if not indices_possible:
                            return False
                if not success:
                    protein.acids = previous_folding
        return True

    @classmethod
    def random_folding(cls, protein):
        """Return a random folding of the protein."""
        acids_sequence = protein.acids_sequence
        previous_acid = {'x': 0,
                         'y': 0,
                         'z': 0,
                         'acid_type': acids_sequence[0]}
        result = [previous_acid]
        for type_acid in acids_sequence[1::]:
            possible = [('x', -1), ('x',  1),
                        ('y', -1), ('y',  1),
                        ('z', -1), ('z',  1)]
            random.shuffle(possible)
            success = False
            while not success:
                axis, direction = possible[0]
                acid = previous_acid.copy()
                acid[axis] += direction
                acid['acid_type'] = type_acid
                success = not cls._same_position(result, acid)
                if not success:
                    del possible[0]
                    if not possible:
                        return False
            previous_acid = acid
            result.append(acid)
        protein.acids = result
        return True

    @classmethod
    def cube_folding(cls, protein, d3=True, shift=''):
        """Fold the given protein in a zigzagging manner into a rectangular
        shape (if d3 is false) or a cube (if d3 is true). With the shift
        keyword the foldings of the protein shift the length of the shift
        argument."""
        acid_types = list(shift + protein.acids_sequence)
        if d3:
            d = len(acid_types)**(1/3)
        else:
            d = math.sqrt(len(acid_types))
        width = list(range(math.floor(d)))
        length = list(range(math.ceil(d)))
        z = 0
        result = []
        while True:
            for y in length:
                for x in width:
                    acid = {'x': x,
                            'y': y,
                            'z': z,
                            'acid_type': acid_types.pop(0)}
                    result.append(acid)
                    if not acid_types:
                        protein.acids = result[len(shift):]
                        return True
                width = width[::-1]
            length = length[::-1]
            z += 1
        return False

    @staticmethod
    def _same_position(folded_part, new_acid):
        """Return True iff the position of some new acid is already taken."""
        for previous_acid in folded_part:
            if previous_acid['x'] == new_acid['x'] and\
               previous_acid['y'] == new_acid['y'] and\
               previous_acid['z'] == new_acid['z']:
                return True
        return False

    @staticmethod
    def _distance(point1, point2):
        """Return euclidean distance between two points."""
        sum_squares = 0
        for dimension in 'xyz':
            sum_squares += abs(point1[dimension] - point2[dimension])**2
        return math.sqrt(sum_squares)
