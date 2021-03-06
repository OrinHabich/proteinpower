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
            if not cls._neighbors(acids, acid1, acid2) and\
               cls._distance(acid1, acid2) == 1:
                if acid1['acid_type'] == 'H' or acid2['acid_type'] == 'H':
                    weak_bonds.append([acid1, acid2])
                else:
                    strong_bonds.append([acid1, acid2])
        return weak_bonds, strong_bonds

    @classmethod
    def fold_n_times(cls, n, protein):
        """Fold the given protein n times. Stop if folding is not possible.
        Note that this is basically a hillclimber.
        """
        highscore = cls.score(protein) or 1
        for i in range(n):
            indices_possible = list(range(1, len(protein.acids)-1))
            random.shuffle(indices_possible)
            success = False
            while not success:
                previous_folding = [acid.copy() for acid in protein.acids]
                index = indices_possible[0]
                if protein.fold(index):
                    score = cls.score(protein)
                    if score > highscore:
                        protein.acids = previous_folding
                    else:
                        highscore = score
                    success = True
                else:
                    protein.acids = previous_folding
                    del indices_possible[0]
                    if not indices_possible:
                        return False
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

    @classmethod
    def _same_position(cls, folded_part, new_acid):
        """Return True iff the position of some new acid is already taken."""
        coordinates_new_acid = cls._coordinates(new_acid)
        for acid in folded_part:
            if cls._coordinates(acid) == coordinates_new_acid:
                return True
        return False

    @classmethod
    def _neighbors(cls, acids, acid1, acid2):
        """Check if two acids are consecutive in the acids sequence."""
        return abs(acids.index(acid1) - acids.index(acid2)) == 1

    @classmethod
    def _distance(cls, acid1, acid2):
        """Euclidean distance betweeen two points."""
        return math.dist(cls._coordinates(acid1), cls._coordinates(acid2))

    @staticmethod
    def _coordinates(acid):
        """Returns list with coordinates of given acid."""
        return [v for (k, v) in acid.items() if k in 'xyz']
