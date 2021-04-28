import random
import itertools
import numpy as np


class Protein():
    """A string of acids that can be folded on 90 degree angles."""
    def __init__(self, acids_sequence, folding=None):
        self.acids_sequence = acids_sequence
        self.acids = folding or self._no_folding()

    def fold(self, index):
        """Attempt to fold this protein at the given index in a valid manner
        (e.i. preserving injectivity). Return True when success, False when
        possibilities are exhausted.
        """
        collision = True
        attempts = 0
        rotations = self._useful_rotations(index)
        while collision:
            rotation_choice = rotations[0]
            M = self._rotation_matrices()[rotation_choice]
            for acid in self.acids[index+1:]:
                self._rotate(acid, M, self.acids[index])
            collision = not self._injective()
            if collision:
                del rotations[0]
                if not rotations:
                    return False
        return True

    def _useful_rotations(self, index):
        """Obtain all possible rotations and filter out the unuseful."""
        possible_rotations = list(self._rotation_matrices().keys())
        acid = self.acids[index]
        acid_before = self.acids[index - 1]
        not_useful = [c for c in 'xyz' if acid[c] != acid_before[c]][0]
        useful_rotations = list(filter(lambda r: r[0] != not_useful,
                                        possible_rotations))
        random.shuffle(useful_rotations)
        return useful_rotations

    def _injective(self):
        """ Assert this protein is injective."""
        for acid_i, acid_j in itertools.combinations(self.acids, 2):
            if acid_i['x'] == acid_j['x'] and\
               acid_i['y'] == acid_j['y'] and\
               acid_i['z'] == acid_j['z']:
                return False
        return True

    @staticmethod
    def _rotate(acid, M, basepoint):
        """Rotate a part (`acid`) of this protein."""
        v = [acid[c] - basepoint[c] for c in 'xyz']
        v_new = M.dot(v)
        for i, c in enumerate('xyz'):
            acid[c] = v_new[i] + basepoint[c]

    @staticmethod
    def _rotation_matrices():
        """Return a dictionary of 3X3 rotation matrices."""
        return {'x_clockwise': np.array([[1, 0, 0],
                                         [0, 0, 1],
                                         [0, -1, 0]]),
                'x_counterclockwise': np.array([[1, 0, 0],
                                                [0, 0, -1],
                                                [0, 1, 0]]),
                'y_clockwise': np.array([[0, 0, -1],
                                         [0, 1, 0],
                                         [1, 0, 0]]),
                'y_counterclockwise': np.array([[0, 0, 1],
                                                [0, 1, 0],
                                                [-1, 0, 0]]),
                'z_clockwise': np.array([[0, 1, 0],
                                         [-1, 0, 0],
                                         [0, 0, 1]]),
                'z_counterclockwise': np.array([[0, -1, 0],
                                                [1, 0, 0],
                                                [0, 0, 1]])}

    def _no_folding(self):
        """Return a straight protein-string."""
        previous_acid = {'x': 0,
                         'y': 0,
                         'z': 0,
                         'acid_type': self.acids_sequence[0]}
        result = [previous_acid]
        for type_acid in self.acids_sequence[1::]:
            acid = previous_acid.copy()
            acid['acid_type'] = type_acid
            acid['x'] += 1
            previous_acid = acid
            result.append(acid)
        return result


if __name__ == '__main__':
    pass
