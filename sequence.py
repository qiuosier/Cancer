import numpy as np
import editdistance


COMPLEMENTS = {
    "A": "T",
    "T": "A",
    "U": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


# The levenshtein function is copied from
# https://stackabuse.com/levenshtein-distance-and-text-similarity-in-python/
def levenshtein(s1, s2):
    size_x = len(s1) + 1
    size_y = len(s2) + 1
    matrix = np.zeros((size_x, size_y))
    for x in range(size_x):
        matrix[x, 0] = x
    for y in range(size_y):
        matrix[0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if s1[x - 1] == s2[y - 1]:
                matrix[x, y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix[x, y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1] + 1,
                    matrix[x, y-1] + 1
                )
    return matrix[size_x - 1, size_y - 1]


def semi_global(s1, s2):
    pass


def hamming(s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    mismatch = 0
    for i in range(max(l1, l2)):
        c_1 = s1[i] if i < l1 else ""
        c_2 = s2[i] if i < l2 else ""
        if c_1 != c_2:
            mismatch += 1
    return mismatch


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.length = len(sequence)

    def __len__(self):
        return self.length

    @property
    def complements(self):
        return "".join([COMPLEMENTS[c] for c in self.sequence])

    @property
    def reverse_complements(self):
        return self.complements[::-1]

    def match(self, sequence):
        s1 = self.sequence.upper()
        s2 = sequence.upper()
        # return hamming(s1, s2)
        return editdistance.eval(s1, s2)
