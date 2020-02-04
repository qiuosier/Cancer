import numpy as np
import editdistance
import parasail


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


class Distance:
    """Contains static functions for calculating distance between two sequences.

    All functions are static.
    Each function returns a tuple of:
        0: distance, or the number of mismatches
        1: last matched position (0-indexed) in the second (reference) sequence.
            -1 if there is no match at all.
            None if this is not calculated.

    """
    # The levenshtein function is copied from
    # https://stackabuse.com/levenshtein-distance-and-text-similarity-in-python/
    @staticmethod
    def levenshtein(s1, s2):
        """This is a slow version of levenshtein distance.
        Use edit_distance() for faster calculation.
        """
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
        return matrix[size_x - 1, size_y - 1], None

    @staticmethod
    def edit_distance(s1, s2):
        return editdistance.eval(s1, s2), None

    @staticmethod
    def hamming(s1, s2):
        """Calculates the hamming distance of two sequences.
        The two sequences should have the same length.
        Otherwise, the longer one will be truncated (additional sequences are ignored).

        Args:
            s1:
            s2:

        Returns:

        """
        l1 = len(s1)
        l2 = len(s2)
        last_matched = -1
        mismatch = 0
        for i in range(min(l1, l2)):
            c_1 = s1[i] if i < l1 else ""
            c_2 = s2[i] if i < l2 else ""
            if c_1 != c_2:
                mismatch += 1
            else:
                last_matched = i
        return mismatch, last_matched

    dna_alignment_matrix = parasail.matrix_create("ACGT", 0, -1)

    @staticmethod
    def dna_alignment(query, ref, algorithm='sg_de'):
        """

        Args:
            query: query
            ref: reference
            algorithm: Function names in the algorithm table from parassail package.
            See https://github.com/jeffdaily/parasail#standard-function-naming-convention for all function names

        Returns: Alignment results

        """
        func_name = "%s_stats" % algorithm
        func = getattr(parasail, func_name)
        results = func(query, ref, 1, 1, Distance.dna_alignment_matrix)
        return results.score, results.end_ref


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

    def match(self, sequence, algorithm="hamming"):
        s1 = self.sequence.upper()
        s2 = sequence.upper()
        if hasattr(Distance, algorithm):
            func = getattr(Distance, algorithm)
            return func(s1, s2)[0]
        func = Distance.dna_alignment
        return func(s1, s2, algorithm)[0]
