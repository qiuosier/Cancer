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


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence.upper()

    @property
    def complements(self):
        return "".join([COMPLEMENTS[c] for c in self.sequence])

    @property
    def reverse_complements(self):
        return self.complements[::-1]

    def match(self, sequence):
        s1 = self.sequence.upper()
        s2 = sequence.upper()
        l1 = len(s1)
        l2 = len(s2)
        mismatch = 0
        for i in range(max(l1, l2)):
            c_1 = s1[i] if i < l1 else ""
            c_2 = s2[i] if i < l2 else ""
            if c_1 != c_2:
                mismatch += 1
        return mismatch
