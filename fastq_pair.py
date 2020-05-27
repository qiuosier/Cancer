import os
import logging
import math
import dnaio
import gc
import parasail
logger = logging.getLogger(__name__)

match_score = 0
DEFAULT_SCORE_MATRIX = parasail.matrix_create("ACGTN", match_score, -1)


class ReadPair:
    def __init__(self, read1, read2):
        self.read1 = read1
        self.read2 = read2

        arr1 = self.read1.name.split(" ", 1)
        arr2 = self.read2.name.split(" ", 1)

        if arr1[0] != arr2[0]:
            raise ValueError(
                "Identifiers of the read pairs does not match each other.\n"
                "Read1: %s\n"
                "Read2: %s" % (arr1[0], arr2[0])
            )

        self._identifier = arr1[0]

        if len(arr1) > 1 and len(arr2) > 1:
            pair1 = arr1[1][0]
            pair2 = arr2[1][0]
            if pair1 == pair2:
                raise ValueError("Both reads are pair %s" % pair1)

            # Check if read1 and read2 are switched
            if str(pair1) == '2' and str(pair2) == '1':
                tmp = self.read1
                self.read1 = self.read2
                self.read2 = tmp
        else:
            raise ValueError("Invalid read identifiers.\nRead1: %s\nRead2: %s" % (arr1[0], arr2[0]))

    @property
    def identifier(self):
        return self._identifier

    @property
    def barcode(self):
        """Barcode as it appear in the read name

        This is the plain barcode taken from the read name.
        For Illumina sequencer, this is I7 + reverse compliment of I5.

        """
        barcode1 = self.read1.name.strip().split(":")[-1]
        barcode2 = self.read2.name.strip().split(":")[-1]
        if not barcode1 == barcode2:
            raise ValueError("Read1 and Read2 have different barcodes.\nRead1: %s\nRead2: %s" % (barcode1, barcode2))
        return barcode1

    @property
    def reads(self):
        return self.read1, self.read2

    def trim(self, adapters, error_rate):
        # Trim both read1 and read2 with all adapters before return
        reads = [self.read1, self.read2]
        # Indicates whether R1 or R2 matches the adapter.
        matched = [""] * len(reads)
        for i in range(len(reads)):
            read = reads[i]
            for adapter in adapters:
                result = parasail.sg_de_stats(adapter, read.sequence[:20], 1, 1, DEFAULT_SCORE_MATRIX)
                if result.matches <= 5:
                    continue
                distance = match_score * result.matches - result.score
                max_distance = math.floor(len(adapter) * error_rate)
                if distance <= max_distance:
                    matched[i] = adapter
                    read.sequence = read.sequence[result.end_ref + 1:]
                    read.qualities = read.qualities[result.end_ref + 1:]
                    break
        return matched[0], matched[1]

    def check_identifier(self):
        ident1 = self.read1.name.split(" ", 1)[0]
        ident2 = self.read2.name.split(" ", 1)[0]
        if ident1 != ident2:
            raise ValueError(
                "Identifiers of the read pairs does not match each other.\n"
                "Read1: %s\n"
                "Read2: %s" % (ident1, ident2)
            )
        return ident1


class FASTQPair:
    """Represents a pair of FASTQ files
    """

    def __init__(self, r1, r2, ident=None):
        """Initializes an object representing a pair of FASTQ files

        Args:
            r1: FASTQ R1
            r2: FASTQ R2
            ident: An optional identifier for the FASTQ pair.
        """
        self.r1 = r1
        self.r2 = r2
        self.ident = ident

    def print_output(self, s):
        """Prints a message with processor ID
        This is used to identify the message source when using multi-processing.
        """
        if self.ident:
            print("%s: %s" % (self.ident, s))
        else:
            print(s)

    def build_index(self, size=None):
        """Builds an dictionary to index the reads.

        Returns (dict):
            If size is None, A dictionary, where,
            key: An identifier of the read pair, this is the first part of the identifier line up to the first space.
            value: A 2-tuple of read sequence.
        """
        fastq1_dict = dict()
        counter = 0
        with dnaio.open(self.r1, file2=self.r2) as fastq1:
            for read1, read2 in fastq1:
                read_pair = ReadPair(read1, read2)
                ident = read_pair.identifier
                fastq1_dict[ident] = (read_pair.read1.sequence, read_pair.read2.sequence)
                counter += 1
                if size and len(fastq1_dict.keys()) >= size:
                    print("%d reads indexed" % counter)
                    counter = 0
                    yield fastq1_dict
                    # Clear the dictionary so that it will take the next chunk.
                    fastq1_dict.clear()
                    gc.collect()
        print("%d reads indexed" % counter)
        yield fastq1_dict

    # Max index size: 1M keys
    DEFAULT_CHUNK_SIZE = 1 * 1000 * 1000

    def diff(self, r1, r2, output_dir, chunk_size=DEFAULT_CHUNK_SIZE):
        """Compares the reads with another pair of FASTQ file.

        Args:
            r1:
            r2:
            output_dir:
            chunk_size:

        In this method,
        the current file (self.r1 and self.r2) is referred as P1,
        the file we are comparing with (from r1 and r2 in the input args) is referred as P2

        Returns:

        """
        if not chunk_size:
            chunk_size = self.DEFAULT_CHUNK_SIZE

        fastq1_found = dict()

        # Output filenames
        # Stores reads found in FASTQ1 only
        fastq1_only = os.path.join(output_dir, "P1_only.txt")
        # Stores reads found in FASTQ2 only
        fastq2_only = os.path.join(output_dir, "P2_only.txt")

        # Stores reads found in both FASTQs but with difference sequence
        # Reads that are trimmed slightly differently
        diff_trim = os.path.join(output_dir, "diff_trim.txt")
        # Reads
        diff_seq = os.path.join(output_dir, "diff_seq.txt")

        # Remove existing files.
        for file_path in [fastq1_only, fastq2_only, diff_trim, diff_seq]:
            if os.path.exists(file_path):
                os.remove(file_path)

        counter_diff = 0
        # The number of reads that are exactly the same
        counter_same = 0
        # The number of reads that are not exactly the same but matched
        counter_match = 0

        # The number of reads in P1
        counter_1 = 0
        # The number of reads in P2
        counter_2 = 0

        # The number of reads in P1 only
        counter_1_only = 0
        for fastq1_dict in self.build_index(chunk_size):
            counter_2 = 0
            counter_1 += len(fastq1_dict.keys())
            with dnaio.open(r1, file2=r2) as fastq2, \
                    open(diff_trim, 'a') as f_trim, \
                    open(diff_seq, 'a') as f_diff:
                for read1, read2 in fastq2:
                    counter_2 += 1
                    if counter_2 % 100000 == 0:
                        print("%s reads in P2 processed." % counter_2)
                    read_pair = ReadPair(read1, read2)
                    ident = read_pair.identifier
                    f2_seq1 = read_pair.read1.sequence
                    f2_seq2 = read_pair.read2.sequence
                    if ident in fastq1_dict:
                        fastq1_found[ident] = True
                        f1_seq1, f1_seq2 = fastq1_dict.get(ident)

                        if f2_seq1 == f1_seq1 and f2_seq2 == f1_seq2:
                            # Both reads are exactly the same
                            counter_same += 1
                            continue
                        elif self.__match_trimmed_reads(f2_seq1, f1_seq1) and \
                                self.__match_trimmed_reads(f2_seq2, f1_seq2):
                            # Both reads matches but FASTQ 1 is trimmed more.
                            f_trim.write(ident + '\n')
                            f_trim.write("P2R1: " + f2_seq1 + '\n')
                            f_trim.write("P2R2: " + f2_seq2 + '\n')
                            counter_match += 1
                            continue

                        # Reads with different sequences
                        counter_diff += 1
                        f_diff.write(ident + '\n')
                        if not self.__match_trimmed_reads(f2_seq1, f1_seq1):
                            f_diff.write("P1R1: " + f1_seq1 + '\n')
                            f_diff.write("P2R1: " + f2_seq1 + '\n')
                        if not self.__match_trimmed_reads(f2_seq2, f1_seq2):
                            f_diff.write("P1R2: " + f1_seq2 + '\n')
                            f_diff.write("P2R2: " + f2_seq2 + '\n')
                    else:
                        pass
                        # TODO: Read not found in FASTQ1
                        # TODO: The following does not work if there are multiple chunks
                        # f2_only.write(ident + '\n')
                        # f2_only.write("R1: " + read1.sequence + '\n')
                        # f2_only.write("R2: " + read2.sequence + '\n')

            # Reads in FASTQ1 only
            with open(fastq1_only, 'w') as f1_only:
                for ident in fastq1_dict.keys():
                    if ident in fastq1_found:
                        continue
                    counter_1_only += 1
                    seq1, seq2 = fastq1_dict.get(ident)
                    f1_only.write(ident + '\n')
                    f1_only.write("R1: " + seq1 + '\n')
                    f1_only.write("R2: " + seq2 + '\n')

        print("=" * 50)
        print("Pair 1: %s and %s" % (self.r1, self.r2))
        print("Pair 2: %s and %s" % (r1, r2))
        print("=" * 50)
        print("%d reads in Pair 1." % counter_1)
        print("%d reads in Pair 2." % counter_2)
        print("%d reads in Pair 1 only." % counter_1_only)
        both = counter_same + counter_match + counter_diff
        print("%d reads in Pair 2 only." % (counter_2 - both))
        print("%d reads in both pairs." % both)
        print("%d reads in both pairs are exactly the same." % counter_same)
        print("%d reads in Pair 1 are trimmed more." % counter_match)
        print("%d reads have difference sequence or trimmed differently." % counter_diff)

    @staticmethod
    def __match_trimmed_reads(s1, s2):
        """Compare two read sequences and see if s1 contains s2.

        Args:
            s1: Sequence of read 1
            s2: Sequence of read 2

        Returns:

        """
        return s2 in s1
