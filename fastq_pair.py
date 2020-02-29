import os
import logging
import math
import dnaio
import gc
import parasail
logger = logging.getLogger(__name__)

match_score = 2
score_matrix = parasail.matrix_create("ACGT", match_score, -1)


class ReadPairs:
    def __init__(self, read1, read2):
        self.read1 = read1
        self.read2 = read2

    def trim(self, adapters, error_rate):
        # Trim both read1 and read2 with all adapters before return
        reads = [self.read1, self.read2]
        # Indicates whether R1 or R2 matches the adapter.
        matched = [""] * len(reads)
        for i in range(len(reads)):
            read = reads[i]
            for adapter in adapters:
                result = parasail.sg_de_stats(adapter, read.sequence[:20], 1, 1, score_matrix)
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


class FASTQPair:
    """Represents a pair of FASTQ files
    """

    # Output Filenames
    r1_matched_filename = "R1_matched.fastq.gz"
    r2_matched_filename = "R2_matched.fastq.gz"
    r1_unmatched_filename = "R1_unmatched.fastq.gz"
    r2_unmatched_filename = "R2_unmatched.fastq.gz"

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

    def extract_reads_by_adapters(self, adapters, output_dir, error_rate=0.2):
        """Extracts reads from FASTQ files by matching the adapters at the beginning of either R1 or R2.
        The reads are extracted into two groups, matched or unmatched.
        The matched adapters are trimmed from the reads.

        Args:
            adapters (list): A list of strings, each is a barcode adapter.
            output_dir: The output directory for storing the extracted reads.
            error_rate: Max percentage of error allowed.

        """
        print("Adapters: %s" % adapters)
        counter = 0
        counter_matched = 0
        counts = dict()
        counter_unmatched = 0
        r1_match_path = os.path.join(output_dir, "R1_matched.fastq.gz")
        r2_match_path = os.path.join(output_dir, "R2_matched.fastq.gz")
        r1_unmatch_path = os.path.join(output_dir, "R1_unmatched.fastq.gz")
        r2_unmatch_path = os.path.join(output_dir, "R2_unmatched.fastq.gz")

        with dnaio.open(self.r1, file2=self.r2) as fastq_in, \
                dnaio.open(r1_match_path, file2=r2_match_path, mode='w') as out_match, \
                dnaio.open(r1_unmatch_path, file2=r2_unmatch_path, mode='w') as out_unmatch:
            for read1, read2 in fastq_in:
                # adapter, trimmed_read1, trimmed_read2 = self.__match_adapters(read1, read2, adapters, error_rate)
                read_pair = ReadPairs(read1, read2)
                adapter1, adapter2 = read_pair.trim(adapters, error_rate)

                if adapter1:
                    key = "%s_1" % adapter1
                    c = counts.get(key, 0)
                    c += 1
                    counts[key] = c
                if adapter2:
                    key = "%s_2" % adapter2
                    c = counts.get(key, 0)
                    c += 1
                    counts[key] = c

                # The longer adapter has higher priority
                adapter = adapter1 if len(adapter1) > len(adapter2) else adapter2
                if adapter:
                    # Count the number of reads matching the longer adapter
                    c = counts.get(adapter, 0)
                    c += 1
                    counts[adapter] = c

                    trimmed_read1 = read_pair.read1
                    trimmed_read2 = read_pair.read2
                    # Sequence matched a barcode
                    counter_matched += 1
                    out_match.write(trimmed_read1, trimmed_read2)
                else:
                    # Sequence does not match a barcode
                    counter_unmatched += 1
                    out_unmatch.write(read1, read2)
                counter += 1
                if counter % 100000 == 0:
                    self.print_output("%s reads processed." % counter)
        self.print_output("%s reads processed." % counter)
        self.print_output("%s reads matched." % counter_matched)
        self.print_output("%s reads unmatched." % counter_unmatched)
        self.print_output("Output Files:\n%s" % "\n".join([
            r1_match_path, r2_match_path, r1_unmatch_path, r2_unmatch_path
        ]))
        counts['unmatched'] = counter_unmatched
        counts['total'] = counter
        return counts

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
                ident = read1.name.split(" ", 1)[0]
                fastq1_dict[ident] = (read1.sequence, read2.sequence)
                counter += 1
                if size and len(fastq1_dict.keys()) >= size:
                    print("%d reads indexed" % counter)
                    counter = 0
                    yield fastq1_dict
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
        the current file (self.r1 and self.r2) is referred as FASTQ1,
        the file we are comparing with (from r1 and r2 in the input args) is referred as FASTQ2

        Returns:

        """
        if not chunk_size:
            chunk_size = self.DEFAULT_CHUNK_SIZE

        fastq1_found = dict()

        # Output filenames
        # Stores reads found in FASTQ1 only
        fastq1_only = os.path.join(output_dir, "fastq1_only.txt")
        # Stores reads found in FASTQ2 only
        fastq2_only = os.path.join(output_dir, "fastq2_only.txt")

        # Stores reads found in both FASTQs but with difference sequence
        # Reads that are trimmed slightly differently
        diff_trim = os.path.join(output_dir, "diff_trim.txt")
        # Reads
        diff_seq = os.path.join(output_dir, "diff_seq.txt")

        counter_diff = 0
        # The number of reads that are exactly the same
        counter_same = 0
        # The number of reads that are not exactly the same but matched
        counter_match = 0

        # The number of reads in FASTQ1
        counter_1 = 0
        # The number of reads processed, also the number of reads in FASTQ2
        counter_2 = 0

        # The number of reads in FASTQ1 only
        counter_1_only = 0
        for fastq1_dict in self.build_index(chunk_size):
            counter_2 = 0
            counter_1 += len(fastq1_dict.keys())
            with dnaio.open(r1, file2=r2) as fastq2, \
                    open(fastq2_only, 'w') as f2_only, \
                    open(diff_trim, 'w') as f_trim, \
                    open(diff_seq, 'w') as f_diff:
                for read1, read2 in fastq2:
                    counter_2 += 1
                    if counter_2 % 100000 == 0:
                        print("%s reads processed." % counter_2)
                    ident = read1.name.split(" ", 1)[0]
                    f2_seq1 = read1.sequence
                    f2_seq2 = read2.sequence
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
                            f_trim.write("F2R1: " + f2_seq1 + '\n')
                            f_trim.write("F2R2: " + f2_seq2 + '\n')
                            counter_match += 1
                            continue

                        # Reads with different sequences
                        counter_diff += 1
                        f_diff.write(ident + '\n')
                        if not self.__match_trimmed_reads(f2_seq1, f1_seq1):
                            f_diff.write("F1R1: " + f1_seq1 + '\n')
                            f_diff.write("F2R1: " + f2_seq1 + '\n')
                        if not self.__match_trimmed_reads(f2_seq2, f1_seq2):
                            f_diff.write("F1R2: " + f1_seq2 + '\n')
                            f_diff.write("F2R2: " + f2_seq2 + '\n')
                    else:
                        # Read not found in FASTQ1
                        f2_only.write(ident + '\n')
                        f2_only.write("R1: " + read1.sequence + '\n')
                        f2_only.write("R2: " + read2.sequence + '\n')

            # Reads in FASTQ1 only
            with open(fastq1_only, 'w') as f1_only:
                for ident in fastq1_dict.keys():
                    if ident in fastq1_found:
                        continue
                    counter_1_only += 1
                    read1, read2 = fastq1_dict.get(ident)
                    f1_only.write(ident + '\n')
                    f1_only.write("R1: " + read1 + '\n')
                    f1_only.write("R2: " + read2 + '\n')

            # Clear the dictionary so that it will take the next chunk.
            fastq1_dict.clear()
            gc.collect()
        print("......")
        print("%d reads in FASTQ1." % counter_1)
        print("%d reads in FASTQ2." % counter_2)
        print("%d reads in FASTQ1 only." % counter_1_only)
        both = counter_same + counter_match + counter_diff
        print("%d reads in FASTQ2 only." % (counter_2 - both))
        print("%d reads in both pairs." % both)
        print("%d reads in both pairs are exactly the same." % counter_same)
        print("%d reads in FASTQ1 are a substrings of reads in FASTQ2." % counter_match)
        print("%d reads have difference sequence." % counter_diff)

    @staticmethod
    def __match_trimmed_reads(s1, s2):
        """Compare two reads and see if they matches each other.

        The two sequences are considered as matched if
            They are exactly the same, OR
            Read 1 ends with read 2 AND is mo more than 12bp longer than read 2

        Args:
            s1: Sequence of read 1
            s2: Sequence of read 2

        Returns:

        """
        return (s1 == s2) or (len(s1) - len(s2) <= 12 and s1.endswith(s2))
