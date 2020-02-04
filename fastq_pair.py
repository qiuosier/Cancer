import os
import logging
import math
import dnaio
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
        matched = [None] * len(reads)
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

    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

    def extract_reads_by_adapters(self, adapters, output_dir, error_rate=0.1):
        print("Adapters: %s" % adapters)
        counter = 0
        counter_matched = 0
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
                adapter = adapter1 or adapter2
                if adapter:
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
                    print("%s reads processed." % counter)
        print("%s reads processed." % counter)
        print("%s reads matched." % counter_matched)
        print("%s reads unmatched." % counter_unmatched)

    # @staticmethod
    # def __match_adapters(read1, read2, adapters, error_rate=0.1):
    #     for adapter in adapters:
    #         for read in [read1, read2]:
    #             # Align adapter to read 1
    #             result = parasail.sg_de_stats(adapter, read.sequence[:15], 1, 1, FASTQPair.score_matrix)
    #             if result.matches < 3:
    #                 continue
    #             min_score = - math.floor(result.end_ref * error_rate)
    #             if result.score >= min_score:
    #                 read.sequence = read.sequence[result.end_ref + 1:]
    #                 read.qualities = read.qualities[result.end_ref + 1:]
    #                 return adapter, read1, read2
    #     return None, read1, read2

    def build_index(self):
        fastq1_dict = dict()
        with dnaio.open(self.r1, file2=self.r2) as fastq1:
            for read1, read2 in fastq1:
                ident = read1.name.split(" ", 1)[0]
                fastq1_dict[ident] = (read1.sequence, read2.sequence)
        print("%d reads in FASTQ1" % len(fastq1_dict.keys()))
        return fastq1_dict

    def diff(self, r1, r2, output_dir):
        """Compares the reads with another pair of FASTQ file.

        Args:
            r1:
            r2:
            output_dir

        Returns:

        """
        fastq1_dict = self.build_index()
        fastq1_found = dict()
        # Stores reads found in FASTQ1 only
        diff_fastq1 = os.path.join(output_dir, "diff_fastq1.txt")
        # Stores reads found in FASTQ2 only
        diff_fastq2 = os.path.join(output_dir, "diff_fastq2.txt")
        # Stores reads found in both FASTQs but with difference sequence
        diff_seq = os.path.join(output_dir, "diff_seq.txt")
        counter_1 = 0
        counter_2 = 0
        counter_d = 0
        counter = 0
        with dnaio.open(r1, file2=r2) as fastq2, \
                open(diff_fastq2, 'w') as diff2, \
                open(diff_seq, 'w') as diffs:
            for read1, read2 in fastq2:
                counter += 1
                if counter % 100000 == 0:
                    print("%s reads processed." % counter)
                ident = read1.name.split(" ", 1)[0]
                if ident in fastq1_dict:
                    fastq1_found[ident] = True
                    f1_seq1, f1_seq2 = fastq1_dict.get(ident)
                    if self.__match_reads(read1.sequence, f1_seq1) and self.__match_reads(read2.sequence, f1_seq2):
                        continue
                    counter_d += 1
                    diffs.write(ident + '\n')
                    if not self.__match_reads(read1.sequence, f1_seq1):
                        diffs.write("F1R1: " + f1_seq1 + '\n')
                        diffs.write("F2R1: " + read1.sequence + '\n')
                    if not self.__match_reads(read2.sequence, f1_seq2):
                        diffs.write("F1R2: " + f1_seq2 + '\n')
                        diffs.write("F2R2: " + read2.sequence + '\n')
                else:
                    # Read not found in FASTQ1
                    counter_2 += 1
                    diff2.write(ident + '\n')
                    diff2.write("R1: " + read1.sequence + '\n')
                    diff2.write("R2: " + read2.sequence + '\n')
        # Reads in FASTQ1 only
        with open(diff_fastq1, 'w') as diff1:
            for ident in fastq1_dict.keys():
                if ident in fastq1_found:
                    continue
                counter_1 += 1
                read1, read2 = fastq1_dict.get(ident)
                diff1.write(ident + '\n')
                diff1.write("R1: " + read1 + '\n')
                diff1.write("R2: " + read2 + '\n')

        print("%d reads in FASTQ2." % counter)
        print("%d reads in FASTQ1 only." % counter_1)
        print("%d reads in FASTQ2 only." % counter_2)
        print("%d reads have difference sequence." % counter_d)

    @staticmethod
    def __match_reads(s1, s2):
        return s1.endswith(s2)
