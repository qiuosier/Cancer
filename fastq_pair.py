import os
import logging
import math
import dnaio
import parasail
logger = logging.getLogger(__name__)


class FASTQPair:
    alignment_matrix = parasail.matrix_create("ACGT", 0, -1)

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
                adapter, trimmed_read1, trimmed_read2 = self.__match_adapters(read1, read2, adapters, error_rate)
                if adapter:
                    # Sequence matched a barcode
                    counter_matched += 1
                    out_match.write(trimmed_read1, trimmed_read2)
                else:
                    # Sequence does not match a barcode
                    counter_unmatched += 1
                    out_unmatch.write(trimmed_read1, trimmed_read2)
                counter += 1
                if counter % 100000 == 0:
                    print("%s reads processed." % counter)
        print("%s reads processed." % counter)
        print("%s reads matched." % counter_matched)
        print("%s reads unmatched." % counter_unmatched)

    @staticmethod
    def __match_adapters(read1, read2, adapters, error_rate=0.1):
        for adapter in adapters:
            for read in [read1, read2]:
                # Align adapter to read 1
                result = parasail.sg_de_stats(adapter, read.sequence[:15], 1, 1, FASTQPair.alignment_matrix)
                if result.matches < 3:
                    continue
                min_score = - math.floor(result.matches * error_rate)
                if result.score >= min_score:
                    read.sequence = read.sequence[result.end_ref + 1:]
                    read.qualities = read.qualities[result.end_ref + 1:]
                    return adapter, read1, read2
        return None, read1, read2

    # @staticmethod
    # def trim_read(read, adapter):
    #     k = adapter.length - 1
    #     while read.sequence[k] == adapter.sequence[k] and k > 0:
    #         k -= 1
    #     k += 1
    #     read.sequence = read.sequence[k:]
    #     read.qualities = read.qualities[k:]
    #     return read
