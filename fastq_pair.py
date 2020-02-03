import gzip
import os
import logging
import math
import dnaio
from .sequence import Sequence
logger = logging.getLogger(__name__)


class FASTQPair:
    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

    def extract_reads_by_adapters(self, adapters, output_dir):
        print("Adapters: %s" % adapters)
        # adapters is a dictionary of:
        # key: the barcode
        # value: the number of mismatches allowed
        adapters = {Sequence(s): math.floor(len(s)) * 0.2 for s in adapters}
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
                adapter, trimed_read1, trimed_read2 = self.__match_adapters(read1, read2, adapters)
                if adapter:
                    # Sequence matched a barcode
                    counter_matched += 1
                    out_match.write(trimed_read1, trimed_read2)
                else:
                    # Sequence does not match a barcode
                    counter_unmatched += 1
                    out_unmatch.write(trimed_read1, trimed_read2)
                counter += 1
                if counter % 50000 == 0:
                    print("%s reads processed." % counter)
            print("%s reads processed." % counter)
            print("%s reads matched." % counter_matched)
            print("%s reads unmatched." % counter_unmatched)

    def __match_adapters(self, read1, read2, adapters):
        # if not read1.identifier and not read2.identifier:
        #     raise StopIteration
        # if not (read1.identifier and read2.identifier):
        #     logger.error("R1 and R2 does not have the same number of reads:")
        #     logger.error(read1.identifier)
        #     logger.error(read2.identifier)
        #     raise StopIteration

        for adapter, max_mismatch in adapters.items():
            if adapter.match(read1.sequence[:adapter.length]) <= max_mismatch:
                read1.sequence = read1.sequence[adapter.length:]
                return adapter, read1, read2
            if adapter.match(read2.sequence[:adapter.length]) <= max_mismatch:
                read2.sequence = read2.sequence[adapter.length:]
                return adapter, read1, read2
        return None, read1, read2
