import logging
from .processor import FASTQProcessor, FASTQWorker
logger = logging.getLogger(__name__)


class BarcodeWorker(FASTQWorker):

    def __init__(self, start_pos=0, length=10):
        super().__init__()
        self.start_pos = start_pos
        self.length = length

    def process_read_pair(self, read):
        if isinstance(read, tuple):
            read1, read2 = read
            self.add_count("%s_1" % read1.sequence[self.start_pos:self.start_pos + self.length])
            self.add_count("%s_2" % read2.sequence[self.start_pos:self.start_pos + self.length])
        else:
            self.add_count("%s" % read.sequence[self.start_pos:self.start_pos + self.length])


class BarcodeCounter(FASTQProcessor):

    CUMULATIVE_PERCENT_CUTOFF = 80
    SINGLE_PERCENT_CUTOFF = 0.05

    def __init__(self, start=0, length=10):
        """Processor for getting the usage count of inline barcodes in FASTQ files.
        
        Args:
            start: Starting position of the barcode (0-based)
            length: Length of the barcode
        """
        super().__init__(BarcodeWorker, start_pos=start, length=length)

    def collect_results(self, jobs):
        # Collect the statistics.
        results = [job.get() for job in jobs]
        for r in results:
            self.update_counts(r)
        count_list = [(k, v) for k, v in self.counts.items()]
        count_list.sort(key=lambda x: x[1], reverse=True)
        total = sum(self.counts.values())
        print("Total Reads: %s" % total)
        cutoff = total * self.CUMULATIVE_PERCENT_CUTOFF / 100
        cumulative_sum = 0
        for t in count_list:
            if t[0] in ['total', 'matched', 'unmatched']:
                continue
            percent = t[1] / total * 100
            print("%s: %10s, %3.2f%%" % (t[0], t[1], percent))
            cumulative_sum += t[1]
            if cumulative_sum > cutoff or percent < self.SINGLE_PERCENT_CUTOFF:
                break
        return self.counts
