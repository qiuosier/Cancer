import os
import re
import io
import gzip
import logging
from commons.Aries.storage import StorageFile
from commons.Aries.gcp.storage import GSFile
from commons.Aries.visual.plotly import PlotlyFigure
from commons.Aries.collections import sort_lists
from commons.Aries.tasks import ShellCommand
from .sequence import Sequence
logger = logging.getLogger(__name__)


class FASTQRead:
    """Represents a read in FASTQ file."""
    def __init__(self, lines):
        if len(lines) != 4:
            raise ValueError("lines must be a list of 4 strings.")
        self.identifier = lines[0]
        self.sequence = lines[1]
        self.description = lines[2]
        self.quality = lines[3]


class FASTQGzip:
    def __init__(self, uri):
        self.uri = uri
        self.gzip = gzip.GzipFile(fileobj=StorageFile.init(uri, "rb"))
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        read = FASTQRead(
            [self.gzip.readline(),
             self.gzip.readline(),
             self.gzip.readline(),
             self.gzip.readline()]
        )
        if read.identifier:
            self.current += 1
            return read
        else:
            raise StopIteration

    @property
    def read_count(self):
        logger.debug("Counting reads in file %s..." % self.uri)
        self.gzip = gzip.GzipFile(fileobj=StorageFile.init(self.uri, "rb").local())
        return len(list(self))


class IlluminaFASTQ:
    dual_index_pattern = r"[ACGTN]{8}\+[ACGTN]{8}"

    def __init__(self, file_path):
        file_path = str(file_path)
        if file_path.startswith("gs://"):
            if not GSFile(file_path).blob.exists():
                raise FileNotFoundError("File not found at %s." % file_path)
        elif not os.path.exists(file_path):
            raise FileNotFoundError("File not found at %s." % file_path)
        self.file_path = file_path

    def peek_barcode(self):
        if self.file_path.startswith("gs://"):
            cmd = "gsutil cat %s | zcat | head -n 4000" % self.file_path
            job = ShellCommand(cmd)
            job.run()
            f = io.StringIO(job.std_out)
        else:
            f = open(self.file_path, 'r')
        barcode_dict = self.__process_barcode(f, self.__count_barcode)
        f.close()
        return barcode_dict

    def __process_barcode(self, lines, method):
        barcode_dict = {}
        for i, line in enumerate(lines, start=1):
            if not line.startswith("@"):
                continue
            barcode = line.strip().split(":")[-1]
            if re.match(self.dual_index_pattern, barcode):
                idx = barcode.split("+")
                i7 = idx[0]
                i5 = Sequence(idx[1]).reverse_complements
                barcode = "%s+%s" % (i7, i5)
                barcode_dict[barcode] = method(barcode_dict, barcode, i)
        return barcode_dict

    @staticmethod
    def __count_barcode(barcode_dict, barcode, row_number):
        return barcode_dict.get(barcode, 0) + 1

    @staticmethod
    def __group_barcode(barcode_dict, barcode, row_number):
        line_list = barcode_dict.get(barcode, [])
        line_list.append(row_number)
        return line_list

    def group_by_barcode(self, threshold=0):
        with open(self.file_path, 'r') as f:
            barcode_dict = self.__process_barcode(f, self.__group_barcode)
        if threshold > 0:
            barcode_dict = {k: v for k, v in barcode_dict.items() if len(v) > threshold}
        return barcode_dict

    def count_by_barcode(self, threshold=0):
        with open(self.file_path, 'r') as f:
            barcode_dict = self.__process_barcode(f, self.__count_barcode)
        return {k: len(v) for k, v in barcode_dict.items() if v > threshold}


class BarcodeStatistics:
    def __init__(self, barcode_dict):
        self.barcode_dict = barcode_dict

    def filter_by_reads(self, threshold=0):
        self.barcode_dict = {k: v for k, v in self.barcode_dict.items() if v > threshold}
        return self

    def sort_data(self, max_size=0):
        labels = []
        counts = []
        for k, v in self.barcode_dict.items():
            labels.append(k)
            counts.append(v)
        counts, labels = sort_lists(counts, labels, reverse=True)
        if max_size and len(counts) > max_size:
            counts = counts[:max_size]
            labels = labels[:max_size]
        return counts, labels

    def histogram(self, max_bins=20):
        counts, labels = self.sort_data(max_bins)
        return PlotlyFigure().add_trace("Histogram", x=labels, y=counts, histfunc='sum')

    def bar_chart(self, max_size=20):
        counts, labels = self.sort_data(max_size)
        counts.reverse()
        labels.reverse()
        return PlotlyFigure(height=1000, font=dict(
            family='Courier New, monospace'
        )).bar(x=counts, y=labels, orientation='h')
