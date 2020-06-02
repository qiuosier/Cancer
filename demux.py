import multiprocessing
import os
import csv
import math
import dnaio
import parasail
import editdistance
import re
import logging
import tempfile
import datetime
import time
from multiprocessing import Process, Manager
from .fastq_pair import ReadPair
from .fastq import IlluminaFASTQ, BarcodeStatistics
logger = logging.getLogger(__name__)


class DemultiplexWriter(dict):
    """A dictionary subclass holding barcode adapters and
    the corresponding file-like objects for writing fastq.gz files.

    In the dictionary:
    Each key is a barcode.
    The actual filenames are specified by the paired_end_filenames() method.
    Each value is a file-like object returned by opening a pair of fastq.gz files.

    Attributes:
        barcode_dict: A dictionary mapping barcode to filename prefix.
        prefix_dict: A dictionary mapping filename prefix to file-like object.
            prefix_dict can be used to determine the files with certain prefix are opened.

    This class supports context manager, for example:
        with DemultiplexWriter(barcode_dict) as writer:
            ...PROCESSING CODE HERE...

    """
    @staticmethod
    def paired_end_filenames(prefix):
        """Maps a prefix to a 2-tuple of filenames (R1, R2)

        Args:
            prefix (str): Prefix for the filenames, including the full path.

        Returns: A 2-tuple of strings as the filenames for R1 and R2 FASTQ files.

        """
        return prefix + ".R1.fastq.gz", prefix + ".R2.fastq.gz"

    def __init__(self, barcode_dict):
        """Initializes the writer with a dictionary mapping barcode to filename prefix.

        Args:
            barcode_dict: A dictionary mapping barcode to filename prefix
                Each key is a barcode.
                Each value is a prefix for output filename, including the full path.
                The output file will contain the reads corresponds to the barcode.
                If multiple barcodes are mapping to the same prefix,
                    the reads with those barcodes will be written into the same output file pair.

        """
        self.barcode_dict = barcode_dict
        self.prefix_dict = {}
        super().__init__()

    def open(self):
        """Opens the files for writing
        """
        for barcode, prefix in self.barcode_dict.items():
            if not prefix:
                self[barcode] = None
            if prefix in self.prefix_dict.keys():
                self[barcode] = self.prefix_dict[prefix]
            else:
                r1_out, r2_out = DemultiplexWriter.paired_end_filenames(prefix)
                fp = dnaio.open(r1_out, file2=r2_out, mode='w')
                self.prefix_dict[prefix] = fp
                self[barcode] = fp
        return self

    def close(self):
        """Closes the files
        """
        for fp in self.values():
            fp.close()

    def write(self, barcode, read1, read2):
        fp = self.get(barcode)
        if not fp:
            return
        fp.write(read1, read2)

    def __enter__(self):
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self.close()


class DemultiplexWorker:

    DEFAULT_ERROR_RATE = 0.1

    def __init__(self, barcode_dict, error_rate=None, score=1, penalty=10):
        self.barcode_dict = barcode_dict
        self.adapters = list(barcode_dict.keys())
        self.min_match_length = round(min([len(adapter) / 2 for adapter in self.adapters]))

        # Set the default values
        self.error_rate = error_rate if error_rate else self.DEFAULT_ERROR_RATE
        self.score = int(score) if str(score).isdigit() else 1
        self.penalty = int(penalty) if str(penalty).isdigit() else 10
        logger.debug("Process %s, Penalty: %s, Error Rate: %s, Score: %s" % (
            os.getpid(), self.penalty, self.error_rate, self.score
        ))

        self.score_matrix = self.create_score_matrix()

        # Stores the statistics
        # Each sub-class may have different keys for the counts dictionary.
        # The values should be integers.
        # Whether to use this in the sub-class is optional.
        self.counts = dict(matched=0, unmatched=0, total=0)

    def create_score_matrix(self):
        return parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)

    def semi_global_distance(self, s1, s2):
        score_matrix = self.create_score_matrix()
        result = parasail.sg_de_stats(
            s1, s2, self.penalty, self.penalty, score_matrix
        )
        return (self.score * result.matches - result.score) / self.penalty

    def update_counts(self, counts):
        """Updates the demultiplex statistics.

        Args:
            counts: A dictionary storing the statistics.

        This method merges counts from the input argument with self.counts
        Values of the existing keys are added together.
        Keys from counts will be created if they are not already in self.counts

        """
        for k, v in counts.items():
            self.counts[k] = self.counts.get(k, 0) + v
        return self

    def add_count(self, key, val=1):
        """Increments the value of a particular key in counts (dictionary)
        This is a static method and it does NOT modify the self.counts
        """
        c = self.counts.get(key, 0)
        c += val
        self.counts[key] = c
        return self.counts

    def start(self, in_queue, out_queue):
        active_time = datetime.timedelta()
        batch_count = 0
        with DemultiplexWriter(self.barcode_dict) as writer:
            while True:
                reads = in_queue.get()
                # Keep the starting time for each batch processing
                timer_started = datetime.datetime.now()
                if reads is None:
                    logger.debug("Process %s, Active time: %s (%s batches, %s/batch)." % (
                        os.getpid(), active_time, batch_count, active_time / batch_count
                    ))
                    return self.counts
                # results = []
                for read_pair in reads:
                    barcode, read1, read2 = self.process_read_pair(read_pair)
                    writer.write(barcode, read1, read2)
                    # results.append(result)

                self.add_count('total', len(reads))
                batch_count += 1
                # Add processing time for this batch
                active_time += (datetime.datetime.now() - timer_started)
                out_queue.put(len(reads))

    def process_read_pair(self, read_pair):
        raise NotImplementedError


class DemultiplexInlineWorker(DemultiplexWorker):

    DEFAULT_ERROR_RATE = 0.2

    def trim_adapters(self, read1, read2):
        """Checks if the beginning of the reads in a read pair matches any adapter.
        If so, trim the reads to remove the matching adapter.

        A read and an adapter are matched by using semi-global alignment without penalty at the end of the read.
        They are considered as MATCHED if the number of substitutions and gaps

        Args:
            read1: Forward read.
            read2: Reverse Compliment read.

        read1 and read2 are dnaio Sequence objects.
        They are passed into this method as references.
        The modifications on read1 and read2 will be preserved after return.

        The score_matrix will be generated if not specified.
        However, if specified, it must be generated by:
            score_matrix = parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)
        Pass a score_matrix avoid this method to generate it every time, which may speed up the overall processing.

        Returns: A 2-tuple indicating whether any adapters are matching the read pair.
            If a read matched an adapter, the matching adapter will be returned in the tuple.
            Otherwise, the corresponding element in the tuple will be None.
            For example, ("ACTGACT", None) indicates barcode "ACTGACT" is matching read1 (forward read).

        See Also:
             https://github.com/marcelm/dnaio/blob/master/src/dnaio/_core.pyx
             https://github.com/jeffdaily/parasail-python#substitution-matrices

        """
        # Trim both read1 and read2 with all adapters before return
        reads = [read1, read2]
        # Indicates whether R1 or R2 matches the adapter.
        matched = [""] * len(reads)
        for i in range(len(reads)):
            read = reads[i]
            for adapter in self.adapters:
                result = parasail.sg_de_stats(
                    adapter, read.sequence[:20], self.penalty, self.penalty, self.score_matrix
                )
                if result.matches <= self.min_match_length:
                    continue
                distance = (self.score * result.matches - result.score) / self.penalty
                max_distance = math.floor(len(adapter) * self.error_rate)
                if distance <= max_distance:
                    matched[i] = adapter
                    read.sequence = read.sequence[result.end_ref + 1:]
                    read.qualities = read.qualities[result.end_ref + 1:]
                    break
        # read1 and read2 are preserved implicitly
        return matched[0], matched[1]

    def process_read_pair(self, read_pair):
        read1, read2 = read_pair
        # Initialize ReadPair to check if read1 and read2 are valid
        read1, read2 = ReadPair(read1, read2).reads
        # read1 and read2 are references
        # The modifications on read1 and read2 will be returned implicitly
        adapter1, adapter2 = self.trim_adapters(read1, read2)

        if adapter1:
            self.add_count("%s_1" % adapter1)
        if adapter2:
            self.add_count("%s_2" % adapter2)

        # The longer adapter has higher priority
        adapter = adapter1 if len(adapter1) > len(adapter2) else adapter2
        if adapter:
            # Count the number of reads matching the longer adapter
            self.add_count(adapter)
            # Sequence matched a barcode
            self.add_count('matched')

        else:
            # Sequence does not match a barcode
            self.add_count('unmatched')
            adapter = "NO_MATCH"
        return adapter, read1, read2


class DemultiplexProcess:
    def __init__(self, worker_class):
        self.worker_class = worker_class

        self.pool_size = max(os.cpu_count(), 1)
        logger.debug("Pool Size: %s" % self.pool_size)
        self.manager = Manager()
        self.pool = []

        self.adapters = []
        self.counts = dict()

        self.reader_queue = self.manager.Queue(self.pool_size * 100)
        self.worker_queue = self.manager.Queue(self.pool_size * 100)

    def update_counts(self, counts):
        for k, v in counts.items():
            self.counts[k] = self.counts.get(k, 0) + v
        return self.counts

    def start_readers(self, fastq_files):
        readers = []
        for fastq_pair in fastq_files:
            reader = Process(
                target=DemultiplexProcess.read_data,
                args=([fastq_pair], self.reader_queue, self.pool_size)
            )
            reader.start()
            readers.append(reader)
        return readers

    def finalize_readers(self, readers, reader_q_size):
        for reader in readers:
            reader.join()
        logger.debug("Reader's queue max size: %s" % max(reader_q_size))
        self.print_queue_size(reader_q_size)

    def wait_for_jobs(self, jobs):
        # Wait for the jobs to finish and keep track of the queue size
        reader_q_size = []
        # writer_q_size = []
        counter = 0
        while True:
            reader_q_size.append(self.reader_queue.qsize())
            ready = True
            for job in jobs:
                if not job.ready():
                    ready = False
            while not self.worker_queue.empty():
                counter += self.worker_queue.get()
            print("{:,} reads processed.".format(counter))
            if ready:
                break
            time.sleep(5)
        return reader_q_size

    def collect_results(self, jobs):
        # Collect the statistics.
        results = [job.get() for job in jobs]
        for r in results:
            self.update_counts(r)
        logger.debug(self.counts)
        return self.counts

    @staticmethod
    def prepare_concatenation(barcode_dict, output_list):
        # prefix_dict is a dict storing the final output prefix as keys, and
        # each value is a list of file prefixes (paths) to be concatenated.
        prefix_dict = {}
        for barcode, file_prefix in barcode_dict.items():
            # The following code will merge the barcodes pointing to the same file path
            path_list = prefix_dict.get(file_prefix, [])
            path_list.extend([
                output_dict.get(barcode)
                for output_dict in output_list
                if output_dict.get(barcode)
            ])
            prefix_dict[file_prefix] = path_list
        prefix_dict = {k: list(set(v)) for k, v in prefix_dict.items()}
        return prefix_dict

    def start(self, fastq_files, barcode_dict, error_rate=None, score=1, penalty=10):
        self.adapters.extend(barcode_dict.keys())
        pool = multiprocessing.Pool(self.pool_size)

        output_list = []
        jobs = []
        with tempfile.TemporaryDirectory() as temp_dir:
            for i in range(self.pool_size):
                ident = "Process_%s" % i

                output_dict = {k: os.path.join(temp_dir, "%s_%s" % (k, ident)) for k, v in barcode_dict.items()}
                output_list.append(output_dict)

                job = pool.apply_async(
                    self.start_worker,
                    (self.worker_class, self.reader_queue, self.worker_queue, output_dict, error_rate, score, penalty)
                )
                jobs.append(job)

            readers = self.start_readers(fastq_files)

            # Wait for the jobs to finish and keep track of the queue size
            reader_q_size = self.wait_for_jobs(jobs)

            # Collect the statistics.
            self.collect_results(jobs)
            self.finalize_readers(readers, reader_q_size)

            prefix_dict = self.prepare_concatenation(barcode_dict, output_list)
            self.concatenate_fastq(prefix_dict)

        return self

    @staticmethod
    def concatenate_fastq(prefix_dict):
        """Concatenates the FASTQ files in a list of directories.
        """
        for prefix, prefix_list in prefix_dict.items():
            logger.debug("Concatenating %s pairs of files." % len(prefix_list))
            r1, r2 = DemultiplexWriter.paired_end_filenames(prefix)
            pair_list = [DemultiplexWriter.paired_end_filenames(p) for p in prefix_list]

            cmd = "cat %s > %s" % (" ".join(p[0] for p in pair_list), r1)
            os.system(cmd)
            logger.debug(r1)
            cmd = "cat %s > %s" % (" ".join(p[1] for p in pair_list), r2)
            os.system(cmd)
            logger.debug(r2)

    @staticmethod
    def print_queue_size(q_size_array):
        for s in q_size_array:
            print("*" * s + str(s))

    @staticmethod
    def start_worker(worker_class, in_queue, out_queue, *args, **kwargs):
        return worker_class(*args, **kwargs).start(in_queue, out_queue)

    BATCH_SIZE = 5000

    @staticmethod
    def read_data(fastq_files, queue, pool_size):
        enqueue_time = datetime.timedelta()
        batch_count = 0

        process_started = datetime.datetime.now()
        for fastq_pair in fastq_files:
            with dnaio.open(fastq_pair[0], file2=fastq_pair[1]) as fastq_in:
                size = 0
                reads = []
                for read1, read2 in fastq_in:
                    reads.append((read1, read2))
                    size += 1
                    if size > DemultiplexProcess.BATCH_SIZE:
                        batch_count += 1

                        timer_started = datetime.datetime.now()
                        queue.put(reads)
                        enqueue_time += (datetime.datetime.now() - timer_started)

                        size = 0
                        reads = []

                timer_started = datetime.datetime.now()
                queue.put(reads)
                enqueue_time += (datetime.datetime.now() - timer_started)

        print('Finished reading files. Total time: %s, Enqueue time: %s, Batch size: %s, Batch count: %s' % (
            datetime.datetime.now() - process_started, enqueue_time, DemultiplexProcess.BATCH_SIZE, batch_count
        ))
        for i in range(pool_size):
            queue.put(None)

    @staticmethod
    def write_data(queue, barcode_dict):
        process_started = datetime.datetime.now()
        dequeue_time = datetime.timedelta()
        counter = 0
        with DemultiplexWriter(barcode_dict) as fp_dict:
            while True:
                timer_started = datetime.datetime.now()
                results = queue.get()
                dequeue_time += (datetime.datetime.now() - timer_started)
                if results is None:
                    print('Finished writing files. Total time: %s, Dequeue time: %s' % (
                        datetime.datetime.now() - process_started, dequeue_time
                    ))
                    return
                for barcode, read1, read2 in results:
                    counter += 1
                    fp = fp_dict.get(barcode)
                    if not fp:
                        continue
                    fp.write(read1, read2)

                    if counter % 100000 == 0:
                        print("{:,} reads processed.".format(counter))

    @staticmethod
    def parse_barcode_outputs(barcode_outputs):
        """Parses the barcode and output file prefix pairs specified as a list of strings like:
            "BARCODE=PREFIX", or "BARCODE_1 BARCODE_2=PREFIX"

        Args:
            barcode_outputs (list): A list of strings in the format of BARCODE=PREFIX

        Returns: A dictionary, where each key is a barcode, each value is the file_prefix.

        """
        barcode_dict = {}
        for output in barcode_outputs:
            arr = str(output).split("=", 1)
            barcode_list = arr[0].strip().split()
            file_prefix = arr[1] if len(arr) > 1 else None
            for barcode in barcode_list:
                barcode_dict[barcode] = file_prefix
        return barcode_dict

    @staticmethod
    def pair_fastq_files(r1_list, r2_list):
        """Pairs the FASTQ files from two lists by the identifier of the first read
        """
        if not r1_list or not r2_list:
            raise ValueError("Invalid R1 or R2.\nR1: %s\nR2: %s" % (r1_list, r2_list))
        # Convert r1_list and r2_list to lists if needed.
        if not isinstance(r1_list, list):
            r1_list = [r1_list]
        if not isinstance(r2_list, list):
            r2_list = [r2_list]

        r1_dict = {DemultiplexInline.get_identifier(f): f for f in r1_list}
        r2_dict = {DemultiplexInline.get_identifier(f): f for f in r2_list}

        if r1_dict.keys() != r2_dict.keys():
            raise ValueError("Unable to pair R1 and R2.\nR1: %s\nR2: %s" % (r1_dict.keys(), r2_dict.keys()))

        # Sort the keys, so that the results are deterministic
        keys = list(r1_dict.keys())
        keys.sort()
        fastq_files = []
        for key in keys:
            fastq_files.append((r1_dict[key], r2_dict[key]))
        return fastq_files

    def save_statistics(self, csv_file_path, sample_name=None, header=None):
        """Saves the demultiplex statistics into a CSV file

        """

        if not header:
            header = ""

        if not sample_name:
            sample_name = ""

        with open(csv_file_path, "w") as f:
            writer = csv.writer(f)
            writer.writerow([
                "sample", "barcode",
                "read1_percent", "read2_percent", "total_percent_%s" % header,
                "total_reads", "%s_reads" % header, "non%s_reads" % header
            ])
            if "total" not in self.counts:
                raise KeyError("Total read count is missing.")
            total = self.counts.get("total")
            if "unmatched" not in self.counts:
                raise KeyError("Unmatched read count is missing.")
            unmatched = self.counts.get("unmatched")
            for adapter in self.adapters:
                r1 = self.counts.get("%s_1" % adapter, 0)
                r2 = self.counts.get("%s_2" % adapter, 0)
                matched = self.counts.get(adapter, 0)
                writer.writerow([
                    sample_name, adapter, r1 / total, r2 / total, matched / total, total, matched, unmatched
                ])
        print("Statistics saved to %s" % csv_file_path)
        return csv_file_path


class Demultiplex:
    """Base class for demultiplexing FASTQ files.
    """

    DEFAULT_ERROR_RATE = 0.1

    @staticmethod
    def parse_barcode_outputs(barcode_outputs):
        """Parses the barcode and output file prefix pairs specified as a list of strings like:
            "BARCODE=PREFIX", or "BARCODE_1 BARCODE_2=PREFIX"

        Args:
            barcode_outputs (list): A list of strings in the format of BARCODE=PREFIX

        Returns: A dictionary, where each key is a barcode, each value is the file_prefix.

        """
        barcode_dict = {}
        for output in barcode_outputs:
            arr = str(output).split("=", 1)
            barcode_list = arr[0].strip().split()
            file_prefix = arr[1] if len(arr) > 1 else None
            for barcode in barcode_list:
                barcode_dict[barcode] = file_prefix
        return barcode_dict

    def __init__(self, adapters, error_rate=None, score=1, penalty=10):
        self.adapters = adapters
        self.min_match_length = round(min([len(adapter) / 2 for adapter in adapters]))

        # Set the default values
        self.error_rate = error_rate if error_rate else self.DEFAULT_ERROR_RATE
        logger.debug("Error Rate: %s" % self.error_rate)
        self.score = int(score) if str(score).isdigit() else 1
        logger.debug("Score: %s" % self.score)
        self.penalty = int(penalty) if str(penalty).isdigit() else 10
        logger.debug("Penalty: %s" % self.penalty)

        # Stores the statistics
        # Each sub-class may have different keys for the counts dictionary.
        # The values should be integers.
        # Whether to use this in the sub-class is optional.
        self.counts = {}

    def create_score_matrix(self):
        return parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)

    def semi_global_distance(self, s1, s2):
        score_matrix = self.create_score_matrix()
        result = parasail.sg_de_stats(
            s1, s2, self.penalty, self.penalty, score_matrix
        )
        return (self.score * result.matches - result.score) / self.penalty

    def demultiplex_fastq_pair(self, r1, r2, barcode_dict, ident=None):
        raise NotImplementedError()

    def update_counts(self, counts):
        """Updates the demultiplex statistics.

        Args:
            counts: A dictionary storing the statistics.

        This method merges counts from the input argument with self.counts
        Values of the existing keys are added together.
        Keys from counts will be created if they are not already in self.counts

        """
        for k, v in counts.items():
            self.counts[k] = self.counts.get(k, 0) + v
        return self

    @staticmethod
    def add_count(counts, key):
        """Increments the value of a particular key in counts (dictionary)
        This is a static method and it does NOT modify the self.counts
        """
        c = counts.get(key, 0)
        c += 1
        counts[key] = c
        return counts

    @staticmethod
    def print_output(message, ident=None):
        """Prints a message with an identifier
        This is used for identifying the source of the message when running multi-processing.
        """
        if ident:
            print("%s: %s" % (ident, message))
        else:
            print(message)

    @staticmethod
    def get_identifier(fastq_file):
        """Gets the identifier of the first read in a FASTQ file
        """
        with dnaio.open(fastq_file) as f:
            for read in f:
                return read.name.split(" ", 1)[0]

    @staticmethod
    def pair_fastq_files(r1_list, r2_list):
        """Pairs the FASTQ files from two lists by the identifier of the first read
        """
        if not r1_list or not r2_list:
            raise ValueError("Invalid R1 or R2.\nR1: %s\nR2: %s" % (r1_list, r2_list))
        # Convert r1_list and r2_list to lists if needed.
        if not isinstance(r1_list, list):
            r1_list = [r1_list]
        if not isinstance(r2_list, list):
            r2_list = [r2_list]

        r1_dict = {DemultiplexInline.get_identifier(f): f for f in r1_list}
        r2_dict = {DemultiplexInline.get_identifier(f): f for f in r2_list}

        if r1_dict.keys() != r2_dict.keys():
            raise ValueError("Unable to pair R1 and R2.\nR1: %s\nR2: %s" % (r1_dict.keys(), r2_dict.keys()))

        # Sort the keys, so that the results are deterministic
        keys = list(r1_dict.keys())
        keys.sort()
        fastq_files = []
        for key in keys:
            fastq_files.append((r1_dict[key], r2_dict[key]))
        return fastq_files

    @staticmethod
    def concatenate_fastq(prefix_dict):
        """Concatenates the FASTQ files in a list of directories.
        """
        logger.debug("Concatenating files: %s" % prefix_dict)
        for prefix, prefix_list in prefix_dict.items():
            r1, r2 = DemultiplexWriter.paired_end_filenames(prefix)
            pair_list = [DemultiplexWriter.paired_end_filenames(p) for p in prefix_list]

            cmd = "cat %s > %s" % (" ".join(p[0] for p in pair_list), r1)
            os.system(cmd)
            cmd = "cat %s > %s" % (" ".join(p[1] for p in pair_list), r2)
            os.system(cmd)

    def multi_processing(self, fastq_files, barcode_dict):
        """Demultiplex a list of FASTQ file pairs. The corresponding FASTQ files will be concatenated.

        Args:
            fastq_files: A list of 2-tuples, each stores the paths of a pair of FASTQ files.
            barcode_dict (dict):

        Returns:
            A dictionary containing the statistics of the demultiplex process (see self.counts).

        """
        pairs_count = len(fastq_files)
        # Determine the worker pool
        pool_size = min(pairs_count, os.cpu_count())
        pool = multiprocessing.Pool(pool_size)
        jobs = []
        output_list = []
        with tempfile.TemporaryDirectory() as temp_dir:
            logger.debug("Demultiplexing using %s processes..." % pool_size)
            for i in range(pairs_count):
                fastq_pair = fastq_files[i]
                ident = "Process_%s" % i

                output_dict = {k: os.path.join(temp_dir, "%s_%s" % (k, ident)) for k in barcode_dict.keys()}
                output_list.append(output_dict)

                job = pool.apply_async(
                    self.demultiplex_fastq_pair,
                    (fastq_pair[0], fastq_pair[1], output_dict, ident)
                )
                jobs.append(job)

            # Wait for the jobs to be finished and collect the statistics.
            results = [job.get() for job in jobs]
            counts = dict()
            for r in results:
                for k, v in r.items():
                    counts[k] = counts.get(k, 0) + v

            # path_dict is a dict storing the final output paths as keys, and
            # each value is a list of files (paths) to be concatenated.
            path_dict = {}
            for barcode, file_path in barcode_dict.items():
                # The following code will merge the barcodes pointing to the same file path
                path_list = path_dict.get(file_path, [])
                path_list.extend([
                    output_dict.get(barcode)
                    for output_dict in output_list
                    if output_dict.get(barcode)
                ])
                path_dict[file_path] = path_list

            self.concatenate_fastq(path_dict)
        return counts

    def run_demultiplex(self, r1_list, r2_list, barcode_dict):
        """Demultiplex by barcode adapters and concatenate FASTQ files.

        Args:
            r1_list (list): A list of file paths, each is a FASTQ with forward reads.
            r2_list (list): A list of file paths, each is a FASTQ with reverse compliment reads.
            barcode_dict:

        Returns:

        """
        # Pair the FASTQ files from r1_list and r2_list by the identifier of the first read
        fastq_files = self.pair_fastq_files(r1_list, r2_list)

        pairs_count = len(fastq_files)
        if pairs_count == 1:
            counts = self.demultiplex_fastq_pair(fastq_files[0][0], fastq_files[0][1], barcode_dict)
        else:
            counts = self.multi_processing(fastq_files, barcode_dict)
        self.update_counts(counts)
        logger.info(self.counts)
        return counts


class DemultiplexInline(Demultiplex):
    """Demultiplex FASTQ files by Inline barcode at the beginning of the reads

    self.counts will be a dictionary storing the demultiplex statistics, which includes the following keys:
        matched: the number of read pairs matched at least ONE of the adapters
        unmatched: the number of read pairs matching NONE of the adapters
        total: total number of read pairs processed
        Three keys for each adapter, i.e. BARCODE, BARCODE_1 and BARCODE_2
        BARCODE stores the number of reads matching the corresponding adapter.
        BARCODE_1 stores the number of forward reads (pair 1) matching the corresponding adapter.
        BARCODE_2 stores the number of reverse-compliment reads (pair 2) matching the corresponding adapter.

        For the values corresponding to BARCODE keys,
            each READ PAIR will only be counted once even if it is matching multiple adapters.
            The longer barcode will be used if the two reads in a pair is matching different adapters.

        For the values corresponding to BARCODE_1 and BARCODE_2 keys,
            each READ will be counted once.
            The first barcode in self.adapters will be used if the read is matching multiple adapters.

    """

    DEFAULT_ERROR_RATE = 0.2

    def trim_adapters(self, read1, read2, score_matrix=None):
        """Checks if the beginning of the reads in a read pair matches any adapter.
        If so, trim the reads to remove the matching adapter.

        A read and an adapter are matched by using semi-global alignment without penalty at the end of the read.
        They are considered as MATCHED if the number of substitutions and gaps

        Args:
            read1: Forward read.
            read2: Reverse Compliment read.
            score_matrix: Score/Substitution matrix for semi-global alignment.

        read1 and read2 are dnaio Sequence objects.
        They are passed into this method as references.
        The modifications on read1 and read2 will be preserved after return.

        The score_matrix will be generated if not specified.
        However, if specified, it must be generated by:
            score_matrix = parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)
        Pass a score_matrix avoid this method to generate it every time, which may speed up the overall processing.

        Returns: A 2-tuple indicating whether any adapters are matching the read pair.
            If a read matched an adapter, the matching adapter will be returned in the tuple.
            Otherwise, the corresponding element in the tuple will be None.
            For example, ("ACTGACT", None) indicates barcode "ACTGACT" is matching read1 (forward read).

        See Also:
             https://github.com/marcelm/dnaio/blob/master/src/dnaio/_core.pyx
             https://github.com/jeffdaily/parasail-python#substitution-matrices

        """
        if not score_matrix:
            score_matrix = self.create_score_matrix()

        # Trim both read1 and read2 with all adapters before return
        reads = [read1, read2]
        # Indicates whether R1 or R2 matches the adapter.
        matched = [""] * len(reads)
        for i in range(len(reads)):
            read = reads[i]
            for adapter in self.adapters:
                result = parasail.sg_de_stats(
                    adapter, read.sequence[:20], self.penalty, self.penalty, score_matrix
                )
                if result.matches <= self.min_match_length:
                    continue
                distance = (self.score * result.matches - result.score) / self.penalty
                max_distance = math.floor(len(adapter) * self.error_rate)
                if distance <= max_distance:
                    matched[i] = adapter
                    read.sequence = read.sequence[result.end_ref + 1:]
                    read.qualities = read.qualities[result.end_ref + 1:]
                    break
        # read1 and read2 are preserved implicitly
        return matched[0], matched[1]

    def __process_read_pair(self, read1, read2, score_matrix, fp_dict, counts):
        # Initialize ReadPair to check if read1 and read2 are valid
        read1, read2 = ReadPair(read1, read2).reads
        # read1 and read2 are references
        # The modifications on read1 and read2 will be returned implicitly
        adapter1, adapter2 = self.trim_adapters(read1, read2, score_matrix)

        if adapter1:
            self.add_count(counts, "%s_1" % adapter1)
        if adapter2:
            self.add_count(counts, "%s_2" % adapter2)

        # The longer adapter has higher priority
        adapter = adapter1 if len(adapter1) > len(adapter2) else adapter2
        if adapter:
            # Count the number of reads matching the longer adapter
            self.add_count(counts, adapter)
            # Sequence matched a barcode
            self.add_count(counts, 'matched')
            out_match = fp_dict.get(adapter)
            if out_match:
                out_match.write(read1, read2)
        else:
            # Sequence does not match a barcode
            self.add_count(counts, 'unmatched')
            out_unmatch = fp_dict.get('NO_MATCH')
            if out_unmatch:
                out_unmatch.write(read1, read2)

    def demultiplex_fastq_pair(self, r1, r2, barcode_dict, ident=None):
        """Demultiplex a single pair of FASTQ files.

        Reads matching any barcode in self.adapters will be extracted to a pair of FASTQ files.
        While the other reads will be extracted to another pair.

        Args:
            r1: FASTQ R1 path.
            r2: FASTQ R2 path.
            barcode_dict (dict): A dictionary where each key is a barcode and each value is the output file path.
            ident: Identifier for the demultiplex process.

        Returns:
            A dictionary containing the statistics of the demultiplex process (see self.counts).

        """
        logger.debug("Adapters: %s" % self.adapters)
        counter = 0
        counts = dict(matched=0, unmatched=0, total=0)

        with DemultiplexWriter(barcode_dict) as fp_dict:
            # score_matrix is initialized here because it cannot be pickled
            score_matrix = self.create_score_matrix()

            self.print_output("Demultiplexing %s and %s..." % (r1, r2), ident)
            with dnaio.open(r1, file2=r2) as fastq_in:
                for read1, read2 in fastq_in:
                    self.__process_read_pair(read1, read2, score_matrix, fp_dict, counts)
                    counter += 1
                    if counter % 100000 == 0:
                        self.print_output("%s reads processed." % counter, ident)
            self.print_output("%s reads processed." % counter, ident)
            self.print_output("%s reads matched." % counts.get("matched", 0), ident)
            self.print_output("%s reads unmatched." % counts.get("unmatched", 0), ident)
            self.print_output("Output Prefixes:\n%s" % "\n".join(fp_dict.prefix_dict.keys()), ident)

        counts['total'] = counter
        return counts

    def save_statistics(self, csv_file_path, name=""):
        """Saves the demultiplex statistics into a CSV file
        """
        with open(csv_file_path, "w") as f:
            writer = csv.writer(f)
            writer.writerow([
                "sample", "barcode",
                "read1_percent", "read2_percent", "total_percent_rna",
                "total_reads", "rna_reads", "nonrna_reads"
            ])
            if "total" not in self.counts:
                raise KeyError("Total read count is missing.")
            total = self.counts.get("total")
            if "unmatched" not in self.counts:
                raise KeyError("Unmatched read count is missing.")
            unmatched = self.counts.get("unmatched")
            for adapter in self.adapters:
                r1 = self.counts.get("%s_1" % adapter, 0)
                r2 = self.counts.get("%s_2" % adapter, 0)
                matched = self.counts.get(adapter, 0)
                writer.writerow([name, adapter, r1 / total, r2 / total, matched / total, total, matched, unmatched])
        self.print_output("Statistics saved to %s" % csv_file_path)
        return csv_file_path


class DemultiplexBarcode(Demultiplex):
    def __init__(self, adapters, error_rate=None, score=1, penalty=10):
        super().__init__(adapters, error_rate, score, penalty)
        self.max_error = {adapter: math.floor(len(adapter) * self.error_rate) for adapter in self.adapters}

    def match_adapters(self, barcode):
        # barcode_i7, barcode_i5 = barcode.split("+", 1)

        for adapter in self.adapters:
            mismatch = editdistance.eval(barcode, adapter)
            # adapter_i7, adapter_i5 = adapter.split("+", 1)

            # mismatch = self.semi_global_distance(barcode_i7, adapter_i7) + \
            #     self.semi_global_distance(barcode_i5, adapter_i5)

            # mismatch = editdistance.eval(barcode_i7, adapter_i7) + editdistance.eval(barcode_i5, adapter_i5)
            if mismatch < self.max_error.get(adapter, 0):
                return barcode
        return None

    @staticmethod
    def determine_adapters(r1):
        barcode_dict = dict()
        counter = 0
        with dnaio.open(r1) as fastq_in:
            for read1 in fastq_in:
                barcode = read1.name.strip().split(":")[-1]
                if re.match(IlluminaFASTQ.dual_index_pattern, barcode):
                    barcode = IlluminaFASTQ.convert_barcode(barcode)
                DemultiplexBarcode.add_count(barcode_dict, barcode)
                counter += 1
                if counter > 3000:
                    break
        barcode_list = BarcodeStatistics(barcode_dict).major_barcodes()
        return barcode_list

    def demultiplex_fastq_pair(self, r1, r2, barcode_dict, ident=None):
        counts = dict()
        counter = 0
        with DemultiplexWriter(barcode_dict) as fp_dict:
            with dnaio.open(r1, file2=r2) as fastq_in:
                for read1, read2 in fastq_in:
                    counter += 1
                    if counter % 100000 == 0:
                        self.print_output("%s reads processed." % counter, ident)
                    barcode = ReadPair(read1, read2).barcode
                    if re.match(IlluminaFASTQ.dual_index_pattern, barcode):
                        barcode = IlluminaFASTQ.convert_barcode(barcode)
                    barcode = self.match_adapters(barcode)
                    if not barcode:
                        # TODO: Write unmatched reads.
                        self.add_count(counts, "unmatched")
                        continue
                    self.add_count(counts, barcode)
                    file_obj = fp_dict.get(barcode)
                    if file_obj:
                        # Write the barcode line
                        file_obj.write(read1, read2)

        counts["total"] = counter
        return counts
