import multiprocessing
import os
import csv
import math
import dnaio
import parasail
import editdistance
import re
import logging
from .fastq_pair import ReadPair
from .fastq import IlluminaFASTQ, BarcodeStatistics
logger = logging.getLogger(__name__)


class Demultiplex:
    """Base class for demultiplexing FASTQ files.
    """

    DEFAULT_ERROR_RATE = 0.1

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

    @property
    def output_filenames(self):
        raise NotImplementedError()

    def demultiplex_fastq_pair(self, r1, r2, output_dir, ident=None):
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

    def concatenate_fastq(self, output_dir_list, output_dir):
        """Concatenates the FASTQ files in a list of directories.
        """
        logger.debug("Concatenating %s pairs of output files..." % len(output_dir_list))
        for filename in self.output_filenames:
            cmd = "cat %s > %s" % (
                " ".join([
                    os.path.join(out, filename) for out in output_dir_list
                    if os.path.exists(os.path.join(out, filename))
                ]),
                os.path.join(output_dir, filename)
            )
            os.system(cmd)

    def multi_processing(self, fastq_files, output_dir):
        """Demultiplex a list of FASTQ file pairs. The corresponding FASTQ files will be concatenated.

        Args:
            fastq_files: A list of 2-tuples, each stores the paths of a pair of FASTQ files.
            output_dir: Output directory.

        Returns:
            A dictionary containing the statistics of the demultiplex process (see self.counts).

        """
        pairs_count = len(fastq_files)
        # Determine the worker pool
        pool_size = min(pairs_count, os.cpu_count())
        pool = multiprocessing.Pool(pool_size)
        jobs = []
        output_dir_list = []

        logger.debug("Demultiplexing using %s processes..." % pool_size)
        for i in range(pairs_count):
            fastq_pair = fastq_files[i]
            sub_output_dir = os.path.join(output_dir, str(i))
            if not os.path.exists(sub_output_dir):
                os.mkdir(sub_output_dir)
            output_dir_list.append(sub_output_dir)
            ident = "Process %s" % i
            job = pool.apply_async(
                self.demultiplex_fastq_pair,
                (fastq_pair[0], fastq_pair[1], sub_output_dir, ident)
            )
            jobs.append(job)

        # Wait for the jobs to be finished and collect the statistics.
        results = [job.get() for job in jobs]
        counts = dict()
        for r in results:
            for k, v in r.items():
                counts[k] = counts.get(k, 0) + v

        self.concatenate_fastq(output_dir_list, output_dir)
        # TODO: Remove temp outputs
        return counts

    def run_demultiplex(self, r1_list, r2_list, output_dir):
        """Demultiplex by inline barcode adapters and concatenate FASTQ files.

        Args:
            r1_list (list): A list of file paths, each is a FASTQ with forward reads.
            r2_list (list): A list of file paths, each is a FASTQ with reverse compliment reads.
            output_dir: Output directory

        Returns:

        """
        # Pair the FASTQ files from r1_list and r2_list by the identifier of the first read
        fastq_files = self.pair_fastq_files(r1_list, r2_list)

        # Creates output directory if it does not exist.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        logger.info("Demultiplexing to %s..." % output_dir)

        pairs_count = len(fastq_files)
        if pairs_count == 1:
            counts = self.demultiplex_fastq_pair(fastq_files[0][0], fastq_files[0][1], output_dir)
        else:
            counts = self.multi_processing(fastq_files, output_dir)
        self.update_counts(counts)
        logger.info(self.counts)
        return counts


class DemultiplexInline(Demultiplex):
    """Demultiplex FASTQ files by Inline barcode at the beginning of the reads

    Attributes:
        counts: A dictionary storing the demultiplex statistics, which includes the following keys:
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

    # Output Filenames
    r1_matched_filename = "R1_matched.fastq.gz"
    r2_matched_filename = "R2_matched.fastq.gz"
    r1_unmatched_filename = "R1_unmatched.fastq.gz"
    r2_unmatched_filename = "R2_unmatched.fastq.gz"

    DEFAULT_ERROR_RATE = 0.2

    @property
    def output_filenames(self):
        return [
            self.r1_matched_filename,
            self.r2_matched_filename,
            self.r1_unmatched_filename,
            self.r2_unmatched_filename
        ]

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
            score_matrix = parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)

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

    def __process_read_pair(self, read1, read2, counts, out_match, out_unmatch, score_matrix):
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
            out_match.write(read1, read2)
        else:
            # Sequence does not match a barcode
            self.add_count(counts, 'unmatched')
            out_unmatch.write(read1, read2)

    def demultiplex_fastq_pair(self, r1, r2, output_dir, ident=None):
        """Demultiplex a single pair of FASTQ files.

        Reads matching any barcode in self.adapters will be extracted to a pair of FASTQ files.
        While the other reads will be extracted to another pair.

        Args:
            r1: FASTQ R1 path.
            r2: FASTQ R2 path.
            output_dir: Output directory for the output files.
            ident: Identifier for the demultiplex process.

        Returns:
            A dictionary containing the statistics of the demultiplex process (see self.counts).

        """
        logger.debug("Adapters: %s" % self.adapters)
        counter = 0
        counts = dict(matched=0, unmatched=0, total=0)
        r1_match_path = os.path.join(output_dir, self.r1_matched_filename)
        r2_match_path = os.path.join(output_dir, self.r2_matched_filename)
        r1_unmatch_path = os.path.join(output_dir, self.r1_unmatched_filename)
        r2_unmatch_path = os.path.join(output_dir, self.r2_unmatched_filename)

        # score_matrix is initialized here because it cannot be pickled
        score_matrix = parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)

        self.print_output("Demultiplexing %s and %s..." % (r1, r2), ident)
        with dnaio.open(r1, file2=r2) as fastq_in, \
                dnaio.open(r1_match_path, file2=r2_match_path, mode='w') as out_match, \
                dnaio.open(r1_unmatch_path, file2=r2_unmatch_path, mode='w') as out_unmatch:
            for read1, read2 in fastq_in:
                self.__process_read_pair(read1, read2, counts, out_match, out_unmatch, score_matrix)
                counter += 1
                if counter % 100000 == 0:
                    self.print_output("%s reads processed." % counter, ident)
        self.print_output("%s reads processed." % counter, ident)
        self.print_output("%s reads matched." % counts.get("matched", 0), ident)
        self.print_output("%s reads unmatched." % counts.get("unmatched", 0), ident)
        self.print_output("Output Files:\n%s" % "\n".join([
            r1_match_path, r2_match_path, r1_unmatch_path, r2_unmatch_path
        ]), ident)

        counts['total'] = counter
        return counts

    def save_statistics(self, output_dir, name=""):
        """Saves the demultiplex statistics into a CSV file
        """
        output_csv = os.path.join(output_dir, "all_barcode_stats.csv")
        with open(output_csv, "w") as f:
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
        return output_csv


class DemultiplexBarcode(Demultiplex):
    def __init__(self, adapters, error_rate=None, score=1, penalty=10):
        super().__init__(adapters, error_rate, score, penalty)
        self.max_error = {adapter: math.floor(len(adapter) * self.error_rate) for adapter in self.adapters}

    @property
    def output_filenames(self):
        r1_list = ["R1_%s.fastq.gz" % adapter for adapter in self.adapters]
        r2_list = ["R2_%s.fastq.gz" % adapter for adapter in self.adapters]
        return r1_list + r2_list

    def match_adapters(self, barcode):
        for adapter in self.adapters:
            if editdistance.eval(barcode, adapter) < self.max_error.get(adapter, 0):
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

    def demultiplex_fastq_pair(self, r1, r2, output_dir, ident=None):
        counts = dict()
        # A dictionary to hold the opened file obj
        file_obj_dict = dict()
        counter = 0
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
                    self.add_count(counts, "unmatched")
                    continue
                self.add_count(counts, barcode)
                file_obj = file_obj_dict.get(barcode)
                if not file_obj:
                    r1_file_path = os.path.join(output_dir, "R1_%s.fastq.gz" % barcode)
                    r2_file_path = os.path.join(output_dir, "R2_%s.fastq.gz" % barcode)
                    logger.debug("Creating files: %s" % [r1_file_path, r2_file_path])
                    file_obj = dnaio.open(r1_file_path, file2=r2_file_path, mode='w')
                    file_obj_dict[barcode] = file_obj
                # Write the barcode line
                file_obj.write(read1, read2)
        for file_obj in file_obj_dict.values():
            file_obj.close()
        counts["total"] = counter
        return counts
