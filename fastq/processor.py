import multiprocessing
import os
import dnaio
import logging
import datetime
import time
from multiprocessing import Process, Manager
logger = logging.getLogger(__name__)


class FASTQWorker:
    """Represents a worker process for processing FASTQ reads
    """

    def __init__(self):
        """Initialize a worker process.

        """
        # Stores the statistics
        # Each sub-class may have different keys for the counts dictionary.
        # The values should be integers.
        # Whether to use this in the sub-class is optional.
        self.counts = dict(matched=0, unmatched=0, total=0)

    def add_count(self, key, val=1):
        """Increments the value of a particular key in counts (dictionary)
        This is a static method and it does NOT modify the self.counts
        """
        c = self.counts.get(key, 0)
        c += val
        self.counts[key] = c
        return self.counts

    def start(self, in_queue, out_queue):
        """Starts processing reads from in_queue.
        The number of reads processed are put into the out_queue for counting purpose.

        Args:
            in_queue: A queue holding list of reads to be processed.
                Each item in the in_queue is a list reads so that the frequency of access the queue are reduced.
            out_queue: A queue holding integers for counting purpose.

        Returns:

        """
        active_time = datetime.timedelta()
        batch_count = 0

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
                self.process_read_pair(read_pair)

            self.add_count('total', len(reads))
            batch_count += 1
            # Add processing time for this batch
            active_time += (datetime.datetime.now() - timer_started)
            out_queue.put(len(reads))

    def process_read_pair(self, read_pair):
        """Process the read pair
        """
        raise NotImplementedError


class FASTQProcessor:

    # The number of read pairs to be packed into each item in the processing queue.
    BATCH_SIZE = 5000

    @staticmethod
    def read_data(fastq_files, queue, pool_size):
        """Reads read pairs from FASTQ files in to a queue.
        Once finished reading the files, a number of None values will be put into the queue,
        which can be used by other processes to determine if there will be more items for processing.

        Args:
            fastq_files: A list of 2-tuples containing paired-end FASTQ filenames.
                The reads from the list of files will be read into the queue sequentially.
                To read files in parallel, start a read_data process for each pair of files.
            queue: A queue for holding reads to be processed.
                Each item in the queue will be a list of read pairs (2-tuples).
            pool_size: The number of processors for processing items in the queue.
                The same number of None values will be put into the queue at the end of reading the files.

        """
        enqueue_time = datetime.timedelta()
        batch_count = 0

        process_started = datetime.datetime.now()
        for fastq in fastq_files:
            if isinstance(fastq, str):
                r1 = fastq
                r2 = None
            else:
                r1 = fastq[0]
                r2 = fastq[1]
            with dnaio.open(r1, file2=r2) as fastq_in:
                size = 0
                reads = []
                for read in fastq_in:
                    reads.append(read)
                    size += 1
                    if size > FASTQProcessor.BATCH_SIZE:
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
            datetime.datetime.now() - process_started, enqueue_time, FASTQProcessor.BATCH_SIZE, batch_count
        ))
        for i in range(pool_size):
            queue.put(None)

    @staticmethod
    def get_identifier(fastq_file):
        """Gets the identifier of the first read in a FASTQ file
        """
        with dnaio.open(fastq_file) as f:
            for read in f:
                return read.name.split(" ", 1)[0]

    @classmethod
    def pair_fastq_files(cls, r1_list, r2_list):
        """Pairs the FASTQ files from two lists by the identifier of the first read
        """
        if not r1_list or not r2_list:
            raise ValueError("Invalid R1 or R2.\nR1: %s\nR2: %s" % (r1_list, r2_list))
        # Convert r1_list and r2_list to lists if needed.
        if not isinstance(r1_list, list):
            r1_list = [r1_list]
        if not isinstance(r2_list, list):
            r2_list = [r2_list]

        r1_dict = {cls.get_identifier(f): f for f in r1_list}
        r2_dict = {cls.get_identifier(f): f for f in r2_list}

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
    def print_queue_size(q_size_array):
        for s in q_size_array:
            print("*" * s + str(s))

    @staticmethod
    def start_worker(worker_class, in_queue, out_queue, *args, **kwargs):
        return worker_class(*args, **kwargs).start(in_queue, out_queue)

    def __init__(self, worker_class, *args, **kwargs):
        self.worker_class = worker_class
        self.worker_args = args
        self.worker_kwargs = kwargs

        self.pool_size = max(os.cpu_count(), 1)
        logger.debug("Pool Size: %s" % self.pool_size)
        self.manager = Manager()
        self.pool = []
        self.counts = dict()

        self.reader_queue = self.manager.Queue(self.pool_size * 100)
        self.worker_queue = self.manager.Queue()

    def update_counts(self, counts):
        """Updates the demultiplex statistics.

        Args:
            counts: A dictionary storing the statistics. All values should be numeric.

        This method merges counts from the input argument with self.counts
        Values of the existing keys are added together.
        Keys from counts will be created if they are not already in self.counts

        """
        for k, v in counts.items():
            self.counts[k] = self.counts.get(k, 0) + v
        return self.counts

    def start_readers(self, fastq_files):
        readers = []
        for fastq_pair in fastq_files:
            reader = Process(
                target=FASTQProcessor.read_data,
                args=([fastq_pair], self.reader_queue, self.pool_size)
            )
            reader.start()
            readers.append(reader)
        return readers

    def finalize_readers(self, readers, reader_q_size):
        for reader in readers:
            reader.join()
        self.print_queue_size(reader_q_size)
        logger.debug("Reader's queue max size: %s" % max(reader_q_size))

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
            print("{:,} reads/pairs processed.".format(counter))
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

    def start(self, fastq_files):
        pool = multiprocessing.Pool(self.pool_size)
        jobs = []

        for i in range(self.pool_size):

            job = pool.apply_async(
                self.start_worker,
                (self.worker_class, self.reader_queue, self.worker_queue, *self.worker_args),
                self.worker_kwargs
            )
            jobs.append(job)

        readers = self.start_readers(fastq_files)

        # Wait for the jobs to finish and keep track of the queue size
        reader_q_size = self.wait_for_jobs(jobs)

        # Collect the statistics.
        self.collect_results(jobs)
        self.finalize_readers(readers, reader_q_size)

        return self
