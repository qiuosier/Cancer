import multiprocessing
import os
import shutil
from .fastq_pair import FASTQPair


def demux_fastq_pair(r1, r2, adapters, output_dir, error_rate=0.2, ident=None):
    p = FASTQPair(r1, r2, ident)
    p.extract_reads_by_adapters(adapters, output_dir, error_rate)


def demux_inline(fastq_files, adapters, output_dir, error_rate=0.2):
    """

    Args:
        fastq_files (list): A list of 2-tuples, each is a pair of R1 and R2.
        adapters:
        output_dir:
        error_rate:

    Returns:

    """
    if not isinstance(fastq_files, list):
        fastq_files = [fastq_files]
    pairs_count = len(fastq_files)
    if pairs_count == 1:
        demux_fastq_pair(fastq_files[0][0], fastq_files[0][1], adapters, output_dir, error_rate)
        return

    pool_size = min(pairs_count, os.cpu_count())
    pool = multiprocessing.Pool(pool_size)
    jobs = []

    output_dir_list = []

    print("Demultiplexing using %s processes..." % pool_size)
    for i in range(pairs_count):
        fastq_pair = fastq_files[i]
        sub_output_dir = os.path.join(output_dir, str(i))
        if not os.path.exists(sub_output_dir):
            os.mkdir(sub_output_dir)
        output_dir_list.append(sub_output_dir)
        ident = "Process %s" % i
        job = pool.apply_async(
            demux_fastq_pair,
            (fastq_pair[0], fastq_pair[1], adapters, sub_output_dir, error_rate, ident)
        )
        jobs.append(job)

    print("Demultiplexing to %s..." % output_dir)
    results = [job.get() for job in jobs]

    print("Concatenating %s pairs of output files..." % pairs_count)

    filenames = [
        FASTQPair.r1_matched_filename,
        FASTQPair.r2_matched_filename,
        FASTQPair.r1_unmatched_filename,
        FASTQPair.r2_unmatched_filename
    ]

    for filename in filenames:
        cmd = "cat %s > %s" % (
            " ".join([os.path.join(out, filename) for out in output_dir_list]),
            os.path.join(output_dir, filename)
        )
        os.system(cmd)
    # for sub_output_dir in output_dir_list:
    #     shutil.rmtree(sub_output_dir)
