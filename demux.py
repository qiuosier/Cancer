import multiprocessing
import os
import csv
import math
import dnaio
import parasail
from threading import Lock
from .fastq_pair import FASTQPair


class Demultiplex:
    def __init__(self, adapters, error_rate=0.2, score=1, penalty=10):
        self.adapters = adapters
        self.min_match_length = round(min([len(adapter) / 2 for adapter in adapters]))

        self.error_rate = error_rate if error_rate else 0.2
        self.score = int(score) if str(score).isdigit() else 1
        print("Score: %s" % self.score)
        self.penalty = int(penalty) if str(penalty).isdigit() else 10
        print("Penalty: %s" % self.penalty)

        self.counts = {}


class DemultiplexInline(Demultiplex):
    def update_counts(self, counts):
        for k, v in counts.items():
            self.counts[k] = self.counts.get(k, 0) + v

    @staticmethod
    def print_output(message, ident=None):
        if ident:
            print("%s: %s" % (ident, message))
        else:
            print(message)

    def trim_adapters(self, read1, read2, score_matrix):
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
        return matched[0], matched[1]

    def demultiplex_fastq_pair(self, r1, r2, output_dir, ident=None):
        print("Adapters: %s" % self.adapters)
        counter = 0
        counter_matched = 0
        counts = dict()
        counter_unmatched = 0
        r1_match_path = os.path.join(output_dir, "R1_matched.fastq.gz")
        r2_match_path = os.path.join(output_dir, "R2_matched.fastq.gz")
        r1_unmatch_path = os.path.join(output_dir, "R1_unmatched.fastq.gz")
        r2_unmatch_path = os.path.join(output_dir, "R2_unmatched.fastq.gz")

        score_matrix = parasail.matrix_create("ACGTN", self.score, -1 * self.penalty)

        with dnaio.open(r1, file2=r2) as fastq_in, \
                dnaio.open(r1_match_path, file2=r2_match_path, mode='w') as out_match, \
                dnaio.open(r1_unmatch_path, file2=r2_unmatch_path, mode='w') as out_unmatch:
            for read1, read2 in fastq_in:
                # adapter, trimmed_read1, trimmed_read2 = self.__match_adapters(read1, read2, adapters, error_rate)
                adapter1, adapter2 = self.trim_adapters(read1, read2, score_matrix)

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

                    # Sequence matched a barcode
                    counter_matched += 1
                    out_match.write(read1, read2)
                else:
                    # Sequence does not match a barcode
                    counter_unmatched += 1
                    out_unmatch.write(read1, read2)
                counter += 1
                if counter % 100000 == 0:
                    self.print_output("%s reads processed." % counter, ident)
        self.print_output("%s reads processed." % counter, ident)
        self.print_output("%s reads matched." % counter_matched, ident)
        self.print_output("%s reads unmatched." % counter_unmatched, ident)
        self.print_output("Output Files:\n%s" % "\n".join([
            r1_match_path, r2_match_path, r1_unmatch_path, r2_unmatch_path
        ]), ident)
        counts['unmatched'] = counter_unmatched
        counts['total'] = counter
        return counts

    def multi_processing(self, fastq_files, output_dir):
        pairs_count = len(fastq_files)
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
                self.demultiplex_fastq_pair,
                (fastq_pair[0], fastq_pair[1], sub_output_dir, ident)
            )
            jobs.append(job)

        results = [job.get() for job in jobs]
        counts = dict()
        for r in results:
            for k, v in r.items():
                counts[k] = counts.get(k, 0) + v

        self.concatenate_fastq(output_dir_list, output_dir)
        # TODO: Remove temp outputs
        return counts

    @staticmethod
    def concatenate_fastq(output_dir_list, output_dir):
        print("Concatenating %s pairs of output files..." % len(output_dir_list))
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

    def run_demultiplex(self, fastq_files, output_dir):
        """

        Args:
            fastq_files (list): A list of 2-tuples, each is a pair of R1 and R2.
            output_dir:

        Returns:

        """
        if not isinstance(fastq_files, list):
            fastq_files = [fastq_files]
        # Creates output directory if it does not exist.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        print("Demultiplexing to %s..." % output_dir)

        pairs_count = len(fastq_files)
        if pairs_count == 1:
            counts = self.demultiplex_fastq_pair(fastq_files[0][0], fastq_files[0][1], output_dir)
        else:
            counts = self.multi_processing(fastq_files, output_dir)
        self.update_counts(counts)
        return counts

    def save_statistics(self, output_dir, name=""):
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
        print(self.counts)
        return output_csv
