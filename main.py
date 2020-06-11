import datetime
import argparse
import os
import logging
import json
from .fastq_file import ReadIdentifier
from .fastq_pair import FASTQPair
from .fastq.demux import DemultiplexInline, DemultiplexDualIndex, DemultiplexWriter
from .fastq.barcode import BarcodeCounter
from .variants import files
from Aries.outputs import LoggingConfig
logger = logging.getLogger(__name__)


class Program:
    """Contains static methods for sub-commands to start the processing program with args
    """
    @staticmethod
    def demux_inline(args):
        """Demultiplex FASTQ files with inline barcode adapters.
        """
        if len(args.r1) != len(args.r2):
            raise ValueError("R1 and R2 must have the same number of files.")
        barcode_dict = DemultiplexInline.parse_barcode_outputs(args.barcode)
        barcode_dict[DemultiplexWriter.BARCODE_NOT_MATCHED] = args.unmatched

        demux_inline = DemultiplexInline(
            barcode_dict, error_rate=args.error_rate, score=args.score, penalty=args.penalty
        )
        fastq_files = demux_inline.pair_fastq_files(args.r1, args.r2)
        demux_inline.start(fastq_files)

        if args.stats:
            demux_inline.save_statistics(args.stats, sample_name=args.name, header=args.header)

    @staticmethod
    def demux_barcode(args):
        if len(args.r1) != len(args.r2):
            raise ValueError("R1 and R2 must have the same number of files.")

        if args.barcode:
            adapters = [s.strip() for s in args.barcode]
        else:
            if isinstance(args.r1, list):
                r1 = args.r1[0]
            else:
                r1 = args.r1
            logger.debug("Determining the barcodes in %s..." % r1)
            adapters = DemultiplexDualIndex.determine_adapters(r1)
        logger.debug("Barcodes: %s" % adapters)
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        # Use barcode as output file prefix
        barcode_dict = {adapter: os.path.join(args.output, adapter) for adapter in adapters}
        demux_dual = DemultiplexDualIndex(
            barcode_dict, error_rate=args.error_rate, score=args.score, penalty=args.penalty
        )
        fastq_files = demux_dual.pair_fastq_files(args.r1, args.r2)
        demux_dual.start(fastq_files)

    @staticmethod
    def compare_fastq(args):
        """Compare FASTQ files
        """
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        FASTQPair(*args.FASTQ).diff(args.compare[0], args.compare[1], args.output, args.chunk_size)

    @staticmethod
    def count_inline_barcode(args):
        """Counts the inline barcode
        """
        if args.r2:
            if len(args.r1) != len(args.r2):
                raise ValueError("R1 and R2 must have the same number of files.")
            fastq_files = BarcodeCounter.pair_fastq_files(args.r1, args.r2)
        else:
            fastq_files = args.r1

        BarcodeCounter(args.start, args.length).start(fastq_files)
        return

    @staticmethod
    def filter_whitelist(args):
        whitelist_path = str(args.whitelist)
        if whitelist_path.endswith(".vcf"):
            whitelist = files.VCFVariants(whitelist_path)
            out_whitelist = os.path.join(args.output, files.WhitelistFilter.whitelist_output_filename)
        elif whitelist_path.endswith(".csv") or whitelist_path.endswith(".tsv"):
            whitelist = files.CSVVariants(whitelist_path)
            out_whitelist = os.path.join(args.output, files.WhitelistFilter.whitelist_output_filename)
        else:
            raise TypeError("Whitelist file type is not supported: %s" % whitelist_path)
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        out_vcf = os.path.join(args.output, files.WhitelistFilter.vcf_output_filename)
        vcf = files.VCFVariants(args.vcf)
        whitelist_filter = files.WhitelistFilter(whitelist)
        print("%s variants in the white list." % len(whitelist_filter.index.keys()))
        description = "In Whitelist: %s" % whitelist_path
        vcf.apply_filter(out_vcf, "Whitelist", description, whitelist_filter.filter_variant, passed_only=True)
        print("%d whitelist variants found in VCF." % len(whitelist_filter.in_whitelist))
        whitelist_filter.print_passed()
        print("Saving whitelist subset to %s" % out_whitelist)
        if out_whitelist:
            whitelist_filter.save_passed(out_whitelist)

    @staticmethod
    def parse_read_identifier(args):
        read_id = ReadIdentifier(args.line)
        print(json.dumps(read_id.info, indent=4))


def main():
    parser = argparse.ArgumentParser(description="Command line entry points to Cancer package.")
    # parser.add_argument("program", nargs=1, type=str, help="Program name")
    subparsers = parser.add_subparsers(title="Program", help="Program", dest='program')

    sub_parser = subparsers.add_parser("demux_inline", help="Demultiplex FASTQ files using Inline Barcodes")
    sub_parser.add_argument('--r1', nargs='+', required=True, help="FASTQ R1 files")
    sub_parser.add_argument('--r2', nargs='+', required=True, help="FASTQ R2 files")
    sub_parser.add_argument('--barcode', nargs='+',
                            help="Barcode and the output file prefix in the format of BARCODE=PREFIX")
    sub_parser.add_argument('--unmatched', type=str, help="File path for saving the unmatched reads.")
    sub_parser.add_argument('--error_rate', type=float, help="Max Error Allowed, defaults to 20%%")
    sub_parser.add_argument('--score', type=int, help="Score for each bp matched, defaults to 1")
    sub_parser.add_argument('--penalty', type=int, help="Penalty for each bp mis-matched, defaults to 10")
    sub_parser.add_argument('--stats', type=str, help="Specify a CSV file path to save the statistics.")
    sub_parser.add_argument('--name', type=str, help="For statistics, Sample name")
    sub_parser.add_argument('--header', type=str, help="For statistics, the column header for matched read counts")

    sub_parser = subparsers.add_parser("demux_barcode", help="Demultiplex FASTQ files using Read Barcodes")
    sub_parser.add_argument('--r1', nargs='+', required=True, help="FASTQ R1 files")
    sub_parser.add_argument('--r2', nargs='+', required=True, help="FASTQ R2 files")
    sub_parser.add_argument('--barcode', nargs='+',
                            help="Barcode and the output file prefix in the format of BARCODE=PREFIX")
    sub_parser.add_argument('--output', required=True, help="Output Directory")
    sub_parser.add_argument('--name', type=str, help="Sample Name for statistics")
    sub_parser.add_argument('--error_rate', type=float, help="Max Error Allowed")
    sub_parser.add_argument('--score', type=int, help="Score for each bp matched")
    sub_parser.add_argument('--penalty', type=int, help="Penalty for each bp mis-matched")

    sub_parser = subparsers.add_parser("compare_fastq", help="Compare reads in two pairs of FASTQ files.")
    sub_parser.add_argument('FASTQ', nargs=2, help="First pair of FASTQ R1 and R2 files (P1)")
    sub_parser.add_argument('--compare', nargs=2, help="Second pair of FASTQ R1 and R2 files (P2)")
    sub_parser.add_argument('--output', required=True, help="Directory for storing output files")
    sub_parser.add_argument('--chunk_size', type=int, help="Chunks size in number of reads")

    sub_parser = subparsers.add_parser("filter_whitelist", help="Filter the variants using a whitelist.")
    sub_parser.add_argument('whitelist', help="VCF/CSV/TSV File containing whitelist variants.")
    sub_parser.add_argument('vcf', help="VCF file containing the variants to be filtered.")
    sub_parser.add_argument('--output', required=True, help="Output Directory")

    sub_parser = subparsers.add_parser("parse_read_identifier", help="Parses the identifier of a read from FASTQ file.")
    sub_parser.add_argument('line', help="Identifier line of a read from FASTQ file.")

    sub_parser = subparsers.add_parser("count_inline_barcode", help="Count the usage of Inline Barcodes")
    sub_parser.add_argument('--r1', nargs='+', required=True, help="FASTQ R1 files")
    sub_parser.add_argument('--r2', nargs='+', help="FASTQ R2 files")
    sub_parser.add_argument('-s', '--start', type=int, help="Starting position of the barcode (0-based)")
    sub_parser.add_argument('-l', '--length', type=int, required=True, help="Length of the barcode")

    args = parser.parse_args()
    # Show help if no subparser matched.
    if not vars(args).keys() or not args.program or not hasattr(Program, args.program):
        parser.parse_args(["-h"])
        return
    func = getattr(Program, args.program)
    start = datetime.datetime.now()
    print("Starting %s at %s" % (args.program, start))
    func(args)
    print("Total Time: %s" % (datetime.datetime.now() - start))


if __name__ == '__main__':
    with LoggingConfig():
        main()
