import datetime
import argparse
from .fastq_pair import FASTQPair
from .demux import demux_inline


def main():
    parser = argparse.ArgumentParser(description="Command line entry points to Cancer package.")
    # parser.add_argument("program", nargs=1, type=str, help="Program name")
    subparsers = parser.add_subparsers(title="Program", help="Program", dest='program')

    parser_demux = subparsers.add_parser("demux_inline", help="Demultiplex FASTQ files using Inline Barcodes")
    parser_demux.add_argument('--r1', nargs='+', required=True, help="FASTQ R1 files")
    parser_demux.add_argument('--r2', nargs='+', required=True, help="FASTQ R2 files")
    parser_demux.add_argument('--barcode', nargs='+', required=True, help="Inline Barcodes")
    parser_demux.add_argument('--output', required=True, help="Output Directory")

    parser_demux = subparsers.add_parser("compare_fastq", help="Compare reads in two pairs of FASTQ files.")
    parser_demux.add_argument('FASTQ', nargs=2, help="FASTQ R1 and R2 files")
    parser_demux.add_argument('--compare', nargs=2, help="FASTQ R1 and R2 files")
    parser_demux.add_argument('--output', required=True, help="Output Directory")

    args = parser.parse_args()
    # Show help if no subparser matched.
    if not vars(args).keys():
        parser.parse_args(["-h"])
        return

    program = args.program
    start = datetime.datetime.now()
    print("Starting %s at %s" % (program, start))
    if program == "demux_inline":
        if len(args.r1) != len(args.r2):
            raise ValueError("R1 and R2 must have the same number of files.")
        fastq_files = []
        for i in range(len(args.r1)):
            fastq_files.append((args.r1[i], args.r2[i]))
        adapters = [s.strip() for s in args.barcode]
        output_dir = args.output
        demux_inline(fastq_files, adapters, output_dir, error_rate=0.2)
    elif program == "compare_fastq":
        FASTQPair(*args.FASTQ).diff(args.compare[0], args.compare[1], args.output)

    print("Total Time: %s" % (datetime.datetime.now() - start))


if __name__ == '__main__':
    main()
