import datetime
import argparse
import os
from .fastq_pair import FASTQPair
from .demux import demux_inline
from .variants import files


class Program:
    @staticmethod
    def demux_inline(args):
        if len(args.r1) != len(args.r2):
            raise ValueError("R1 and R2 must have the same number of files.")
        fastq_files = []
        for i in range(len(args.r1)):
            fastq_files.append((args.r1[i], args.r2[i]))
        adapters = [s.strip() for s in args.barcode]
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        demux_inline(fastq_files, adapters, args.output, error_rate=args.error_rate, name=args.name)

    @staticmethod
    def compare_fastq(args):
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        FASTQPair(*args.FASTQ).diff(args.compare[0], args.compare[1], args.output)

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


def main():
    parser = argparse.ArgumentParser(description="Command line entry points to Cancer package.")
    # parser.add_argument("program", nargs=1, type=str, help="Program name")
    subparsers = parser.add_subparsers(title="Program", help="Program", dest='program')

    sub_parser = subparsers.add_parser("demux_inline", help="Demultiplex FASTQ files using Inline Barcodes")
    sub_parser.add_argument('--r1', nargs='+', required=True, help="FASTQ R1 files")
    sub_parser.add_argument('--r2', nargs='+', required=True, help="FASTQ R2 files")
    sub_parser.add_argument('--barcode', nargs='+', required=True, help="Inline Barcodes")
    sub_parser.add_argument('--output', required=True, help="Output Directory")
    sub_parser.add_argument('--error_rate', type=float, help="Max Error Allowed")
    sub_parser.add_argument('--name', type=str, help="Sample Name for statistics")

    sub_parser = subparsers.add_parser("compare_fastq", help="Compare reads in two pairs of FASTQ files.")
    sub_parser.add_argument('FASTQ', nargs=2, help="FASTQ R1 and R2 files")
    sub_parser.add_argument('--compare', nargs=2, help="FASTQ R1 and R2 files")
    sub_parser.add_argument('--output', required=True, help="Output Directory")

    sub_parser = subparsers.add_parser("filter_whitelist", help="Filter the variants using a whitelist.")
    sub_parser.add_argument('whitelist', help="VCF/CSV/TSV File containing whitelist variants.")
    sub_parser.add_argument('vcf', help="VCF file containing the variants to be filtered.")
    sub_parser.add_argument('--output', required=True, help="Output Directory")

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
    main()
