import sys
import datetime
from .fastq_pair import FASTQPair


def main(program, *args):
    start = datetime.datetime.now()
    if program == "demux_inline":
        r1 = args[0]
        r2 = args[1]
        adapters = [s.strip() for s in str(args[2]).split(",")]
        output_dir = args[3]
        p = FASTQPair(r1, r2)
        p.extract_reads_by_adapters(adapters, output_dir, error_rate=0.2)
    print("Total Time: %s" % (datetime.datetime.now() - start))


if __name__ == '__main__':
    main(*sys.argv[1:])
