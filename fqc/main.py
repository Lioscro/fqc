import argparse
import logging
import sys

from . import __version__
from .bam import BAM
from .config import N_READS, SKIP_READS
from .fqc import fqc

logger = logging.getLogger(__name__)


def main():
    """Command-line entrypoint.
    """
    # Main parser
    parser = argparse.ArgumentParser(description='fqc {}'.format(__version__))
    parser._actions[0].help = parser._actions[0].help.capitalize()

    parser.add_argument(
        'files',
        metavar='FILES',
        help='Input files (FASTQs or a single BAM)',
        nargs='+'
    )
    fastq_args = parser.add_argument_group('optional arguments for FASTQ files')
    fastq_args.add_argument(
        '-s',
        metavar='SKIP',
        help=f'Number of reads to skip at the beginning (default: {SKIP_READS})',
        type=int,
        default=SKIP_READS
    )
    fastq_args.add_argument(
        '-n',
        metavar='READS',
        help=f'Number of reads to use after skip (default: {N_READS})',
        type=int,
        default=N_READS
    )
    bam_args = parser.add_argument_group('optional arguments for BAM files')
    bam_args.add_argument(
        '-p',
        metavar='PREFIX',
        help='Prefix for FASTQ files to generate',
        type=str,
        default=''
    )
    parser.add_argument(
        '-t',
        metavar='THREADS',
        help='Number of threads (default: 4)',
        type=int,
        default=4
    )
    parser.add_argument(
        '--verbose', help='Print debugging information', action='store_true'
    )
    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)7s %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    logger.debug('Printing verbose output')
    logger.debug(args)

    if len(args.files) == 1 and args.files[0].endswith('.bam'):
        logger.info('Running in mode: BAM')
        bam = BAM(args.files[0])
        fastqs = bam.to_fastq(prefix=args.p, threads=args.t)
        fqc(fastqs, args.s, args.n)
    elif all(file.endswith(('.fastq.gz', '.fastq')) for file in args.files):
        logger.info('Running in mode: FASTQ')
        fqc(args.files, args.s, args.n)
    else:
        parser.error(
            'All input files must be FASTQ (either .fastq.gz or .fastq) or a single BAM (.bam)'
        )
