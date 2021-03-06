import argparse
import logging
import sys

from . import __version__
from .config import N_READS, SKIP_READS
from .fqc import fqc_bam, fqc_fastq

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
        help=(
            'Prefix for FASTQ files to generate. Ignored if `--split-bam` is not used.'
        ),
        type=str,
        default=''
    )
    bam_args.add_argument(
        '-t',
        metavar='THREADS',
        help='Number of threads (default: 4)',
        type=int,
        default=4
    )
    bam_args.add_argument(
        '--split-bam',
        help=(
            'Revert the BAM file into its constituent FASTQ files. '
            'The FASTQ files will be named PREFIX_i.fastq.gz if `-p` is '
            'provided, i.fastq.gz otherwise, where i is a positive read index.'
        ),
        action='store_true'
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
        result = fqc_bam(
            args.files[0], split=args.split_bam, prefix=args.p, threads=args.t
        )
        if not args.split_bam:
            logger.info((
                'Use `--split-bam` to revert the BAM file to its constituent FASTQ '
                'files, which can then be used for pre-processing.'
            ))
            logger.info(f'Detected technology: {result}')
            print(result)
            return
    elif all(file.endswith(('.fastq.gz', '.fastq')) for file in args.files):
        logger.info('Running in mode: FASTQ')
        result = fqc_fastq(args.files, args.s, args.n)

    else:
        parser.error(
            'All input files must be FASTQ (either .fastq.gz or .fastq) or a single BAM (.bam)'
        )

    fastqs, technologies = result
    if not technologies:
        logger.error('Failed to detect technology')
    elif len(technologies) == 1:
        technology = technologies[0]
        logger.info(f'Detected technology: {technology}')
        print(technology.technology)
        print(' '.join(fastqs[i] for i in technology.permutation))
    else:
        logger.warning(
            f'Ambiguous technologies {", ".join(str(technology) for technology in technologies)}'
        )
