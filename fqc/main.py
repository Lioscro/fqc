import argparse
import logging
import sys

from . import __version__
from .config import N_READS, SKIP_READS
from .fqc import fqc

logger = logging.getLogger(__name__)


def main():
    """Command-line entrypoint.
    """
    # Main parser
    parser = argparse.ArgumentParser(description='fqc {}'.format(__version__))
    parser._actions[0].help = parser._actions[0].help.capitalize()

    parser.add_argument('fastqs', help='FASTQ files', nargs='+')
    parser.add_argument(
        '-s', help='Number of reads to skip', type=int, default=SKIP_READS
    )
    parser.add_argument(
        '-n',
        help='Number of reads to use after skip',
        type=int,
        default=N_READS
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
    fqc(args.fastqs, args.s, args.n)
