import gzip

import logging
from urllib.parse import urlparse
from urllib.request import urlopen

from tqdm import tqdm


class TqdmLoggingHandler(logging.Handler):
    """Custom logging handler so that logging does not affect progress bars.
    """

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            self.handleError(record)


def open_as_text(path, mode):
    """Open a textfile or gzip file in text mode.

    :param path: path to textfile or gzip
    :type path: str
    :param mode: mode to open the file, either `w` for write or `r` for read
    :type mode: str

    :return: file object
    :rtype: file object
    """
    return gzip.open(path,
                     f'{mode}t') if path.endswith('.gz') else open(path, mode)


def fastq_reads(path):
    """Generator for reads in a local or remote FASTQ file.

    :param path: path or url to FASTQ file
    :type path: str

    :return: generator for each read
    :rtype: generator
    """
    parse = urlparse(path)
    open_func = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'

    with open_func(urlopen(path) if parse.scheme and parse.netloc else path,
                   mode) as f:
        for n, line in enumerate(f):
            if (n + 3) % 4 == 0:
                yield line.strip()


def sequence_equals(seq1, seq2, distance=0):
    """Determine if two DNA string sequences are equal, also considering two
    sequences with <= the specified hamming distance apart as equal. Deletions
    are not allowed.

    :param seq1: first sequence
    :type seq1: str
    :param seq2: second sequence
    :type seq2: str
    :param distance: allowed hamming distance, defaults to `0`
    :type distance: int, optional

    :return: whether or not the two sequences are "equal"
    :rtype: bool
    """
    if len(seq1) != len(seq2):
        return False

    seq1 = seq1.upper()
    seq2 = seq2.upper()
    return sum(
        s1 != 'N' and s2 != 'N' and s1 != s2 for s1, s2 in zip(seq1, seq2)
    ) <= distance
