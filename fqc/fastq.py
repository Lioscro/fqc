import gzip
from urllib.parse import urlparse
from urllib.request import urlopen


class Fastq:
    """Class that represents a single FASTQ file.

    :param path: path to FASTQ file, may be remote
    :type path: str
    """

    def __init__(self, path):
        self.path = path

    def open(self, mode='r'):
        open_func = gzip.open if self.path.endswith('.gz') else open
        mode = f'{mode}t' if self.path.endswith('.gz') else mode
        parse = urlparse(self.path)
        return open_func(
            urlopen(self.path) if parse.scheme else self.path, mode
        )

    def __getitem__(self, index):
        if not isinstance(index, slice):
            raise NotImplementedError('Indexing is not supported. Use slice.')
        if index.stop is None or index.start < 0 or index.stop < 0 or (
                index.step is not None and index.step < 0):
            raise NotImplementedError(
                'Slices must only contain non-negative integers.'
            )
        indices = set(range(index.start or 0, index.stop, index.step or 1))

        reads = []
        with self.open('r') as f:
            for l, line in enumerate(f):  # noqa
                i = (l - 1) / 4
                if i > index.stop - 1:
                    break

                if i in indices:
                    reads.append(line.strip())
        return reads
