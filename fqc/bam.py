import logging

import pysam
from tqdm import tqdm

from .technologies import TECHNOLOGIES
from .utils import open_as_text

logger = logging.getLogger(__name__)


def extract_10x(alignment):
    """Given a 10x BAM entry as a pysam.AlignedSegment, extract the barcode,
    UMI, and sequence.

    :param alignment: a single BAM entry
    :type alignment: pysam.AlignedSegment

    :return: a 3-tuple containing (barcodes, umis, sequence)
    :rtype: tuple
    """
    barcode = alignment.get_tag('CR')
    umi = alignment.get_tag('UR')
    # The quality scores for the sequence is returned as an array.array of
    # unsigned chars.
    sequence = alignment.query_sequence

    return ([barcode], [umi], sequence)


class BAM:
    """Class to work with BAM files.

    :param path: path to BAM file
    :type path: str
    """
    # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
    TAGS_10X = (
        'CR',  # (uncorrected) barcode
        'CY',  # barcode quality score
        'UR',  # (uncorrected) UMI
        'UY',  # UMI quality score
    )
    EXTRACT_FUNCTIONS = {
        '10xv1': extract_10x,
        '10xv2': extract_10x,
        '10xv3': extract_10x,
    }

    def __init__(self, path):
        self.path = path
        self.technology = self.detect_technology()
        logger.info(f'Detected technology {self.technology.name}')

    def detect_technology(self):
        """Detect what technology was used to generate this BAM.

        :return: a Technology object
        :rtype: Technology
        """
        with pysam.AlignmentFile(self.path, 'rb') as f:
            # Check first read of file to see the headers.
            for item in f.fetch(until_eof=True):
                # This is a 10x BAM
                if all(item.has_tag(tag) for tag in BAM.TAGS_10X):
                    # Construct a dictionary of 10x technologies for fast lookup.
                    technologies_10x = {(
                        sum(
                            substring.stop - substring.start
                            for substring in t.barcode_positions
                        ),
                        sum(
                            substring.stop - substring.start
                            for substring in t.umi_positions
                        )
                    ): t
                                        for t in TECHNOLOGIES
                                        if t.name.startswith('10x')}

                    # Extract barcode and UMI lengths
                    barcode_length = len(item.get_tag('CR'))
                    umi_length = len(item.get_tag('UR'))
                    key = (barcode_length, umi_length)

                    if key not in technologies_10x:
                        raise Exception((
                            'There is no 10x technology with barcode length '
                            f'{barcode_length} and UMI length {umi_length}'
                        ))

                    return technologies_10x[key]
                break

        raise Exception(f'Failed to detect technology for BAM {self.path}')

    def to_fastq(self, prefix='', threads=1):
        """Split the BAM into FASTQs.

        :param path: path to BAM file
        :type path: str
        :param prefix: prefix to output FASTQ files, defaults to empty string
        :type prefix: str, optional
        :param threads: number of threads to use to read the BAM file, defaults to `1`
        :type threads: int, optional

        :return: list of paths to generated FASTQs
        :rtype: list
        """
        # Write to uncompressed FASTQ for speed
        fastqs = [
            f'{prefix}_{i+1}.fastq.gz' if prefix else f'{i+1}.fastq.gz'
            for i in range(self.technology.n_files)
        ]
        logger.info(f'Splitting BAM file into FASTQs {", ".join(fastqs)}')
        logger.warning('All quality scores will be converted to F')
        files = []
        lengths = [0, 0, 0]
        for substring in self.technology.barcode_positions + self.technology.umi_positions:
            lengths[substring.file
                    ] = max(lengths[substring.file], substring.stop)

        try:
            for fastq in fastqs:
                files.append(open_as_text(fastq, 'w'))

            with pysam.AlignmentFile(self.path, 'rb', threads=threads) as f:
                count = f.count(until_eof=True)
            logger.info(f'Detected {count} BAM entries')

            with pysam.AlignmentFile(self.path, 'rb', threads=threads) as f,\
                tqdm(total=count) as pbar:
                for item in f.fetch(until_eof=True):
                    reads = ['N' * l for l in lengths]
                    barcodes, umis, sequence = BAM.EXTRACT_FUNCTIONS[
                        self.technology.name](item)  # noqa

                    # Set sequence.
                    reads[self.technology.reads_file.file] = sequence

                    # Barcode and UMI
                    for barcode, substring in zip(
                            barcodes, self.technology.barcode_positions):
                        bc = reads[substring.file]
                        reads[
                            substring.file
                        ] = f'{bc[:substring.start]}{barcode}{bc[substring.stop:]}'
                    for umi, substring in zip(umis,
                                              self.technology.umi_positions):
                        u = reads[substring.file]
                        reads[
                            substring.file
                        ] = f'{u[:substring.start]}{umi}{u[substring.stop:]}'

                    # Write to each file.
                    for file, read in zip(files, reads):
                        file.write(f'@{item.query_name}\n')
                        file.write(f'{read.upper()}\n')
                        file.write(f'+\n')
                        file.write(f'{"F" * len(read)}\n')

                    pbar.update(1)

        finally:
            for file in files:
                file.close()

        return fastqs
