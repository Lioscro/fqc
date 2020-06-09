import os
from collections import namedtuple

from .config import WHITELIST_DIR


class Technology(
        namedtuple('Technology',
                   ['name', 'description', 'n_files', 'reads_file',
                    'umi_positions', 'barcode_positions', 'whitelist_path'])):

    def __str__(self):
        return str(self.name)


# If the permutation is [1, 0, 2], then read 0 is from fastq 1,
# read 1 is from fastq 0, and read 2 is from fastq 2
class OrderedTechnology(namedtuple('OrderedTechnology',
                                   ['technology', 'permutation'])):

    def __str__(self):
        return f'{self.technology}{self.permutation}'


ReadSubstring = namedtuple('ReadSubstring', ['file', 'start', 'stop'])

TECHNOLOGIES = [
    Technology(
        '10xv1',
        '10x version 1',
        3,
        ReadSubstring(2, None, None),
        [ReadSubstring(1, 0, 10)],
        [ReadSubstring(0, 0, 14)],
        os.path.join(WHITELIST_DIR, '10xv1_whitelist.txt.gz'),
    ),
    Technology(
        '10xv2',
        '10x version 2',
        2,
        ReadSubstring(1, None, None),
        [ReadSubstring(0, 16, 26)],
        [ReadSubstring(0, 0, 16)],
        os.path.join(WHITELIST_DIR, '10xv2_whitelist.txt.gz'),
    ),
    Technology(
        '10xv3',
        '10x version 3',
        2,
        ReadSubstring(1, None, None),
        [ReadSubstring(0, 16, 28)],
        [ReadSubstring(0, 0, 16)],
        os.path.join(WHITELIST_DIR, '10xv3_whitelist.txt.gz'),
    ),
    # Technology(
    #     'dropseq',
    #     'DropSeq',
    #     2,
    #     ReadSubstring(1, None, None),
    #     [ReadSubstring(0, 12, 20)],
    #     [ReadSubstring(0, 0, 12)],
    #     None,
    # ),
    # Technology(
    #     'indropsv1',
    #     'inDrops version 1',
    #     2,
    #     ReadSubstring(1, None, None),
    #     [ReadSubstring(0, 42, 48)],
    #     [ReadSubstring(0, 0, 11),
    #      ReadSubstring(0, 30, 38)],
    #     None,
    # ),
    # Technology(
    #     'indropsv2',
    #     'inDrops version 2',
    #     2,
    #     ReadSubstring(0, None, None),
    #     [ReadSubstring(1, 42, 48)],
    #     [ReadSubstring(1, 0, 11),
    #      ReadSubstring(1, 30, 38)],
    #     None,
    # ),
    Technology(
        'indropsv3',
        'inDrops version 3',
        3,
        ReadSubstring(2, None, None),
        [ReadSubstring(1, 8, 14)],
        [ReadSubstring(0, 0, 8), ReadSubstring(1, 0, 8)],
        os.path.join(WHITELIST_DIR, 'indropsv3_whitelist.txt.gz'),
    ),
]
TECHNOLOGIES_MAPPING = {t.name: t for t in TECHNOLOGIES}
