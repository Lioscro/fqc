import os
from unittest import TestCase


class TestMixin(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.bam_10xv2_path = os.path.join(cls.fixtures_dir, '10xv2.bam')
        cls.fastq_10xv2_paths = [
            os.path.join(cls.fixtures_dir, '10xv2_1.fastq.gz'),
            os.path.join(cls.fixtures_dir, '10xv2_2.fastq.gz')
        ]
