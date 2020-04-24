import gzip
import os
import tempfile
from unittest import TestCase

import fqc.bam as bam
from fqc.technologies import OrderedTechnology, TECHNOLOGIES_MAPPING
from tests.mixins import TestMixin


class TestBAM(TestMixin, TestCase):

    def test_detect_technology(self):
        b = bam.BAM(self.bam_10xv2_path)
        self.assertEqual(TECHNOLOGIES_MAPPING['10xv2'], b.technology)

    def test_to_fastq_10x(self):
        b = bam.BAM(self.bam_10xv2_path)
        fastqs, technologies = b.to_fastq(
            os.path.join(tempfile.mkdtemp(), '10xv2')
        )

        self.assertEqual([
            OrderedTechnology(TECHNOLOGIES_MAPPING['10xv2'], (0, 1))
        ], technologies)
        for fastq1, fastq2 in zip(self.fastq_10xv2_paths, fastqs):
            with gzip.open(fastq1, 'rt') as f1, gzip.open(fastq2, 'rt') as f2:
                self.assertEqual(f1.read(), f2.read())
