from unittest import TestCase

import fqc.fastq as fastq
from tests.mixins import TestMixin


class TestFastq(TestMixin, TestCase):

    def test_getitem(self):
        f = fastq.Fastq(self.fastq_10xv2_paths[0])
        self.assertEqual(['TTCTACAGTGTGGTTTTGGACAGGTG'], f[0:1])
