from unittest import TestCase

import fqc.fqc as fqc
from fqc.technologies import TECHNOLOGIES_MAPPING


class TestFqc(TestCase):

    def test_filter_files(self):
        fastqs = [1, 2, 3]
        result = fqc.filter_files(fastqs)
        self.assertListEqual([
            TECHNOLOGIES_MAPPING['10xv1'], TECHNOLOGIES_MAPPING['indropsv3']
        ], result)

    def test_fqc(self):
        fqc.fqc([])
