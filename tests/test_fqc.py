from collections import OrderedDict
from unittest import mock, TestCase

import fqc.fqc as fqc
from fqc.technologies import OrderedTechnology, TECHNOLOGIES_MAPPING


class TestFqc(TestCase):

    def test_all_ordered_technologies(self):
        self.assertEqual([
            OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (0, 1)),
            OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (1, 0))
        ], fqc.all_ordered_technologies([TECHNOLOGIES_MAPPING['10xv1']], 2))

    def test_extract_barcodes_umis(self):
        reads = OrderedDict()
        reads['1'] = ['2' * 50]
        reads['2'] = ['4' * 50]

        self.assertEqual(
            ({
                '10xv1': {(0, 1): [['22222222222222']]}
            }, {
                '10xv1': {(0, 1): [['4444444444']]}
            }, {
                '10xv1': set()
            }),
            fqc.extract_barcodes_umis(
                reads,
                [OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (0, 1))]
            )
        )

    def test_extract_barcodes_umis_invalid(self):
        reads = OrderedDict()
        reads['1'] = ['2' * 5]
        reads['2'] = ['4' * 5]

        self.assertEqual(
            ({
                '10xv1': {}
            }, {
                '10xv1': {}
            }, {
                '10xv1': {(0, 1)}
            }),
            fqc.extract_barcodes_umis(
                reads,
                [OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (0, 1))]
            )
        )

    def test_filter_files(self):
        fastqs = [1, 2, 3]
        result = fqc.filter_files(
            fastqs, [
                OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (0, 1, 2)),
                OrderedTechnology(TECHNOLOGIES_MAPPING['10xv2'], (0, 1))
            ]
        )
        self.assertListEqual([
            OrderedTechnology(TECHNOLOGIES_MAPPING['10xv1'], (0, 1, 2))
        ], result)

    def test_filter_barcodes_umis(self):
        pass

    def test_fqc_fastq(self):
        with mock.patch('fqc.fqc.Fastq') as Fastq,\
            mock.patch('fqc.fqc.filter_files') as filter_files,\
            mock.patch('fqc.fqc.filter_barcodes_umis') as filter_barcodes_umis:
            Fastq.return_value = ['r1', 'r2']
            fqc.fqc_fastq(['f1', 'f2'], 0, 2)

            reads = OrderedDict()
            reads['f1'] = ['r1', 'r2']
            reads['f2'] = ['r1', 'r2']
            filter_files.assert_called_once_with(reads)
            filter_barcodes_umis.assert_called_once_with(reads, filter_files())
