import unittest
import pandas as pd
import numpy as np

from cns.utils.conversions import *
from cns.utils.files import *


class TestConversions(unittest.TestCase):
    def test_genome_to_segments(self):
        assembly = type('Assembly', (object,), {
            'chr_lens':{'chr1': 100, 'chr2': 200 }
        })
        exp = [('chr1', 0, 100), ('chr2', 0, 200)]
        self.assertEqual(genome_to_segments(assembly), exp)

    def test_breaks_to_segments(self):
        breakpoints = [('chr1', [0, 50, 100]), ('chr2', [0, 125, 150, 176])]
        exp = [('chr1', 0, 50), ('chr1', 50, 100), ('chr2', 0, 125), ('chr2', 125, 150), ('chr2', 150, 176)]
        self.assertEqual(breaks_to_segments(breakpoints), exp)

    def test_regions_to_segments(self):
        regions = pd.DataFrame({
            'chrom': ['chr1', 'chr2'],
            'start': [0, 50],
            'end': [100, 150]
        })
        exp = [('chr1', 0, 100), ('chr2', 50, 150)]
        self.assertEqual(cns_to_segments(regions), exp)

    def test_tuples_to_segments(self):
        gaps = [('chr1', 1, 5, 'deletion', False), ('chr2', 10, 15, 'duplication', True)]
        exp = [('chr1', 1, 5), ('chr2', 10, 15)]
        self.assertEqual(tuples_to_segments(gaps), exp)

    def test_cytobands_to_df(self):
        cytobands = [['chr1', 0, 2300000, 'p36.33', 'gneg']]
        df = cytobands_to_df(cytobands)
        self.assertTrue(isinstance(df, pd.DataFrame))
        self.assertEqual(df.columns.tolist(), ["chrom", "start", "end", "name", "stain"])
        self.assertEqual(df.values.tolist(), cytobands)

    def test_gaps_to_df(self):
        gaps = [['chr1', 0, 2300000, 'p36.33', 'gneg']]
        df = gaps_to_df(gaps)
        self.assertTrue(isinstance(df, pd.DataFrame))
        self.assertEqual(df.columns.tolist(), ["chrom", "start", "end", "type", "bridge"])
        self.assertEqual(df.values.tolist(), gaps)

    def test_chrom_to_sortable(self):
        self.assertEqual(chrom_to_sortable('chr1'), 1)
        self.assertEqual(chrom_to_sortable('chrX'), 23)
        self.assertEqual(chrom_to_sortable('chrY'), 24)
        self.assertEqual(chrom_to_sortable('chrM'), 25)

    def test_sortable_to_chrom(self):
        self.assertEqual(sortable_to_chrom(1), 'chr1')
        self.assertEqual(sortable_to_chrom(23), 'chrX')
        self.assertEqual(sortable_to_chrom(24), 'chrY')
        self.assertEqual(sortable_to_chrom(25), 'chrM')

    def test_column_to_label(self):
        self.assertEqual(column_to_label('total_cn'), 'Total CN')
        self.assertEqual(column_to_label('major_cn'), 'Major CN')
        self.assertEqual(column_to_label('minor_cn'), 'Minor CN')
        self.assertEqual(column_to_label('other'), 'other')

    def test_canonize_cns_df(self):
        cns_df = pd.DataFrame([["sample1", "chr1", 0, 20, 0,  0, 1], ["sample1", "chr1", 0, 20, 0, 0, 1]])
        print(cns_df)
        cns_df.columns = ["", "", "", "", "foo", "major_cn", "minor_cn"]
        print(cns_df)                            
        cns_df, cols = canonize_cns_df(cns_df)
        self.assertEqual(cns_df.columns.tolist(), ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"])
        self.assertEqual(cols, ["major_cn", "minor_cn"])

    def test_segs_to_chrom_dict(self):
        segs = [('chr1', 1, 5), ('chr1', 4, 8), ('chr2', 10, 15)]
        chrom_dict = segs_to_chrom_dict(segs)
        expected =  {'chr1': [(1, 5), (4, 8)], 'chr2': [(10, 15)]}
        self.assertEqual(chrom_dict, expected)


class TestFiles(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4'],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0],
            'minor_cn': [0, 2, 0, 4, 3, 1, 0, 0],
            'chrom': ['chr1', 'chrX', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99],
            'end': [100, 200, 300, 400, 500, 50, 99, 100],
        })        
        self.samples = pd.DataFrame({
            'sex': ['xy', 'NA', 'xx', 'NA']
        }, index=['s1', 's2', 's3', 's4'])
        self.samples.index.name = "sample_id"

    def test_fill_sex_if_missing(self):
        result = fill_sex_if_missing(self.cns, self.samples)
        self.assertEqual(result.loc['s1', 'sex'], 'xy')
        self.assertEqual(result.loc['s2', 'sex'], 'xy')
        self.assertEqual(result.loc['s4', 'sex'], 'xx')


if __name__ == '__main__':
    unittest.main()