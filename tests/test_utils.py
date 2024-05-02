import unittest
from cns.utils.files import rename_columns
from cns.utils.conversions import *

import unittest
import pandas as pd
import numpy as np

class TestConversions(unittest.TestCase):
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

    def test_sampleid_to_numeric_df(self):
        cna_df = pd.DataFrame({'sample_id': ['sample1', 'sample2']})
        samples, df = sampleid_to_numeric_df(cna_df)
        self.assertEqual(samples, ['sample1', 'sample2'])
        self.assertTrue((df['sample_id'] == [0, 1]).all())

    def test_numeric_to_sampleid_df(self):
        samples = ['sample1', 'sample2']
        cna_df = pd.DataFrame({'sample_id': [0, 1]})
        df = numeric_to_sampleid_df(samples, cna_df)
        self.assertTrue((df['sample_id'] == ['sample1', 'sample2']).all())

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

    def test_rename_columns(self):
        cna_df = pd.DataFrame(np.zeros((5, 6)))
        df = rename_columns(cna_df)
        self.assertEqual(df.columns.tolist(), ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"])

    def test_segs_to_chrom_dict(self):
        segs = [('chr1', 1, 5), ('chr1', 4, 8), ('chr2', 10, 15)]
        chrom_dict = segs_to_chrom_dict(segs)
        expected =  {'chr1': [(1, 5), (4, 8)], 'chr2': [(10, 15)]}
        self.assertEqual(chrom_dict, expected)

if __name__ == '__main__':
    unittest.main()