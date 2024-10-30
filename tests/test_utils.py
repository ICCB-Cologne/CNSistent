import unittest
import pandas as pd
import numpy as np
from cns.utils.canonization import canonize_cns_df, get_cn_cols, rename_cn_cols
from cns.utils.conversions import *
from cns.utils.files import *
from cns.utils.anomaly import *


class TestConversions(unittest.TestCase):
    def test_tuples_to_segments(self):
        gaps = [('chr1', 1, 5, 'deletion', False), ('chr2', 10, 15, 'duplication', True)]
        exp =  {'chr1': [(1, 5, 'deletion')], 'chr2': [(10, 15, 'duplication')]}
        self.assertEqual(tuples_to_segments(gaps), exp)

    def test_segs_to_df(self):
        segments = {'chr1': [(1, 5)], 'chr2': [(10, 15)]}
        df = segs_to_df(segments)
        self.assertTrue(isinstance(df, pd.DataFrame))
        self.assertEqual(df.columns.tolist(), ["chrom", "start", "end"])
        self.assertEqual(df.iloc[0].tolist(), ['chr1', 1, 5])

    def test_df_to_segs(self):
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'start': [0, 10, 20],
            'end': [5, 15, 25],
            'name': ['seg1', 'seg2', 'seg3']
        })
        exp = {'chr1': [(0, 5, 'seg1'), (10, 15, 'seg2')], 'chr2': [(20, 25, 'seg3')]}
        self.assertEqual(df_to_segs(df), exp)

    def test_genome_to_segments(self):
        assembly = type('Assembly', (object,), {
            'chr_lens':{'chr1': 100, 'chr2': 200 }
        })
        exp = {'chr1': [(0, 100, 'chr1')], 'chr2': [(0, 200, 'chr2')]}
        self.assertEqual(genome_to_segments(assembly), exp)

    def test_breaks_to_segments(self):
        breakpoints = {'chr1': [0, 50, 100], 'chr2': [0, 125, 150, 176]}
        exp = {'chr1': [(0, 50, 0), (50, 100, 1)], 'chr2': [(0, 125, 2), (125, 150, 3), (150, 176, 4)]}
        self.assertEqual(breaks_to_segments(breakpoints), exp)

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

    def test_segs_to_chrom_dict(self):
        segs = [('chr1', 1, 5), ('chr1', 4, 8), ('chr2', 10, 15)]
        chrom_dict = segs_to_chrom_dict(segs)
        expected =  {'chr1': [(1, 5), (4, 8)], 'chr2': [(10, 15)]}
        self.assertEqual(chrom_dict, expected)


class TestFiles(unittest.TestCase):
    def setUp(self):
        self.cns_df = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4'],
            'chrom': ['chrY', 'chrX', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99],
            'end': [100, 200, 300, 400, 500, 50, 99, 100],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0],
            'minor_cn': [0, 2, 0, 4, 3, 1, 0, 0],
        })        
        self.samples = pd.DataFrame({
            'sex': ['xy', 'NA', 'xx', 'NA']
        }, index=['s1', 's2', 's3', 's4'])
        self.samples.index.name = "sample_id"

    def test_fill_sex_if_missing(self):
        result = fill_sex_if_missing(self.cns_df, self.samples)
        self.assertEqual(result.loc['s1', 'sex'], 'xy')
        self.assertEqual(result.loc['s2', 'sex'], 'xy')
        self.assertEqual(result.loc['s4', 'sex'], 'xx')

    def test_chry_col(self):
        res = find_y_column(self.cns_df, self.samples, ["major_cn", "minor_cn"])
        self.assertEqual(res.loc['s1', 'y_col'], 'major_cn')
        self.assertEqual(res.loc['s2', 'y_col'], 'NA')


class TestCanonization(unittest.TestCase):
    def setUp(self):
        self.cns_df = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4'],
            'chrom': ['chr1', 'chrX', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99],
            'end': [100, 200, 300, 400, 500, 50, 99, 100],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0],
            'minor_cn': [0, 2, 0, 4, 3, 1, 0, 0],
        })        
        self.samples = pd.DataFrame({
            'sex': ['xy', 'NA', 'xx', 'NA']
        }, index=['s1', 's2', 's3', 's4'])
        self.samples.index.name = "sample_id"

    def test_canonize_cns_df(self):                      
        cns_df = canonize_cns_df(self.cns_df, print_info=False)
        self.assertEqual(cns_df.columns.tolist(), ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"])
        pd.testing.assert_frame_equal(cns_df, self.cns_df, check_dtype=False)

        cns_df = self.cns_df.rename(columns={"major_cn": "cn_a", "minor_cn": "cn_b", "chrom": "Chr", "start": "Start", "end": "End", "sample_id": ""})
        cns_df = canonize_cns_df(cns_df, order_columns=True, print_info=False)
        pd.testing.assert_frame_equal(cns_df, self.cns_df, check_dtype=False)

        cns_df = self.cns_df.rename(columns={"major_cn": "cn_a", "minor_cn": "cn_b", "chrom": "Chr", "start": "Start", "end": "End", "sample_id": ""})
        cns_df["sumCn"] = cns_df["cn_a"] + cns_df["cn_b"]
        cns_df = canonize_cns_df(cns_df, cn_columns=["sumCn"], print_info=False)
        self.assertEqual(cns_df.columns.tolist(), ["sample_id", "chrom", "start", "end", "total_cn"])
        self.assertEqual(cns_df.shape, (8, 5))

    def test_find_cn_cols_if_none(self):
        cols = get_cn_cols(self.cns_df, None)
        self.assertEqual(cols, ["major_cn", "minor_cn"])

    def test_rename_cn_cols(self):
        cns, cols = rename_cn_cols(self.cns_df.copy())
        self.assertEqual(cns.columns.tolist(), ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"])
        
        cns = self.cns_df.copy().rename(columns={"major_cn": "cn_a"})
        cns, cols = rename_cn_cols(cns)
        self.assertEqual(cols, ["major_cn", "minor_cn"])
                
        cns = self.cns_df.copy().rename(columns={"major_cn": "cn_X", "minor_cn": "cn_Y"})
        cns.loc[3, "cn_X"] = 0
        cns.loc[3, "cn_Y"] = 1
        cns, cols = rename_cn_cols(cns)
        self.assertEqual(cols, ["hap1_cn", "hap2_cn"])

        cns = self.cns_df.copy().rename(columns={"major_cn": "cn_X", "minor_cn": "cn_Y"})
        cns = cns.query("chrom != 'chrY'").copy()
        cns.loc[1, "cn_X"] = 0
        cns, cols = rename_cn_cols(cns)
        self.assertEqual(cols, ["hap1_cn", "hap2_cn"])


class TestCutoff(unittest.TestCase):
    def test_count_below_lim(self):
        vals = np.array([1, 2, 3, 4, 5])
        x, y = count_below_lim(vals, min_val=0, max_val=5, steps=5)
        self.assertTrue(np.array_equal(x, np.array([0., 1., 2., 3., 4., 5.])))
        self.assertTrue(np.array_equal(y, np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))

    def test_calculate_signed_angle(self):
        angle = calculate_signed_angle(0, 1)
        self.assertEqual(angle, 45)
        angle = calculate_signed_angle(1, 0)
        self.assertEqual(angle, -45)

    def test_find_knee(self):
        test_y = [0, 2, 4, 5, 6, 7, 8, 9, 10, 12, 14]
        test_x = list(range(len(test_y)))
        knee, _ = find_knee(test_x, test_y, knee=True)
        self.assertEqual(knee, 2)
        elbow, _ = find_knee(test_x, test_y, knee=False)
        self.assertEqual(elbow, 8)

    def test_find_bends(self):
        test_vals = [0.1, 0.2, .5, .7, .8]
        res = find_bends(test_vals, steps=5)
        self.assertEqual(len(res), 6)
        self.assertEqual(len(res[0]), 6)
        self.assertEqual(len(res[1]), 6)
        self.assertEqual(res[2], 4)
        self.assertEqual(res[4], 3)
    

if __name__ == '__main__':
    unittest.main()
