import unittest
import numpy as np
import pandas as pd

from cns.analyze import *
from cns.process import *
from cns.pipelines import main_coverage

class TestCoverage(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4'],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0],
            'minor_cn': [0, 2, np.nan, np.nan, 3, 1, 0, 0],
            'chrom': ['chr1', 'chrX', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99],
            'end': [100, 200, 300, 400, 500, 50, 99, 100],
        })
        self.samples = pd.DataFrame({
            'sex': ['xy', 'NA', 'xx', 'NA']
        }, index=['s1', 's2', 's3', 's4'])
        self.samples.index.name = "sample_id"
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        pd.set_option('display.max_columns', 10)

    def test_get_missing_chroms(self):
        res = get_missing_chroms(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.loc['s1', 'chrom_count'], 2)

    def test_get_covered_bases(self):
        res = get_covered_bases(self.cns, self.samples, True)
        self.assertEqual(res.loc['s1', 'cover_any_aut'], 100)

    def test_get_base_frac(self):
        samples_df = get_covered_bases(self.cns, self.samples, True)
        norm_sizes = get_norm_sizes(None, self.assembly)
        res = normalize_feature(samples_df, "cover_any", norm_sizes)
        self.assertEqual(res.loc['s1', 'cover_any_aut'], 1/3)
        self.assertEqual(res.loc['s2', 'cover_any_sex'], 1)
        self.assertEqual(res.loc['s4', 'cover_any_all'], 1/4)

    def test_calculate_coverage(self):
        res = main_coverage(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res['cover_any_sex']['s1'], 1/2)
        self.assertEqual(res['cover_any_all']['s1'], 2/5)
        self.assertEqual(res['cover_both_aut']['s2'], 0)
        self.assertEqual(res['cover_any_aut']['s2'], 1/3)
        self.assertEqual(res['chrom_count']['s2'], 2)
        self.assertEqual(res['chrom_missing']['s3'].size, 3)


class TestBreakage(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, 0, 4, 3, 1, 0, 0, 1, 0, 1],
            'chrom': ['chr1', 'chr1', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [100, 200, 300, 400, 500, 50, 99, 100, 100, 120, 130],
        })       
        self.samples = pd.DataFrame({
            'sample_id': ['s1', 's2', 's3', 's4'],
            'sex': ['xx', 'xy', 'xx', 'xy']
        }).set_index('sample_id')
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_names': ['chr1', 'chr2', 'chr3', 'chrX', 'chrY'],
            'sex_names': ['chrX', 'chrY'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 600,
            'chr_x': 'chrX',
            'chr_y': 'chrY',
            'chr_count': 5,
            'aut_count': 3,
            'sex_count': 2
        })
        pd.set_option('display.max_columns', 10)

        self.prepare_segments = lambda cns_df, cn_col: merge_cns_df(
            cns_df[["sample_id", "chrom", "start", "end", cn_col]], cn_col, False
        )

    def test_calc_breaks_per_chr(self):
        segs_df = self.prepare_segments(self.cns, "major_cn")
        result = calc_breaks_per_chr(segs_df)
        self.assertEqual(result.query('sample_id == "s1" and chrom == "chr1"')['breaks'].values[0], 1)
        self.assertEqual(result.query('sample_id == "s2" and chrom == "chr2"')['breaks'].values[0], 0)
        self.assertEqual(result.query('sample_id == "s4" and chrom == "chr1"')['breaks'].values[0], 2)

    def test_calc_breaks_per_sample(self):
        segs_df = self.prepare_segments(self.cns, "major_cn")
        res = calc_breaks_per_sample(segs_df, self.samples, "major_cn", self.assembly)
        self.assertEqual(res.query('sample_id == "s1"')['breaks_major_cn_aut'].values[0], 1)
        self.assertEqual(res.query('sample_id == "s4"')['breaks_major_cn_all'].values[0], 3)      

    def test_calc_step_per_chr(self):
        segs_df = self.prepare_segments(self.cns, "major_cn")
        res = calc_step_per_chr(segs_df, "major_cn")
        self.assertEqual(res.query('sample_id == "s1" and chrom == "chr1"')['step'].values[0], 1)
        self.assertEqual(res.query('sample_id == "s1" and chrom == "chr1"')['count'].values[0], 1)

    def test_calc_step_per_sample(self):
        segs_df = self.prepare_segments(self.cns, "major_cn")
        res = calc_step_per_sample(segs_df, self.samples, "major_cn", self.assembly)
        self.assertEqual(res.loc['s1', 'step_major_cn_aut'], 1)    
        segs_df = self.prepare_segments(self.cns, "minor_cn") 
        res = calc_step_per_sample(segs_df, self.samples, "minor_cn", self.assembly)
        self.assertEqual(res.loc['s1', 'step_minor_cn_aut'], 2)     


class TestAneuploidy(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'major_cn': [1, 2, 3, 4, np.nan, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, np.nan, 4, 3, 1, 0, 0, 1, 0, 1],
            'chrom': ['chr1', 'chr1', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 100, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [100, 200, 300, 400, 500, 50, 99, 100, 100, 120, 130],
        }).fillna(0)  
        self.samples = pd.DataFrame({
            'sample_id': ['s1', 's2', 's3', 's4'],
            'sex': ['xx', 'xy', 'xx', 'xy']
        }).set_index('sample_id')
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_names': ['chr1', 'chr2', 'chr3', 'chrX', 'chrY'],
            'sex_names': ['chrX', 'chrY'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'chr_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 600,
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        self.cols = ["major_cn", "minor_cn"]        
        pd.set_option('display.max_columns', 10)

    def test_get_ane_bases(self):
        res_df = calc_ane_bases(self.samples, self.cns, self.cols, "both", self.assembly)
        res_df = calc_ane_bases(res_df, self.cns, self.cols, "any", self.assembly)
        norm_sizes = get_norm_sizes(None, self.assembly)
        self.assertEqual(res_df.shape, (4, 7))
        self.assertEqual(res_df.loc['s4', 'ane_any_aut'], 170)
        self.assertEqual(res_df.loc['s2', 'ane_any_sex'], 100)
        self.assertEqual(res_df.loc['s2', 'ane_any_all'], 200)
        res_df = normalize_feature(res_df, "ane_any", norm_sizes)
        self.assertEqual(res_df.loc['s2', 'ane_any_sex'], 1/2)
        res_df = normalize_feature(res_df, "ane_both", norm_sizes)
        self.assertEqual(res_df.loc['s1', 'ane_both_aut'], 1/6)

    def test_get_loh_bases(self):
        res_df = calc_loh_bases(self.samples, self.cns, self.cols, "both", self.assembly)
        res_df = calc_loh_bases(res_df, self.cns, self.cols, "any", self.assembly)
        self.assertEqual(res_df.shape, (4, 7))
        self.assertEqual(res_df.loc['s1', 'loh_any_aut'], 100)
        self.assertEqual(res_df.loc['s4', 'loh_any_aut'], 70)
        self.assertEqual(res_df.loc['s4', 'loh_both_all'], 1)
        
    def test_imb_score(self):
        res = calc_imb_bases(self.cns, self.samples, self.cols, 0, self.assembly)
        self.assertEqual(res.shape, (4, 4))
        self.assertEqual(res.loc['s1', 'imb_major_cn_aut'], 100)
        self.assertEqual(res.loc['s2', 'imb_major_cn_sex'], 0)
        self.assertEqual(res.loc['s4', 'imb_major_cn_all'], 169)

class TestDistance(unittest.TestCase):
    def setUp(self):
        self.cns = {
            'sample_id': ['s1', 's1', 's2', 's2'],
            'major_cn': [1, 1, 0, 2],
            'minor_cn': [0, 1, 1, 0],
            'chrom': ['chr1', 'chr1', 'chr1', 'chr1'],	
            'start': [0, 100, 0, 100],
            'end': [100, 200, 100, 200],
            'name': ["seg_1", "seg_2", "seg_1", "seg_2"]
        }
    
    def test_sample_dist(self):
        cns_df = pd.DataFrame(self.cns)
        res = calc_distances(cns_df, 'major_cn')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], 1)
        res = calc_distances(cns_df, 'major_cn', 'euclidean')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], np.sqrt(2)/2)
        res = calc_distances(cns_df, 'major_cn', 'wasserstein')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], 0.5)

        res = calc_distances(cns_df, 'minor_cn')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], 2)
        res = calc_distances(cns_df, 'minor_cn', 'euclidean')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], np.sqrt(2))
        res = calc_distances(cns_df, 'minor_cn', 'wasserstein')
        self.assertEqual(res.loc['s1', 's1'], 0)
        self.assertEqual(res.loc['s2', 's1'], 1)

    
