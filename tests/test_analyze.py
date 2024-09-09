import unittest
import numpy as np
import pandas as pd

from cns.analyze.aneuploidy import get_expected_ploidy, get_ane_for_cols, get_ane_for_samples, get_ane_bases
from cns.analyze.coverage import normalize_feature, get_covered_bases, get_missing_chroms
from cns.analyze.breakpoints import calc_breaks_per_sample, calc_breaks_per_chr
from cns.process.binning import add_cns_loc, sum_cns
from cns.process.pipelines import main_coverage


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
        self.cns["length"] = self.cns["end"] - self.cns["start"]
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
        self.assertEqual(res.loc['s1', 'cover_het_aut'], 100)

    def test_get_base_frac(self):
        samples_df = get_covered_bases(self.cns, self.samples, True)
        res = normalize_feature(samples_df, "cover_het", assembly=self.assembly)
        self.assertEqual(res.loc['s1', 'cover_het_aut'], 1/3)
        self.assertEqual(res.loc['s2', 'cover_het_sex'], 1)
        self.assertEqual(res.loc['s4', 'cover_het_tot'], 1/4)

    def test_calculate_coverage(self):
        res = main_coverage(self.cns, self.samples, assembly=self.assembly)
        print(res)
        self.assertEqual(res['cover_het_sex']['s1'], 1/2)
        self.assertEqual(res['cover_het_tot']['s1'], 2/5)
        self.assertEqual(res['cover_hom_aut']['s2'], 0)
        self.assertEqual(res['cover_het_aut']['s2'], 1/3)
        self.assertEqual(res['chrom_count']['s2'], 2)
        self.assertEqual(res['chrom_missing']['s3'].size, 3)


class TestSignatures(unittest.TestCase):
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
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })

    def test_calc_breaks_per_chr(self):
        result = calc_breaks_per_chr(self.cns)
        self.assertEqual(result.query('sample_id == "s1" and chrom == "chr1"')['breaks'].values[0], 1)
        self.assertEqual(result.query('sample_id == "s2" and chrom == "chr2"')['breaks'].values[0], 0)
        self.assertEqual(result.query('sample_id == "s4" and chrom == "chr1"')['breaks'].values[0], 2)

    def test_add_breaks_per_sample(self):
        res = calc_breaks_per_sample(self.cns, self.samples, "major_cn", self.assembly)
        print(res)
        self.assertEqual(res.query('sample_id == "s1"')['breaks_major_cn_aut'].values[0], 1)
        self.assertEqual(res.query('sample_id == "s4"')['breaks_major_cn_tot'].values[0], 3)        
        res = calc_breaks_per_sample(self.cns, self.samples, "major_cn", self.assembly)
        print(res)


class TestAneuploidy(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'major_cn': [1, 2, 3, 4, np.nan, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, np.nan, 4, 3, 1, 0, 0, 1, 0, 1],
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
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'chr_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        self.ane_cols = ["major_cn", "minor_cn"]

    def test_get_expected_ploidy(self):
        self.assertEqual(get_expected_ploidy("minor_cn", "chrX", True), 0)
        self.assertEqual(get_expected_ploidy("major_cn", "chrX", True), 1)
        self.assertEqual(get_expected_ploidy("total_cn", "chrX", False), 2)
        self.assertEqual(get_expected_ploidy("major_cn", "chrX", False), 1)
        self.assertEqual(get_expected_ploidy("minor_cn", "chrY", True), 0)
        self.assertEqual(get_expected_ploidy("major_cn", "chrY", True), 1)
        self.assertEqual(get_expected_ploidy("total_cn", "chrY", False), 0)
        self.assertEqual(get_expected_ploidy("total_cn", "chr1", True), 2)
        self.assertEqual(get_expected_ploidy("major_cn", "chr1", True), 1)

    def test_get_ane_for_cols(self):
        res = get_ane_for_cols(self.cns, self.samples, self.ane_cols, self.assembly)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(res[0]), len(self.cns["major_cn"]))
        self.assertEqual(res[0][0], False)
        self.assertEqual(res[1][0], True)

    def test_calc_ane_per_sample(self):
        cns_df = add_cns_loc(self.cns.copy())
        res = get_ane_for_samples(cns_df, self.samples, self.ane_cols, True, self.assembly)
        self.assertEqual(len(res), 4)
        self.assertEqual(res.values[3], 170)
        res = get_ane_for_samples(cns_df, self.samples, self.ane_cols, False, self.assembly)
        self.assertEqual(res.values[3], 1)

    def test_get_ane_bases(self):
        res = get_ane_bases(self.cns, self.samples, self.ane_cols, self.assembly)
        self.assertEqual(res.shape, (4, 7))
        self.assertEqual(res.loc['s4', 'ane_het_aut'], 170)
        self.assertEqual(res.loc['s2', 'ane_het_sex'], 100)
        self.assertEqual(res.loc['s2', 'ane_het_tot'], 200)
        norm = normalize_feature(res, "ane_het", self.assembly)
        norm = normalize_feature(res, "ane_hom", self.assembly)
        print(norm)