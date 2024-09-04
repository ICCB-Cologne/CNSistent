import unittest
import numpy as np
import pandas as pd

from cns.analyze.aneuploidy import calc_ane_per_chrom, calc_aut_aneuploidy, calc_sex_aneuploidy, get_expected_ploidy, norm_aut_aneuploidy, norm_sex_aneuploidy
from cns.analyze.coverage import get_base_frac, get_covered_bases, get_missing_chroms
from cns.analyze.signatures import add_breaks_per_sample, calc_breaks_per_chr
from cns.process.binning import add_cns_loc, sum_cns
from cns.process.pipelines import main_coverage


class TestCoverage(unittest.TestCase):
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
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300
        })

    def test_get_missing_chroms(self):
        result = get_missing_chroms(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(result.loc['s1', 'chrom_count'], 2)

    def test_get_covered_bases(self):
        result = get_covered_bases(self.cns, self.samples)
        self.assertEqual(result.loc['s1', 'cover_bases_aut'], 100)

    def test_get_base_frac(self):
        samples = get_covered_bases(self.cns, self.samples)
        result = get_base_frac(samples, assembly=self.assembly)
        self.assertEqual(result.loc['s1', 'cover_frac_aut'], 1/3)
        self.assertEqual(result.loc['s2', 'cover_frac_sex'], 1)
        self.assertEqual(result.loc['s4', 'cover_frac_tot'], 1/4)

    def test_calculate_coverage(self):
        result = main_coverage(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(result['cover_bases_aut']['s1'], 100)
        self.assertEqual(result['cover_bases_sex']['s1'], 100)
        self.assertEqual(result['cover_bases_tot']['s1'], 200)
        self.assertEqual(result['chrom_count']['s2'], 2)
        self.assertEqual(result['chrom_missing']['s3'].size, 3)


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
            'sex_names': ['chrX', 'chrY']
        })

    def test_calc_breaks_per_chr(self):
        result = calc_breaks_per_chr(self.cns)
        self.assertEqual(result.query('sample_id == "s1" and chrom == "chr1"')['breaks'].values[0], 1)
        self.assertEqual(result.query('sample_id == "s2" and chrom == "chr2"')['breaks'].values[0], 0)
        self.assertEqual(result.query('sample_id == "s4" and chrom == "chr1"')['breaks'].values[0], 2)

    def test_add_breaks_per_sample(self):
        result = add_breaks_per_sample(self.cns, self.samples, self.assembly)
        self.assertEqual(result.query('sample_id == "s1"')['breaks_aut'].values[0], 1)
        self.assertEqual(result.query('sample_id == "s4"')['breaks_tot'].values[0], 4)


class TestAneuploidy(unittest.TestCase):
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
            'chr_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY']
        })
        self.ane_cols = ["ane_major_cn", "ane_minor_cn", "ane_total_cn"]

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

    def test_calc_ane_per_chrom(self):
        cns_df = sum_cns(add_cns_loc(self.cns, self.assembly))
        samples = self.samples
        res = calc_ane_per_chrom(cns_df, samples, "major_cn")
        self.assertEqual(len(res), 6)
        self.assertTrue("ane_major_cn" in res.columns)
        test_row = res.query('sample_id == "s4" and chrom == "chr1"')
        self.assertEqual(test_row['ane_major_cn'].values[0], 51)
        
        res = calc_ane_per_chrom(cns_df, samples, "minor_cn")
        test_row = res.query('sample_id == "s4" and chrom == "chr1"')
        self.assertEqual(test_row['ane_minor_cn'].values[0], 50)

    def test_calc_ane_per_sample(self):
        cns_df = sum_cns(add_cns_loc(self.cns, self.assembly))
        autosomes_sum = calc_aut_aneuploidy(cns_df, self.samples, assembly=self.assembly)
        print(autosomes_sum)
        sex_chrom_sum = calc_sex_aneuploidy(cns_df, self.samples, assembly=self.assembly)
        self.assertEqual(len(autosomes_sum), 4)
        test_row = autosomes_sum.query('sample_id == "s4"')
        self.assertEqual(test_row['ane_major_cn'].values[0], 101)
        self.assertEqual(test_row['ane_minor_cn'].values[0], 70)
        self.assertEqual(test_row['ane_total_cn'].values[0], 170)
        self.assertEqual(len(sex_chrom_sum), 4)

    def test_norm_aut_aneuploidy(self):
        cns_df = sum_cns(add_cns_loc(self.cns, self.assembly))
        autosomes_sum = calc_aut_aneuploidy(cns_df, self.samples, assembly=self.assembly)
        result = norm_aut_aneuploidy(autosomes_sum, assembly=self.assembly)
        # TODO: add more
        self.assertEqual(len(result.columns), 6)

    def test_norm_sex_aneuploidy(self):
        cns_df = sum_cns(add_cns_loc(self.cns, self.assembly))
        sex_chrom_sum = calc_sex_aneuploidy(cns_df, self.samples, assembly=self.assembly)
        result = norm_sex_aneuploidy(self.samples, sex_chrom_sum, assembly=self.assembly)
        # TODO: add more
        self.assertEqual(len(result.columns), 6)
