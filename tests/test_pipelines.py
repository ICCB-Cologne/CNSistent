import unittest
import numpy as np
import pandas as pd

from cns.process.pipelines import main_segment, main_signatures, main_coverage, main_ploidy, main_fill, main_impute
from cns.process.segments import get_genome_segments, regions_remove, regions_select
from cns.utils.selection import only_aut
    
class TestPipelines(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 50, 100, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [50, 100, 150, 300, 400, 500, 50, 99, 100, 100, 120, 130],
            'major_cn': [1, 2, 1, 3, 4, 5, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, np.NaN, 0, 4, 3, 1, 0, 0, 1, 0, 1],
        })       
        self.samples = pd.DataFrame({
            'sample_id': ['s1', 's2', 's3', 's4'],
            'sex': ['xx', 'xy', 'xx', 'xy']
        }).set_index('sample_id')
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_names': ['chr1', 'chr2', 'chr3', 'chrX', 'chrY'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 600,
            'sex_names': ['chrX', 'chrY'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        pd.set_option('display.max_columns', 10)

    def test_main_fill(self):
        res = main_fill(self.cns, self.samples, assembly=self.assembly)
        res["length"] = res["end"] - res["start"]
        # assert that length sum is equal to aut_len for each sample
        self.assertTrue(np.allclose(only_aut(res, self.assembly).groupby("sample_id")["length"].sum(), self.assembly.aut_len))
        # assert that chrY exists in chrom column where index is s4 and not in s3
        self.assertTrue("chrY" in res.query("sample_id == 's4'")['chrom'].values)
        self.assertFalse("chrY" in res.query("sample_id == 's3'")['chrom'].values)

    def test_main_impute(self):
        res = main_fill(self.cns, self.samples, assembly=self.assembly)
        res["length"] = res["end"] - res["start"]
        res = main_impute(res, self.samples)
        res["length"] = res["end"] - res["start"]
        self.assertTrue(np.allclose(only_aut(res, self.assembly).groupby("sample_id")["length"].sum(), self.assembly.aut_len))
        # assert that chrY exists in chrom column where index is s4 and not in s3
        self.assertTrue("chrY" in res.query("sample_id == 's4'")['chrom'].values)
        self.assertFalse("chrY" in res.query("sample_id == 's3'")['chrom'].values)

    def test_main_coverage(self):
        res = main_coverage(self.cns, self.samples, assembly=self.assembly)     
        self.assertEqual(res.shape, (4, 9))
        self.assertEqual(res.loc['s1', 'chrom_missing'][-1], "chrX")
        self.assertEqual(res.loc['s1', 'chrom_count'], 1)
        # for all rows, cov_hom_aut is lower than cov_het_aut
        self.assertTrue(np.all(res['cover_hom_aut'] <= res['cover_het_aut']))        
        self.assertEqual(res.loc['s1', 'cover_het_sex'], 0)
        self.assertEqual(res.loc['s2', 'cover_het_sex'], 0.5)
        self.assertEqual(res.loc['s4', 'cover_het_aut'], 0.3)
    
    def test_main_ploidy(self):
        res = main_ploidy(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.shape, (4, 16))
        self.assertTrue(np.all(res['ane_hom_aut'] <= res['ane_het_aut']))       
        self.assertEqual(res.loc['s1', 'ane_het_sex'], 0)
        self.assertEqual(res.loc['s2', 'ane_het_sex'], 1/2)
        self.assertEqual(res.loc['s4', 'ane_het_sex'], 0)
        self.assertEqual(res.loc['s2', 'loh_het_all'], 1/8)
        self.assertEqual(res.loc['s4', 'loh_hom_aut'], 1/self.assembly.aut_len)
    
    def test_main_signatures(self):
        res = main_signatures(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.shape, (4, 28))        
        self.assertEqual(res.loc['s1', 'breaks_minor_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_major_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_total_cn_aut'], 1)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_aut'], 4)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_sex'], 0)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_all'], 4)

    def test_main_segment(self):
        select = {'chr1': [(0, 100)], 'chr2': [(50, 150)], 'chr3': [(100, 200), (250, 300)]}
        remove = {'chr2': [(0, 75)], 'chr3': [(150, 175)], 'chrX': [(0, 100)]}
        res = main_segment(self.cns, select, remove, assembly=self.assembly)
        self.assertEqual(len(res), 5)
        res = main_segment(self.cns, select, remove, filter_size=50, assembly=self.assembly)
        self.assertEqual(len(res), 4)
        res = main_segment(self.cns, select, remove, merge_dist=25, filter_size=50, assembly=self.assembly)
        self.assertEqual(res.iloc[0].tolist(), ["chr1", 0, 50])
        self.assertEqual(res.iloc[2].tolist(), ["chr2", 75, 117])        
        res = main_segment(self.cns, select, remove, 25, 25, 50, self.assembly)
        self.assertEqual(res.iloc[0].tolist(), ["chr1", 0, 25])
        self.assertEqual(res.iloc[4].tolist(), ["chr2", 75, 96])    
