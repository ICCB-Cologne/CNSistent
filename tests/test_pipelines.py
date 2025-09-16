import io
import unittest
import numpy as np
import pandas as pd

from cns import *

class TestPipelines(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 50, 100, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [50, 100, 150, 300, 400, 500, 50, 99, 100, 100, 120, 130],
            'major_cn': [1, 2, 1, 3, 4, 5, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, np.nan, 0, 4, 3, 1, 0, 0, 1, 0, 1],
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
        self.segs = {'chr1': [(0, 100, 0)], 'chr2': [(50, 150, 1)], 'chrY': [(0, 100, 2)], 'chr3': [(100, 200, 3), (250, 300, 4)]}
        self.aut_segs ={'chr1': [(0, 100, "0")]}
        pd.set_option('display.max_columns', 10)

    def test_main_align(self):
        res = main_align(self.cns, self.samples, assembly=self.assembly)
        # assert that length sum is equal to aut_len for each sample
        auts = only_aut(res, self.assembly)
        self.assertTrue(np.allclose(calc_lengths(auts).groupby(auts["sample_id"]).sum(),self.assembly.aut_len))
        # assert that chrY exists in chrom column where index is s4 and not in s3
        self.assertTrue("chrY" in res.query("sample_id == 's4'")['chrom'].values)
        self.assertFalse("chrY" in res.query("sample_id == 's3'")['chrom'].values)

    def test_main_impute(self):
        res = main_align(self.cns, self.samples, assembly=self.assembly)
        res = main_infer(res, self.samples)
        auts = only_aut(res, self.assembly)
        self.assertTrue(np.allclose(calc_lengths(auts).groupby(auts["sample_id"]).sum(), self.assembly.aut_len))
        # assert that chrY exists in chrom column where index is s4 and not in s3
        self.assertTrue("chrY" in res.query("sample_id == 's4'")['chrom'].values)
        self.assertFalse("chrY" in res.query("sample_id == 's3'")['chrom'].values)

    def test_main_coverage(self):
        res = main_coverage(self.cns, self.samples, assembly=self.assembly)    
        self.assertEqual(res.shape, (4, 9))
        self.assertEqual(res.loc['s1', 'chrom_missing'][-1], "chrX")
        self.assertEqual(res.loc['s1', 'chrom_count'], 1)
        # for all rows, cov_both_aut is lower than cov_any_aut
        self.assertTrue(np.all(res['cover_both_aut'] <= res['cover_any_aut']))        
        self.assertEqual(res.loc['s1', 'cover_any_sex'], 0)
        self.assertEqual(res.loc['s2', 'cover_any_sex'], 0.5)
        self.assertEqual(res.loc['s4', 'cover_any_aut'], 0.3)

    def test_main_coverage_segs(self):
        res = main_coverage(self.cns, self.samples, segs=self.segs, assembly=self.assembly)         
        self.assertEqual(res.loc['s2', 'cover_both_aut'], 0)    
        self.assertEqual(res.loc['s2', 'cover_any_aut'], 50/350)
    
    def test_main_ploidy(self):
        res_df = main_ploidy(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res_df.shape, (4, 21))
        self.assertTrue(np.all(res_df['ane_both_aut'] <= res_df['ane_any_aut']))   
        self.assertEqual(res_df.loc['s1', 'ane_any_sex'], 0)
        self.assertEqual(res_df.loc['s2', 'ane_any_sex'], 1/2)
        self.assertEqual(res_df.loc['s4', 'ane_any_sex'], 0)
        self.assertEqual(res_df.loc['s2', 'loh_any_all'], 1/8)
        self.assertEqual(res_df.loc['s2', 'imb_major_cn_aut'], 1/6)        
        self.assertEqual(res_df.loc['s4', 'imb_major_cn_sex'], 0)

    def test_main_ploidy_segs(self):
        res = main_ploidy(self.cns, self.samples, segs=self.segs, assembly=self.assembly)
        self.assertEqual(res.loc['s2', 'ane_any_aut'], 0) # check that NaN and out of scope are not considered
        self.assertEqual(res.loc['s4', 'loh_both_aut'], 1/350)
    
    def test_main_breakage(self):
        with self.assertRaises(Exception):
            main_breakage(self.cns, self.samples, assembly=self.assembly)
        self.cns.loc[2, "minor_cn"] = 0
        res = main_breakage(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.shape, (4, 19))        
        self.assertEqual(res.loc['s1', 'breaks_minor_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_major_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_total_cn_aut'], 1)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_aut'], 4)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_sex'], 0)
        self.assertEqual(res.loc['s4', 'breaks_total_cn_all'], 4)

    def test_main_breakage_segs(self):
        res = main_breakage(self.cns, self.samples, segs=self.aut_segs, assembly=self.assembly)
        self.assertEqual(res.loc['s1', 'breaks_major_cn_aut'], 1)

    def test_main_segment(self):
        select = {'chr1': [(0, 100, "0")], 'chr2': [(50, 150, "1")], 'chr3': [(100, 200, "2"), (250, 300, "3")]}
        remove = {'chr2': [(0, 75, "0")], 'chr3': [(150, 175, "1")], 'chrX': [(0, 100, "3")]}
        res = main_segment(select, remove)
        self.assertEqual(len(res), 3)
        res = main_segment(select, remove, filter_size=50)
        self.assertEqual(len(res), 3)
        res = main_segment(select, remove, merge_dist=25, filter_size=50)
        self.assertEqual(res["chr1"][0], (0, 100, "chr1_0"))
        self.assertEqual(res["chr2"][0], (75, 150, "chr2_0")) 

    def test_main_segment_cns(self):
        remove = {'chr2': [(0, 75, "0")], 'chr3': [(150, 175, "1")], 'chrX': [(0, 100, "3")]}
        cns_segs = cns_df_to_segments(self.cns, process="unify")
        res = main_segment(cns_segs, remove, 25, 25, 50)
        self.assertEqual(res["chr1"][1], (25, 50, "chr1_0_1"))
        self.assertEqual(res["chr2"][0], (75, 110, "chr2_0_0"))   

    def test_get_norm_sizes(self):
        res = get_norm_sizes(self.segs)
        self.assertEqual(res['aut'], 350)
        self.assertEqual(res['sexXX'], 0)
        self.assertEqual(res['sexXY'], 100)
        self.assertEqual(res['allXX'], 350)
        self.assertEqual(res['allXY'], 450)
        res = get_norm_sizes(None, assembly=self.assembly)
        self.assertEqual(res['aut'], 600)
        self.assertEqual(res['sexXX'], 100)
        self.assertEqual(res['sexXY'], 200)
        self.assertEqual(res['allXX'], 700)
        self.assertEqual(res['allXY'], 800)

    def test_get_chr_sets(self):
        res = get_chr_sets(self.cns, assembly=self.assembly)
        self.assertEqual(res['aut'], ['chr1', 'chr2', 'chr3'])
        self.assertEqual(res['sex'], ['chrY'])
        self.assertEqual(res['all'], ['chr1', 'chr2', 'chr3', 'chrY'])
        cns_df = aggregate_by_segments(self.cns, self.aut_segs, 'none')
        res = get_chr_sets(cns_df, assembly=self.assembly)
        self.assertTrue('aut' in res)   
        self.assertTrue('sex' not in res)   
        self.assertTrue('all' not in res)   

class TestData(unittest.TestCase):
    def setUp(self):
        cns = """
        sample_id, chrom, start, end, major_cn, minor_cn
        s1, chr19, 1000000, 3000000, 1,
        s1, chr19, 3000000, 11000000, 1, 1
        s1, chr19, 14000000, 21000000, 3, 1
        s1, chr19, 21000000, 25000000, 3, 
        s1, chr19, 28000000, 58500000, 3,
        s2, chr19, 1000000, 24000000, 2,
        s2, chr19, 29000000, 58000000, 0,
        """
        self.cns_df = pd.read_csv(io.StringIO(cns.strip()), sep=',\s*', engine='python')
        pd.set_option('display.max_columns', 10)

    def test_align(self):
        cns_align_df = main_align(self.cns_df, segs={"chr19": [(0, hg19.chr_lens["chr19"], "chr19")]})
        lens = cns_align_df["end"] - cns_align_df["start"]
        self.assertEqual(lens.sum() / 2, hg19.chr_lens["chr19"]) # 2 samples of full length        
        non_nan = cns_align_df.dropna()
        lens = non_nan["end"] - non_nan["start"]
        self.assertGreater(hg19.chr_lens["chr19"], lens.sum() / 2)

    def test_impute(self):        
        cns_align_df = main_align(self.cns_df, segs={"chr19": [(0, hg19.chr_lens["chr19"], "chr19")]})
        cns_imp_df = main_infer(cns_align_df)
        non_nan = cns_imp_df.dropna()
        lens = non_nan["end"] - non_nan["start"]
        self.assertEqual(lens.sum() / 2, hg19.chr_lens["chr19"]) # 2 samples of full length

    def test_coverage(self):
        # cns_cov_df = main_coverage(self.cns_df)
        # print(cns_cov_df)
        cns_cov_gap_df = main_coverage(self.cns_df, segs={ 'chr19': [(1000000, 5000000, "seg1")] })
        self.assertEqual(cns_cov_gap_df.loc['s1', 'cover_any_aut'], 1)
        self.assertEqual(cns_cov_gap_df.loc['s1', 'cover_both_aut'], 0.5)
        self.assertEqual(cns_cov_gap_df.loc['s1', 'cover_any_sex'], 0)	
        self.assertEqual(cns_cov_gap_df.loc['s2', 'cover_any_all'], 1)
        self.assertEqual(cns_cov_gap_df.loc['s2', 'cover_both_all'], 0.0)
