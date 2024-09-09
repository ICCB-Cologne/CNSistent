import unittest
import numpy as np
import pandas as pd

from cns.process.pipelines import get_genome_segments, main_signatures, regions_remove, regions_select, main_coverage, main_ploidy, main_fill, main_impute
    
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
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        pd.set_option('display.max_columns', 10)

    def test_regions_remove(self):
        bin_size = 10000000
        filter_size = 0
        select = regions_select("")
        remove = regions_remove("gaps")
        self.assertGreater(len(remove), 0)
        segs = get_genome_segments(select, bin_size, remove, filter_size)
        self.assertGreater(len(segs), 0)
        self.assertEqual(remove[0][2], segs[0][1]) # check if the first segment is a gap
        
    def test_get_genome_segments(self):
        select = [(1, 0, 10), (1, 20, 30), (2, 0, 5)]
        remove = [(1, 5, 15)]
        
        filter_size = 1
        expected_result = [(1, 15, 20), (1, 20, 30)]        
        result = get_genome_segments(select, filter_size=filter_size, remove=remove)
        
        filter_size = 6
        expected_result = [(1, 20, 30)]        
        result = get_genome_segments(select, filter_size=filter_size, remove=remove)
        
        self.assertEqual(result, expected_result)

    def test_main_fill(self):
        res = main_fill(self.cns, self.samples, assembly=self.assembly)
        print(res)
        raise NotImplementedError()

    def test_main_impute(self):
        res = main_fill(self.cns, self.samples, assembly=self.assembly)
        res = main_impute(res, self.samples)
        print(res)
        raise NotImplementedError()

    def test_main_coverage(self):
        res = main_coverage(self.cns, self.samples, assembly=self.assembly)      
        self.assertEqual(res.shape, (4, 9))
        self.assertEqual(res.loc['s1', 'chrom_missing'][-1], "chrX")
        self.assertEqual(res.loc['s1', 'chrom_count'], 1)
        self.assertEqual(res.loc['s1', 'cover_hom_sex'], 0.0)
        self.assertEqual(res.loc['s1', 'cover_het_tot'], 0.25)
    
    def test_main_ploidy(self):
        res = main_ploidy(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.shape, (4, 7))
        self.assertEqual(res.loc['s1', 'ane_het_sex'], 0)
        self.assertEqual(res.loc['s1', 'ane_het_tot'], 0.25)
        self.assertEqual(res.loc['s4', 'ane_het_tot'], 0.34)
        self.assertEqual(res.loc['s2', 'ane_hom_tot'], 0.4)
        print(res)
    
    def test_main_signatures(self):
        res = main_signatures(self.cns, self.samples, assembly=self.assembly)
        self.assertEqual(res.shape, (4, 10))        
        self.assertEqual(res.loc['s1', 'breaks_minor_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_major_cn_aut'], 1)
        self.assertEqual(res.loc['s1', 'breaks_aut'], 1)
        self.assertEqual(res.loc['s4', 'breaks_aut'], 4)
        self.assertEqual(res.loc['s4', 'breaks_sex'], 0)
        self.assertEqual(res.loc['s4', 'breaks_tot'], 4)
    