import unittest
import numpy as np
import pandas as pd

from cns.process.pipelines import get_genome_segments, regions_remove, regions_select
from cns.process.segments import *
from cns.process.imputation import *
from cns.process.binning import *
from cns.process.breakpoints import *
from cns.process.cluster import *
from cns.utils.assemblies import hg19, hg38
from cns.utils.conversions import tuples_to_segments

class TestSegments(unittest.TestCase):        
    def test_do_segments_overlap(self):
        segs_a = [(1, 0, 5), (1, 4, 8), (2, 10, 15)]
        segs_b = [(1, 5, 10), (2, 14, 20), (3, 25, 30)]
        gaps_hg19_segs = tuples_to_segments(hg19.gaps)
        gaps_hg38_segs = tuples_to_segments(hg38.gaps)
        self.assertTrue(do_segments_overlap(segs_a))
        self.assertFalse(do_segments_overlap(segs_b))
        self.assertFalse(do_segments_overlap(gaps_hg19_segs))
        self.assertFalse(do_segments_overlap(gaps_hg38_segs))

    def test_merge_segments(self):
        segs = [(1, 5, 10), (1, 10, 15), (2, 20, 25), (2, 25, 30), (3, 35, 40)]
        exp = [(1, 5, 15), (2, 20, 30), (3, 35, 40)]
        self.assertEqual(merge_segments(segs), exp)

    def test_segment_union(self):
        segs_a = [(1, 0, 10), (2, 10, 15)]
        segs_b = [(1, 5, 10), (2, 15, 20), (3, 25, 30)]
        exp = [(1, 0, 10), (2, 10, 20), (3, 25, 30)]
        self.assertEqual(segment_union(segs_a, segs_b), exp)

    def test_find_overlaps(self):
        segs = [(1, 1, 3), (1, 7, 9), (2, 10, 15), (3, 20, 25), (1, 3, 4), (1, 8, 10), (2, 12, 20), (3, 22, 30)]
        exp = [(1, 8, 9), (2, 12, 15), (3, 22, 25)]
        self.assertEqual(find_overlaps(segs), exp)

    def test_segment_difference(self):
        segs_a = [(1, 0, 10), (2, 15, 25), (3, 20, 30)]
        segs_b = [(1, 3, 5), (1, 7, 8), (2, 20, 23), (3, 22, 25), (3, 26, 29)]
        exp = [
            (1, 0, 3),
            (1, 5, 7),
            (1, 8, 10),
            (2, 15, 20),
            (2, 23, 25),
            (3, 20, 22),
            (3, 25, 26),
            (3, 29, 30),
        ]
        self.assertEqual(segment_difference(segs_a, segs_b), exp)

        segs_a = [(1, 0, 10)]
        segs_b = [(1, 9, 10)]
        res = segment_difference(segs_a, segs_b)
        self.assertEqual(res, [(1, 0, 9)])

    def test_filter_min_size(self):
        segs = [(1, 0, 10), (2, 15, 20), (3, 20, 30)]
        min_size = 6
        exp = [(1, 0, 10), (3, 20, 30)]
        self.assertEqual(filter_min_size(segs, min_size), exp)

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

    def test_split_segment(self):
        # Test case 1
        segment = (1, 1, 11)
        step_size = 2
        equidisant = True
        expected_output = [(1, 1, 3), (1, 3, 5), (1, 5, 7), (1, 7, 9), (1, 9, 11)]
        actual_output = split_segment(segment, step_size, equidisant)
        self.assertEqual(actual_output, expected_output)

        # Test case 2
        segment = (1, 1, 11)
        step_size = 3
        equidisant = False
        expected_output = [(1, 1, 5), (1, 5, 8), (1, 8, 11)]
        actual_output = split_segment(segment, step_size, equidisant)
        self.assertEqual(actual_output, expected_output)


class TestImputation(unittest.TestCase):
    def setUp(self):
        self.cns_df = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's2'],
            'chrom': ['chr1', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2'],
            'start': [0, 0, 50, 125, 150, 175],
            'end': [100, 150, 100, 150, 175, 200],
            'major_cn': [1, 2, 3, np.nan, 1, 1],
            'minor_cn': [1, 2, 1, 0, 0, 0]
        }) 
        self.chr_lengths = {'chr1': 100, 'chr2': 200}
        self.total_len = sum(self.chr_lengths.values())
        self.samples_df = pd.DataFrame({
            'sex': ['xx', 'xy']
        }, index=['s1', 's2'])

    def test_get_nan_regs(self):
        result = get_nan_segs(self.cns_df)
        self.assertEqual(result.shape[0], 1)

    def test_add_tails(self):
        result = add_tails(self.cns_df, self.chr_lengths)
        self.assertEqual(result.shape[0], 8)
        self.assertEqual(result.at[3, "start"], 0)
        self.assertEqual(result.at[3, "end"], 50)

    def test_fill_gaps(self):
        result = fill_gaps(self.cns_df, print_info=True)
        self.assertEqual(result.shape[0], 7)
        self.assertEqual(result.at[3, "start"], 100)
        self.assertEqual(result.at[3, "end"], 125)

    def test_add_missing(self):
        result = add_missing(self.cns_df, self.samples_df, self.chr_lengths, print_info=False)
        self.assertEqual(result.shape[0], 7)
        self.assertEqual(result.at[2, "start"], 0)
        self.assertEqual(result.at[2, "end"], 100)

    def test_merge_neighbours(self):
        result = merge_neighbours(self.cns_df, print_info=False)
        self.assertEqual(result.shape[0], 5)
        self.assertEqual(result.at[4, "start"], 150)
        self.assertEqual(result.at[4, "end"], 200)

    def test_fill_nans_with_zeros(self):
        result = fill_nans_with_zeros(self.cns_df, print_info=False)
        self.assertEqual(result.major_cn.isnull().sum(), 0)
        self.assertEqual(result.minor_cn.isnull().sum(), 0)

    def test_create_imputed_entries(self):
        result = add_tails(self.cns_df, self.chr_lengths)
        result = fill_gaps(result, print_info=False)    
        result = add_missing(result, self.samples_df, self.chr_lengths, print_info=False)
        result = create_imputed_entries(result, print_info=False)
        result = merge_neighbours(result, print_info=False)
        result = fill_nans_with_zeros(result, print_info=False)    
        self.assertEqual(result.shape[0], 7)        
        self.assertEqual(result.at[5, "end"], 137)
        self.assertEqual(result.at[6, "start"], 137)
        self.assertEqual(result.at[6, "end"], 200)
        self.assertEqual(result.major_cn.isnull().sum(), 0)


class TestBreakpoints(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_arm_breaks(self):
        result = calc_arm_breaks()
        self.assertEqual(list(result.keys())[0], 'chr1')
        self.assertEqual(list(result.values())[0], [0, 125000000, 249250621])
        
    def test_cytoband_breaks(self):
        result = calc_cytoband_breaks()
        self.assertEqual(list(result.keys())[0], 'chr1')
        sum_of_breaks = sum([len(breaks) - 1 for breaks in result.values()])
        sum_of_bands = len(hg19.cytobands)
        self.assertEqual(sum_of_breaks, sum_of_bands)

    def test_bin_breaks(self):
        result = calc_bin_breaks(1000000, True)
        self.assertEqual(list(result.keys())[0], 'chr1')
        self.assertEqual(list(result.values())[0][0], 0)
        self.assertEqual(list(result.values())[0][-1], 249250621)

    def test_dist_breaks_equidistant(self):
        act = create_step_breaks(10, 1)
        exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 2, True)
        exp = [0, 2, 4, 6, 8, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 2.5, True)
        exp = [0, 3, 5, 8, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 3)
        exp = [0, 3, 7, 10]
        self.assertEqual(act, exp)        
        act = create_step_breaks(13, 10)
        exp = [0, 13]
        self.assertEqual(act, exp)
        act = create_step_breaks(17, 10)
        exp = [0, 9, 17]
        self.assertEqual(act, exp)
        act = create_step_breaks(33, 10, True)
        exp = [0, 11, 22, 33]
        self.assertEqual(act, exp)
        act = create_step_breaks(37, 10, True)
        exp = [0, 9, 19, 28, 37]
        self.assertEqual(act, exp)
        
    def test_dist_breaks_unequidistant(self):
        act = create_step_breaks(10, 1, False)
        exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 2, False)
        exp = [0, 2, 4, 6, 8, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 2.5, False)
        exp = [0, 3, 5, 8, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(10, 3, False)
        exp = [0, 4, 7, 10]
        self.assertEqual(act, exp)
        act = create_step_breaks(13, 10, False)
        exp = [0, 13]
        self.assertEqual(act, exp)
        act = create_step_breaks(17, 10, False)
        exp = [0, 9, 17]
        self.assertEqual(act, exp)
        act = create_step_breaks(33, 10, False)
        exp = [0, 12, 22, 33]
        self.assertEqual(act, exp)
        act = create_step_breaks(37, 10, False)
        exp = [0, 9, 19, 29, 37]
        self.assertEqual(act, exp)

    def test_diffs(self):
        bin_breaks = calc_bin_breaks(10_000_000)
        for chrom, chrom_breaks in bin_breaks.items():
            diffs = np.diff(np.diff(chrom_breaks))
            self.assertTrue(np.abs(np.sum(diffs)) <= 1)
            self.assertTrue(np.max(np.abs(diffs) <= 1))

    def test_get_breaks(self):
        cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's2'],
            'chrom': ['chr1', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2'],
            'start': [0, 0, 50, 125, 150, 175],
            'end': [100, 150, 100, 150, 175, 200],
            'major_cn': [1, 2, 3, np.nan, 1, 1],
            'minor_cn': [1, 2, 1, 0, 0, 0]
        }) 
        breaks = get_breaks(cns, keep_ends=True)
        self.assertEqual(breaks['chr1'], [0, 100])
        self.assertEqual(breaks['chr2'], [0, 50, 100, 125, 150, 175, 200])
        self.assertEqual(breaks['chr3'], [])
        breaks = get_breaks(cns, keep_ends=False)
        self.assertEqual(breaks['chr1'], [100])
        # self.assertEqual(breaks['chr2'], [50, 100, 125, 150, 175, 200])
        assembly = type('Assembly', (object,), {
            'chr_lens': {'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'chr_names': ['chr1', 'chr2', 'chr3', 'chrX', 'chrY']
        })
        breaks = get_breaks(cns, assembly=assembly, keep_ends=False)
        self.assertEqual(breaks['chr1'], [])
        self.assertEqual(breaks['chr2'], [50, 100, 125, 150, 175])



class TestBinning(unittest.TestCase):
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
            'sex_names': ['chrX', 'chrY']
        })

    def test_bin_by_breaks(self):
        segments = [('chr1', 0, 100), ('chr2', 100, 200)]
        breaks = {'chr1': [0, 100], 'chr2': [100, 200]}
        seg_bin = bin_by_segments(self.cns, segments, print_progress=False)
        break_bin = bin_by_breaks(self.cns, breaks, print_progress=False)
        self.assertEqual(seg_bin.shape[0], 4)
        pd.testing.assert_frame_equal(seg_bin, break_bin)
    
    def test_bin_by_segments(self):
        segments = [('chr1', 0, 100), ('chr2', 100, 200)]
        result = bin_by_segments(self.cns, segments)
        self.assertEqual(result.shape[0], 4)
        self.assertEqual(result.at[0, "start"], 0)
        self.assertEqual(result.at[0, "end"], 100)
        self.assertEqual(result.at[0, "major_cn"], 1.5)
        self.assertEqual(result.at[0, "minor_cn"], 1.0)
        self.assertEqual(result.at[1, "start"], 100)
        self.assertTrue(np.isnan(result.at[1, "major_cn"]))
    
    def test_bin_none(self):        
        segments = [('chr1', 0, 100), ('chr2', 100, 200)]
        result = bin_by_segments(self.cns, segments, fun_type="none")
        self.assertEqual(result.shape[0], 8)
        for i in range(result.shape[0]):
            if result.at[i, "chrom"] == "chr1":
                self.assertTrue(result.at[i, "start"] >= 0)
                self.assertTrue(result.at[i, "end"] <= 100)
            elif result.at[i, "chrom"] == "chr2":
                self.assertTrue(result.at[i, "start"] >= 100)
                self.assertTrue(result.at[i, "end"] <= 200)


class TestMCS(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'chrom': ['chr1', 'chr1', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 50, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [50, 100, 300, 400, 500, 50, 99, 100, 100, 120, 130],
            'major_cn': [1, 2, 3, 4, 5, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, 0, 4, 3, 1, 0, 0, 1, 0, 1],
        })        
        self.assembly = type('Assembly', (object,), {
            'chr_lens': {'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100}
        })

    def test_prep_clusters(self):
        chrom_breaks = [50, 99]
        clusters = breaks_to_clusters(chrom_breaks)
        self.assertEqual(clusters[0, 0], 50)
        self.assertEqual(clusters[0, 1], 1)

    def test_clusters_to_breaks(self):
        clusters = [(50, 1), (99, 1)]
        breakpoints = clusters_to_breaks(clusters)
        self.assertEqual(breakpoints, [50, 99])

    def test_merge_clusters(self):
        clusters = np.array([[50, 1], [149, 1], [200, 1], [299, 1]], dtype=np.float64)
        threshold = 100
        result = merge_clusters(clusters, threshold)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0, 0], 100)
        self.assertEqual(result[0, 1], 2)
        self.assertEqual(result[1, 0], 250)
        self.assertEqual(result[1, 1], 2)

    def test_created_merged_segs(self):
        dict_start = {'chr1': [50, 99], 'chr2': [200, 300]}
        dist = 100
        assembly = type('Assembly', (object,), {'chr_lens': {'chr1': 100, 'chr2': 400, 'chr3': 300, 'chrX': 100, 'chrY': 100}})
        result = created_merged_segs(dict_start, dist, assembly, False)
        exp = [('chr1', 0, 74), ('chr1', 74, 100), ('chr2', 0, 200), ('chr2', 200, 300), ('chr2', 300, 400)]
        self.assertEqual(result, exp)
        result = created_merged_segs(dict_start, dist, assembly, True)
        exp = [('chr1', 0, 100), ('chr2', 0, 200), ('chr2', 200, 300), ('chr2', 300, 400)]
        self.assertEqual(result, exp)

    
class TestPipelines(unittest.TestCase):
    def test_get_segments(self):
        bin_size = 10000000
        filter_size = 1000000
        select = regions_select("")
        remove = regions_remove("gaps")
        segs = get_genome_segments(select, bin_size, remove, filter_size)
        print(segs)


if __name__ == "__main__":
    unittest.main()