import os
from cns.data_utils import *
from cns.process.pipelines import main_segment
from cns.utils.conversions import segs_to_df
from cns.utils.files import save_segments
from cns.process.segments import regions_select, regions_remove
import argparse

if __name__ == "__main__":
	# require 1 integer parameter
	parser = argparse.ArgumentParser()
	parser.add_argument("dist", type=int, help="distance for clustering")
	args = parser.parse_args()
	dist = args.dist
	if dist <= 0:
		raise ValueError("Distance must be greater than 0")

	samples = load_merged_samples(print_info=False)
	cns = load_merged_cns(samples)

	select = regions_select("")
	remove = regions_remove("gaps")
	clustered = main_segment(cns, select, remove, merge_dist=dist, filter_size=dist//10)
	res_df = segs_to_df(clustered)
	file = os.path.join(out_path, f'segs_merge_{dist}.bed')
	print("Saving to file:", file)
	save_segments(res_df, file, False, False)


