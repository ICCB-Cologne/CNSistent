import os
from cns.data_utils import *
from cns.pipelines import main_segment
from cns.utils.files import save_segments
from cns.process.segments import regions_select
import argparse

if __name__ == "__main__":
	# require 1 integer parameter
	parser = argparse.ArgumentParser()
	parser.add_argument("dist", type=int, help="distance for clustering")
	args = parser.parse_args()
	dist = args.dist
	if dist <= 0:
		raise ValueError("Distance must be greater than 0")

	samples_df, cns_df = main_load()

	select = regions_select("whole")
	remove = regions_select("gaps")
	clustered = main_segment(cns_df, remove, merge_dist=dist, filter_size=dist//10)
	file = os.path.join(out_path, f'segs_merge_{dist}.bed')
	print("Saving to file:", file)
	save_segments(clustered, file)


