import os
from cns.data_utils import *
from cns.pipelines import main_segment
from cns.utils.files import save_segments
from cns.process.segments import regions_select
import argparse

if __name__ == "__main__":
	# require 1 integer parameter
	parser = argparse.ArgumentParser()
	parser.add_argument("dist", type=str, help="distance for clustering")
	args = parser.parse_args()
	dist_int = int(args.dist.replace("KB", "000").replace("MB", "000000"))
	if dist_int <= 0:
		raise ValueError("Distance must be greater than 0")

	samples_df, cns_df = main_load("imp")

	remove = regions_select("gaps")
	clustered = main_segment(cns_df, remove, cluster_dist=dist_int, filter_size=dist_int//10, print_info=True)
	file = os.path.join(out_path, f'segs_merge_{args.dist}.bed')
	print("Saving to file:", file)
	save_segments(clustered, file)


