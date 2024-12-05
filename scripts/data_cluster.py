import os
import argparse
import cns.data_utils as cdu
import cns

if __name__ == "__main__":
	# require 1 integer parameter
	parser = argparse.ArgumentParser()
	parser.add_argument("dist", type=str, help="distance for clustering")
	args = parser.parse_args()
	dist_int = int(args.dist.replace("KB", "000").replace("MB", "000000"))
	if dist_int <= 0:
		raise ValueError("Distance must be greater than 0")

	samples_df, cns_df = cdu.main_load("imp")

	remove = cns.regions_select("gaps")
	clustered = cns.main_segment(cns_df, remove, cluster_dist=dist_int, filter_size=dist_int//10, print_info=True)
	file = os.path.join(cdu.out_path, f'segs_merge_{args.dist}.bed')
	print("Saving to file:", file)
	cns.save_segments(clustered, file)


