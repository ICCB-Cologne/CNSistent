import pandas as pd
import numpy as np
import cns.data_utils as cdu
import cns
import time

seg_size = "1MB"
print(f"Loading data for {seg_size} segments...")
samples_df, cns_df = cdu.main_load(seg_size)
subs_df = pd.concat([cns.select_cns_by_type(cns_df, samples_df, my_type) for my_type in samples_df["type"].value_counts().index[:6]])
subs_df = cns.only_aut(cns.add_total_cn(subs_df))
print("Total samples:", len(subs_df["sample_id"].unique()))
# subs_df = cns.cns_head(subs_df, 100)

cn_column = "total_cn"
print(f"Calculating distances for {cn_column}...")
start = time.time()
arr = cns.calc_wass_distances(subs_df, cn_column, print_progress=True)
print(f"Elapsed time: {time.time() - start:.2f} seconds")
out_path = cdu.pjoin(cdu.out_path, "top_6_wass_dist.tsv")
print(f"Saving distances to {out_path}...")
arr.to_csv(out_path, sep="\t", index=True, header=True)