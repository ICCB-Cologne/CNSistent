import pandas as pd
import numpy as np
import cns.data_utils as cdu
import cns
import time

print("Loading data...")

samples_df, cns_df = cdu.main_load("10MB")
subs_df = cns.only_aut(cns.add_total_cn(cns_df))
subs_df = cns.cns_head(subs_df, 50)

print("Calculating distances...")
start = time.time()
arr = cns.calc_distances(subs_df, True)
print(f"Elapsed time: {time.time() - start:.2f} seconds")
arr.to_csv(cdu.pjoin(cdu.out_path, "LUAD_LUSC_wass_dist_2.tsv"), sep="\t", index=True, header=True)