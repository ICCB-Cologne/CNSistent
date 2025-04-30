import pandas as pd
import numpy as np
import cns.data_utils as cdu
import cns
import time

print("Loading data...")

samples_df, cns_df = cdu.main_load("1MB")
cns.add_total_cn(cns_df, inplace=True)
subs_df = cns.cns_head(cns_df, 100).set_index("sample_id")

# start timer
start = time.time()

ids = subs_df.index.unique()
arr  = np.zeros((len(ids), len(ids)))
total = len(ids) **2
for i1, id1 in enumerate(ids):
    for i2, id2 in enumerate(ids[:i1]):
        print(f"{(i1*len(ids) + i2)/total:.2%}", end="\r")
        cns1 = subs_df.loc[id1]
        cns2 = subs_df.loc[id2]
        dist = cns.calc_wass_distance(cns1, cns2, "total_cn")
        arr[i1, i2] = dist.mean()

print(f"Elapsed time: {time.time() - start:.2f} seconds")

arr = pd.DataFrame(arr, index=ids, columns=ids)
arr.to_csv(cdu.pjoin(cdu.out_path, "LUAD_LUSC_wass_dist.tsv"), sep="\t", index=True, header=True)



