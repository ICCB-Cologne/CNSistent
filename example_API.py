# %%
%load_ext autoreload
%autoreload 2

# %%
import pandas as pd
import cns

# %%
raw_df = cns.load_cns("./data/20220803_TxPri_mphase_by_sample_df.reduced.csv", cn_columns=["nMajor", "nMinor"], sep=",", print_info=True)
cns.fig_heatmap(cns.cns_head(raw_df, 5), max_cn=6)

# %%
sample_df = cns.load_samples("./data/20221109_TRACERx421_all_patient_df.tsv")
imp_df = cns.main_fill_imp(raw_df, print_info=True)
cns.fig_heatmap(cns.cns_head(imp_df, 5), max_cn=6)

# %%
seg_df = cns.main_seg_agg(imp_df, split_size=3_000_000, print_info=True)
features, rows, columns = cns.bins_to_features(seg_df)
print(features.shape)

# %%
type_groups = {c: cns.select_cns_by_type(seg_df, sample_df, c, "histology_multi_full") for c in ["LUAD", "LUSC"]}	
groups_df = pd.concat([cns.group_samples(v, group_name=k) for k, v in type_groups.items()])
cns.fig_lines(groups_df)
