# %%
%load_ext autoreload
%autoreload 2

# %%
import pandas as pd
from cns import main_fill_imp, main_seg_agg
from cns.analyze import fig_lines, fig_heatmap
from cns.process import group_samples
from cns.utils import load_cns, load_samples, cns_head, select_cns_by_type, bins_to_features

# %%
raw_df = load_cns("./data/20220803_TxPri_mphase_by_sample_df.reduced.csv", cn_columns=["nMajor", "nMinor"], sep=",", canonize=True, print_info=True)
fig_heatmap(cns_head(raw_df, 5), max_cn=6)

# %%
sample_df = load_samples("./data/20221109_TRACERx421_all_patient_df.tsv", cananonize=True)
imp_df = main_fill_imp(raw_df, print_info=True)
fig_heatmap(cns_head(imp_df, 5), max_cn=6)

# %%
seg_df = main_seg_agg(imp_df, split_size=3_000_000, print_info=True)
features, rows, columns = bins_to_features(seg_df)
print(features.shape)

# %%
type_groups = {c: select_cns_by_type(seg_df, sample_df, c, "histology_multi_full") for c in ["LUAD", "LUSC"]}	
groups_df = pd.concat([group_samples(v, group_name=k) for k, v in type_groups.items()])
fig_lines(groups_df)
