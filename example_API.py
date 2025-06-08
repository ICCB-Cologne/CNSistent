# %%
import cns

# Load CNS data a display first 5 rows
# %%
raw_df = cns.load_cns("./data/20220803_TxPri_mphase_by_sample_df.reduced.csv", cn_columns=["nMajor", "nMinor"], sep=",", print_info=True)
cns.fig_heatmap(cns.cns_head(raw_df, 5), max_cn=6)

# Add missing segments, display first 5 rows
# %%
imp_df = cns.main_fill_imp(raw_df, print_info=True)
cns.fig_heatmap(cns.cns_head(imp_df, 5), max_cn=6)

# Create 3 mb segments, convert to a 3D feature array
# %%
seg_df = cns.main_seg_agg(imp_df, split_size=3_000_000, print_info=True)
features, rows, columns = cns.bins_to_features(seg_df)
print("Samples: {0}, Alleles: {1}, Bins: {2}.".format(*features.shape))

# Group segments by cancer type, sum the CNs and create mean linear profile
# %%
sample_df = cns.load_samples("./data/20221109_TRACERx421_all_patient_df.tsv")
type_groups = {c: cns.select_cns_by_type(seg_df, sample_df, c, "histology_multi_full") for c in ["LUAD", "LUSC"]}	
groups_df = cns.stack_groups([cns.group_samples(v, group_name=k) for k, v in type_groups.items()])
cns.fig_lines(cns.add_total_cn(groups_df), cn_columns="total_cn")
