# %%
from matplotlib import pyplot as plt
from cns.display.plot import fig_lines, fig_CN_heatmap
from cns.process.aggregation import group_samples, add_total_cn
from cns.utils.selection import cns_head
from cns.utils.files import load_cns, load_samples
from cns.utils.assemblies import get_assembly

hg19 = get_assembly("hg19")

# Load CNS
# %%  
cns_df = load_cns("./out/PCAWG_bin_3MB.tsv")
print(cns_df.head())


# Plot heatmap for the first 50 segments
# %%
head_50 = cns_head(cns_df, 50)
fig_CN_heatmap(head_50, print_info=True, dpi=100)
plt.show()

# Plot mean total_cn per segment along the genome
# %%
groups = group_samples(cns_df)
fig_lines(groups, width=12, column="total_cn", dpi = 200)
plt.show()

# Load sample information
# %%
samples = load_samples("./out/PCAWG_samples.tsv")
print(samples.info())

# Find CN neutral samples at the autosomes
# %%
cn_netural = samples.query("ane_minor_cn_frac_aut == 0.0 & ane_major_cn_frac_aut == 0.0").index
print(len(cn_netural), "samples are CN neutral")

# Find samples with low coverage (below 95% in autosomes)
# %%
low_coverage = samples.query("cover_frac_aut < 0.95").index
print(len(low_coverage), "samples have low coverage")

# Filter out CN neutral and low coverage samples 
# %%
filtered = samples.query("(index not in @cn_netural) & (index not in @low_coverage) & whitelist")
len(filtered.index.unique())

# Find the representation of cancer types
# %%
print(filtered["type"].value_counts())

# Plot mean total_cn per segment along the genome for liver and pancreas
# %% 
liver = filtered.query("type == 'Liver-HCC'").index
liver_group = group_samples(cns_df.query("sample_id in @liver"))
pancreas = filtered.query("type == 'Panc-AdenoCA'").index
pancreas_group = group_samples(cns_df.query("sample_id in @pancreas"))
fig_lines([liver_group, pancreas_group], ["Liver-HCC", "Panc-AdenoCA"], column="total_cn", dpi=200)
plt.show()
# %%
