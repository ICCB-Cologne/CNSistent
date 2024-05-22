# %%
from matplotlib import pyplot as plt
from cns.analyze.genome import fig_genome, fig_genome_groups
from cns.display.heatmap import fig_CN_heatmap
from cns.process.binning import group_bins
from cns.process.binning import add_cns_loc
from cns.utils.selection import cns_head
from cns.utils.files import load_cns, load_samples
from cns.utils.assemblies import get_assembly
hg19 = get_assembly("hg19")

# Load bins
# %%  
bins, cn_cols = load_cns("./out/PCAWG_bin_3MB.tsv")
print(bins.head())

# Add derived columns (length, cumulative mid,...)
# %%
derived = add_cns_loc(bins)
print(derived.head())

# Plot heatmap for the first 50 bins
# %%
head_50 = cns_head(derived, 50)
fig_CN_heatmap(head_50, print_info=True, dpi=100)
plt.show()

# Plot mean total_cn per bin along the genome
# %%
groups = group_bins(derived)
fig_genome(groups, max_cn=10, width=12, dpi = 200, size=.5, colored=True)
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

# Plot mean total_cn per bin along the genome for liver and pancreas
# %% 
liver = filtered.query("type == 'Liver-HCC'").index
liver_group = group_bins(derived.query("sample_id in @liver"))
pancreas = filtered.query("type == 'Panc-AdenoCA'").index
pancreas_group = group_bins(derived.query("sample_id in @pancreas"))
fig_genome_groups([liver_group, pancreas_group], ["Liver-HCC", "Panc-AdenoCA"], hg19, column="major_cn", max_cn=3, width=12, dpi=200, size=.5, colored=False)
plt.show()