# SimCHA aneuploidy plots

# %%
from matplotlib import pyplot as plt
import numpy as np
from cns.process.binning import add_derived
from cns.utils.files import load_cna, load_samples
from cns.utils.assemblies import hg19

# %%  
all_cns_df = load_cna("./out/PCAWG_cna_imp.tsv")
derived = add_derived(all_cns_df, hg19)
samples = load_samples("./out/PCAWG_samples.tsv")
cn_alt_df = samples.query("ane_minor_cn_frac_aut >= 0.01 & ane_major_cn_frac_aut >= 0.01").index
covered_df = samples.query("cover_frac_aut >= 0.95").index
filtered = samples.query("(index in @cn_alt_df) & (index in @covered_df) & whitelist")
print("samples", filtered.index.unique().shape[0])
cns_df = derived.query("sample_id in @filtered.index").query("chrom != 'chrX' & chrom != 'chrY'")

# %%
dists ={}
for column in ["total_cn", "minor_cn", "major_cn"]:
	vals = cns_df[(cns_df[column] > 0)].groupby('sample_id')['length'].sum()
	dists[column] = 1 - (vals / hg19.aut_len)
dists["expected"] = dists["minor_cn"] * dists["major_cn"]

# %%
bins = np.linspace(0, 0.01, 100)
plt.hist(dists["total_cn"], bins=bins, alpha=0.5, label='Total CN', histtype='step')
plt.hist(dists["expected"], bins=bins, alpha=0.5, label='Expected', histtype='step')
plt.legend(loc='upper right')

# %%
bins = np.linspace(0, 1, 100)
plt.hist(dists["major_cn"], bins=bins, alpha=0.5, label='Major CN', histtype='step')
plt.hist(dists["minor_cn"], bins=bins, alpha=0.5, label='Minor CN', histtype='step')
plt.legend(loc='upper right')

# %%
dists["major_cn"].describe()