# %%
from cns.data_utils import load_samples_out
import matplotlib.pyplot as plt

# %%
pcawg_samples = load_samples_out("PCAWG_samples.tsv")

pcawg_samples.columns
# %%
total_ploidy = pcawg_samples["ploidy_major_cn"] + pcawg_samples["ploidy_minor_cn"]

# %%
plt.hist(total_ploidy, bins=100, density=True, histtype='step')
plt.xlim(0, 14)
# %%
