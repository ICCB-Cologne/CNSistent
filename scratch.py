# %%
import cns.data_utils as cdu
import cns
import pandas as pd
# %%
samples_df, cns_df = cdu.main_load("imp", "PCAWG", False)
# %%
cns_df

# %%
for sample_id in samples_df.index:
    single = cns_df.query(f"sample_id == '{sample_id}'")
    second = single.copy()
    second["sample_id"] = "second_sample"
    double = pd.concat([single, second])
    cns.save_cns(double, f"./for_med/{sample_id}_cns.csv")
# %%
