# %%
from cns.data_utils import main_load
import matplotlib.pyplot as plt
from cns.utils import *

# %%
samples_df, cns_df = main_load("10MB", "TRACERx")

cn_columns = ["major_cn", "minor_cn"]
drop_sex = True

# %%
def bins_to_features(cns_df, cn_columns=None, drop_sex=True):
	cn_columns = get_cn_cols(cns_df, cn_columns)
	sel_df = only_aut(cns_df, inplace=False) if drop_sex else cns_df
	groups = sel_df.groupby("sample_id")
	columns = next(iter(groups))[1][['chrom', 'start', 'end']]
	rows = list(groups.groups.keys())

	# Calculate the cumulative count for each group without assigning it to the DataFrame
	cumcount = groups.cumcount()

	if len(cumcount) != len(columns) * len(rows):
		raise ValueError("The number of cumulative counts does not match the number of rows and columns. Make sure that each sample has the same number of bins.")

	arrays = []
	for cn_col in cn_columns:
		# Use the cumulative count directly in the pivot operation
		array = sel_df.pivot_table(index="sample_id", columns=cumcount, values=cn_col)
		arrays.append(array)

	stacked = np.stack(arrays, axis=0)

	return stacked, rows, columns

features, rows, columns = bins_to_features(cns_df, cn_columns, drop_sex)