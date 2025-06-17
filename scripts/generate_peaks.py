import cns

# This script analyzes copy number segment data, computes angle-based scores, visualizes results, and saves gene-level statistics.
import cns.data_utils as cdu
import pandas as pd
import matplotlib.pyplot as plt


color_map = list(plt.cm.get_cmap('tab10').colors[:10]) + [(0,0,0)]
plt.rcParams.update({'font.size': 12})

cns_dfs = {}
for grouping in ["whole", "arms", "20MB", "10MB", "5MB", "3MB", "2MB", "1MB", "500KB", "250KB", "100KB"]:
	print("loading", grouping)
	samples_df, cns_df = cdu.main_load(grouping)
	cns_dfs[grouping] = cns_df

cosmic = cdu.load_COSMIC()
cosmic_df = cns.segments_to_cns_df(cosmic)[["chrom", "start", "end", "name"]].rename(columns={"name": "gene"})
ensembl = cdu.load_ENSEMBL()
val_count = 3

def get_gouping_type(bins):
	if bins == "whole":
		return "sample"
	elif bins == "arms":
		return "chrom"
	else:
		return "cons"

top_types = samples_df["type"].value_counts().index.tolist()[:6] # top 6 types 
for i, cancer_type in enumerate(["all"] + top_types):
	print(i+1, ":", cancer_type)
	fig, axs = plt.subplots(len(cns_dfs), 1, figsize=(14, 20))

	for j, (grouping, cns_df) in enumerate(cns_dfs.items()):
		print("Plot", grouping)
		sel_df = cns.select_cns_by_type(cns_df, samples_df, cancer_type) if cancer_type != "all" else cns_df
		group_df = cns.group_samples(cns.only_aut(cns.add_total_cn(sel_df)))
		group_df["sample_id"] = f"mean {cancer_type} CN"
		group_df["score"] = cns.calc_angles(group_df, "total_cn", group_by=get_gouping_type(grouping))

		if grouping == "whole" or grouping == "arms":
			cns.plot_bars(axs[j], group_df, cn_column="total_cn", color=color_map[j], label=f"{grouping} bins")
		else:
			cns.plot_lines(axs[j], group_df, cn_column="total_cn", color=color_map[j], label=f"{grouping} bins")
		cns.plot_x_lines(axs[j])
		cns.plot_x_ticks(axs[j])

		group_df = group_df.sort_values(by="score")
		group_df = cns.add_cum_mid(group_df)
		axs[i].scatter(group_df["cum_mid"].head(val_count), group_df["total_cn"].head(val_count), color="k", alpha=0.75, s=25, label=f"Top {val_count} peaks", marker="+")
		axs[i].scatter(group_df["cum_mid"].tail(val_count), group_df["total_cn"].tail(val_count), color="k", alpha=0.75, s=25, label=f"Top {val_count} valleys", marker="X")
		axs[j].set_ylim(0, 8)
		axs[j].set_ylabel("Total CN")
		axs[j].legend(loc="upper right")

	axs[-1].set_xlabel("Poisition on linear genome")

	cdu.save_cns_fig(f"peaks_valleys_{cancer_type}")

	score_means = []
	for j, (grouping, cns_df) in enumerate(cns_dfs.items()):
		print("Score", grouping)
		sel_df = cns.select_cns_by_type(cns_df, samples_df, cancer_type) if cancer_type != "all" else cns_df
		sel_df = cns.group_samples(cns.only_aut(cns.add_total_cn(sel_df)))
		sel_df["sample_id"] = f"mean {cancer_type} CN"
		sel_df["score"] = cns.calc_angles(sel_df, "total_cn", get_gouping_type(grouping))	
		score_means.append(cns.mean_value_per_seg(sel_df, ensembl, "score"))

	mean_dfs = {}
	mean_df = score_means[0].copy()
	for vals in score_means[1:]:
		mean_df["score"] += vals["score"]
	mean_df["score"] /= len(score_means)
	mean_df["total_cn"] = cns.mean_value_per_seg(sel_df, ensembl, "total_cn")["total_cn"]
	mean_df = pd.merge(mean_df, cosmic_df, how="left")

	cns.save_cns(mean_df, cdu.pjoin(cdu.out_path, f"peak_scores_{cancer_type}.tsv"))