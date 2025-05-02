import cns
import cns.data_utils as cdu
import pandas as pd
import matplotlib.pyplot as plt

# set color_map to tab10
color_map = plt.get_cmap('tab10').colors[:10]
plt.rcParams.update({'font.size': 12})

cns_dfs = {}
for grouping in ["10MB", "5MB", "3MB", "2MB", "1MB", "500KB", "250KB"]:
	print("loading", grouping)
	samples_df, cns_df = cdu.main_load(grouping)
	cns_dfs[grouping] = cns_df

cosmic = cdu.load_COSMIC()
cosmic_df = cns.segments_to_cns_df(cosmic)[["chrom", "start", "end", "name"]].rename(columns={"name": "gene"})
ensembl = cdu.load_ENSEMBL()
cancer_type = "LUSC"
val_count = 5

top_types = samples_df["type"].value_counts().index.tolist()[:6]
for i, cancer_type in enumerate(["all"] + top_types):
	print(i, cancer_type)
	fig, axs = plt.subplots(len(cns_dfs), 1, figsize=(14, 14))
	for i, (grouping, cns_df) in enumerate(cns_dfs.items()):
		sel_df = cns.select_cns_by_type(cns_df, samples_df, cancer_type) if cancer_type != "all" else cns_df
		sel_df = cns.group_samples(cns.only_aut(cns.add_total_cn(sel_df)))
		sel_df["sample_id"] = f"mean {cancer_type} CN"
		sel_df["score"] = cns.calc_angles(sel_df, "total_cn")

		cns.plot_lines(axs[i], sel_df, cn_column="total_cn", color=color_map[i])
		cns.plot_x_lines(axs[i])
		cns.plot_x_ticks(axs[i])

		sel_df = cns.add_cum_mid(sel_df)
		sel_df = sel_df.sort_values(by="score")
		axs[i].scatter(sel_df["cum_mid"].head(val_count), sel_df["total_cn"].head(val_count), color="k", alpha=0.5, s=15, label=f"Top {val_count} peaks", marker="+")
		axs[i].scatter(sel_df["cum_mid"].tail(val_count), sel_df["total_cn"].tail(val_count), color="k", alpha=0.5, s=15, label=f"Top {val_count} valleys", marker="X")
		axs[i].set_ylim(0, 8)
		axs[i].set_ylabel("Total CN")
		axs[i].legend(title=f"{grouping} bins", loc="upper right", ncol=3)

	axs[-1].set_xlabel("Poisition on linear genome")

	cdu.save_cns_fig(f"peaks_valleys_{cancer_type}")

	score_means = []
	for i, (grouping, cns_df) in enumerate(cns_dfs.items()):
		print(grouping)
		sel_df = cns.select_cns_by_type(cns_df, samples_df, cancer_type) if cancer_type != "all" else cns_df
		sel_df = cns.group_samples(cns.only_aut(cns.add_total_cn(sel_df)))
		sel_df["sample_id"] = f"mean {cancer_type} CN"
		sel_df["score"] = cns.calc_angles(sel_df, "total_cn")	
		score_means.append(cns.mean_value_per_seg(sel_df, ensembl, "score"))

	mean_dfs = {}
	mean_df = score_means[0].copy()
	for vals in score_means[1:]:
		mean_df["score"] += vals["score"]
	mean_df["score"] /= len(score_means)
	mean_df["total_cn"] = cns.mean_value_per_seg(sel_df, ensembl, "total_cn")["total_cn"]
	mean_df = pd.merge(mean_df, cosmic_df, how="left")

	cns.save_cns(mean_df, cdu.pjoin(cdu.out_path, f"gene_scores_{cancer_type}.tsv"))