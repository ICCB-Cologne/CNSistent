import cns

to_process = [
	("../data/COSMIC_consensus_genes.bed", "../out/COSMIC_with_fill.bed"),
	("../data/ENSEMBL_coding_genes.bed", "../out/ENSEMBL_with_fill.bed"),
]

for input_file, output_file in to_process:
	segment_source = cns.load_segments(input_file)
	missing = cns.genome_to_segments(cns.hg19)
	negative = cns.segment_difference(missing, segment_source)
	merged = cns.segment_union(segment_source, negative, False)
	processed = cns.process_segments(merged, cns.regions_select("gaps"), 50001)
	cns.save_segments(processed, output_file)