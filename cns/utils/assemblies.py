from cns.utils.cytobands import hg19_cytobands, hg38_cytobands
from cns.utils.gaps import hg19_gaps, hg38_gaps
from cns.utils.genomes import chr_names, aut_names, human_chr_colors, hg19_chr_lengths, hg19_genome_length, hg19_autosome_length, hg19_chr_cum_starts, hg38_chr_lengths, hg38_genome_length, hg38_autosome_length, hg38_chr_cum_starts

class Assembly:
    def __init__(self, name, chr_names, chr_count, aut_names, aut_count, chr_lens, chr_starts, chr_colors, gen_len, aut_len, cytobands, gaps):
        self.name = name
        self.chr_names = chr_names
        self.chr_count = chr_count
        self.aut_names = aut_names
        self.aut_count = aut_count
        self.sex_names = list(set(chr_names) - set(aut_names))
        self.sex_count = chr_count - aut_count
        self.chr_lens = chr_lens
        self.chr_starts = chr_starts
        self.chr_colors = chr_colors
        self.gen_len = gen_len
        self.aut_len = aut_len
        self.cytobands = cytobands
        self.gaps = gaps


hg19 = Assembly(
    name="hg19",
    chr_names=chr_names,
    chr_count=24,
    aut_names=aut_names,
    aut_count=22,
    chr_lens=hg19_chr_lengths,
    chr_starts=hg19_chr_cum_starts,
    chr_colors=human_chr_colors,
    gen_len=hg19_genome_length,
    aut_len=hg19_autosome_length,
    cytobands=hg19_cytobands,
    gaps=hg19_gaps
)


hg38 = Assembly(
    name="hg38",
    chr_names=chr_names,
    chr_count=24,
    aut_names=aut_names,
    aut_count=22,
    chr_lens=hg38_chr_lengths,
    chr_starts=hg38_chr_cum_starts,
    chr_colors=human_chr_colors,
    gen_len=hg38_genome_length,
    aut_len=hg38_autosome_length,
    cytobands=hg38_cytobands,
    gaps=hg38_gaps
)


def get_assembly(assembly_id):
    if assembly_id == "hg19":
        return hg19
    elif assembly_id == "hg38":
        return hg38
    else:
        raise ValueError(f"Assembly {assembly_id} not found")