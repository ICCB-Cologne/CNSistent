from cns.utils.cytobands import hg19_cytobands, hg38_cytobands
from cns.utils.gaps import hg19_gaps, hg38_gaps
from cns.utils.genomes import aut_names, x_name, y_name, human_chr_colors, hg19_chr_lengths, hg19_genome_length, hg19_autosome_length, hg19_chr_starts, hg38_chr_lengths, hg38_genome_length, hg38_autosome_length, hg38_chr_starts

class Assembly:
    def __init__(self, name, aut_names, x_name, y_name, chr_lens, chr_starts, chr_colors, gen_len, aut_len, cytobands, gaps):
        self.name = name
        self.chr_names = aut_names + [x_name, y_name]
        self.chr_count = len(self.chr_names)
        self.aut_names = aut_names
        self.aut_count = len(self.aut_names)
        self.sex_names = [x_name, y_name]
        self.sex_count = len(self.sex_names)
        self.chr_x = x_name
        self.chr_y = y_name
        self.chr_lens = chr_lens
        self.chr_starts = chr_starts
        self.chr_colors = chr_colors
        self.gen_len = gen_len
        self.aut_len = aut_len
        self.cytobands = cytobands
        self.gaps = gaps


hg19 = Assembly(
    name="hg19",
    aut_names=aut_names,
    x_name=x_name,
    y_name=y_name,
    chr_lens=hg19_chr_lengths,
    chr_starts=hg19_chr_starts,
    chr_colors=human_chr_colors,
    gen_len=hg19_genome_length,
    aut_len=hg19_autosome_length,
    cytobands=hg19_cytobands,
    gaps=hg19_gaps
)


hg38 = Assembly(
    name="hg38",
    aut_names=aut_names,
    x_name=x_name,
    y_name=y_name,
    chr_lens=hg38_chr_lengths,
    chr_starts=hg38_chr_starts,
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