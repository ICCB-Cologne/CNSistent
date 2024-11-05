from .cytobands import hg19_cytobands, hg38_cytobands
from .gaps import hg19_gaps, hg38_gaps
from .genomes import aut_names, x_name, y_name, human_chr_colors, hg19_chr_lengths, hg19_genome_length, hg19_autosome_length, hg19_chr_starts, hg38_chr_lengths, hg38_genome_length, hg38_autosome_length, hg38_chr_starts

class Assembly:
    """
    A class to represent a genomic assembly.

    Attributes
    ----------
    name : str
        The name of the assembly.
    aut_names : list
        The names of the autosomes.
    chr_x : str
        The name of the X chromosome.
    chr_y : str
        The name of the Y chromosome.
    chr_lens : dict
        The lengths of the chromosomes.
    chr_starts : dict
        The start positions of the chromosomes.
    chr_colors : dict
        The colors of the chromosomes.
    gen_len : int
        The total length of the genome.
    aut_len : int
        The total length of the autosomes.
    cytobands : list
        The cytobands of the chromosomes.
    gaps : list
        The gaps in the chromosomes.
    """

    def __init__(self, name, aut_names, x_name, y_name, chr_lens, chr_starts, chr_colors, gen_len, aut_len, cytobands, gaps):
        """
        Constructs all the necessary attributes for the Assembly object.

        Parameters
        ----------
        name : str
            The name of the assembly.
        aut_names : list
            The names of the autosomes.
        x_name : str
            The name of the X chromosome.
        y_name : str
            The name of the Y chromosome.
        chr_lens : dict
            The lengths of the chromosomes.
        chr_starts : dict
            The start positions of the chromosomes.
        chr_colors : dict
            The colors of the chromosomes.
        gen_len : int
            The total length of the genome.
        aut_len : int
            The total length of the autosomes.
        cytobands : list
            The cytobands of the chromosomes.
        gaps : list
            The gaps in the chromosomes.
        """
        self.name = name
        self.aut_names = aut_names
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
"""
An instance of the Assembly class representing the hg19 genomic assembly.
"""

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
"""
An instance of the Assembly class representing the hg38 genomic assembly.
"""


def get_assembly(assembly_id):
    """
    Retrieve an Assembly instance by its ID.

    Parameters
    ----------
    assembly_id : str
        The ID of the assembly to retrieve. Valid values are "hg19" and "hg38".

    Returns
    -------
    Assembly
        The Assembly instance corresponding to the given ID.

    Raises
    ------
    ValueError
        If the assembly_id is not "hg19" or "hg38".
    """
    if assembly_id == "hg19":
        return hg19
    elif assembly_id == "hg38":
        return hg38
    else:
        raise ValueError(f"Assembly {assembly_id} not found")