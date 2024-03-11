""" This module contains code for dealing with amino acid sequence and coding sequence (CDS) data. """
import math
import dataclasses

AA_SINGLE_LETTER = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "End": "*",  # Stop codon
}


def is_valid_aa_letter(aa):
    return aa in AA_SINGLE_LETTER.values() and aa != '*'


class CodonFrequencyTable:
    def __init__(self, table_path):
        file = open(table_path, 'r')
        self.codons = set()
        self.codon_to_aa = {}
        self.aa_to_codons = {}
        self.codon_freq = {}
        self.aa_max_freq = {}
        # The format uses T instead of U
        self.nt_map = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
        for line in file:
            tokens = line.strip(" \n").split()
            if len(tokens) < 3:
                continue
            aa = tokens[0]
            aa = AA_SINGLE_LETTER[aa]
            codon = ''.join([self.nt_map[nt] for nt in tokens[1]])
            freq = round(float(tokens[2]))
            self.codons.add(codon)
            self.codon_to_aa[codon] = aa
            if aa not in self.aa_to_codons:
                self.aa_to_codons[aa] = set()
            self.aa_to_codons[aa].add(codon)
            self.codon_freq[codon] = freq
            if aa not in self.aa_max_freq:
                self.aa_max_freq[aa] = 0
            self.aa_max_freq[aa] = max(self.aa_max_freq[aa], freq)
        file.close()

    def get_codon_freq(self, codon):
        return self.codon_freq[codon]

    def get_aa_max_freq(self, aa):
        return self.aa_max_freq[aa]

    def get_codons(self, aa) -> set[str]:
        return self.aa_to_codons[aa]

    # Maximum number of codons for a single amino acid
    def max_codons(self) -> int:
        return max(len(self.get_codons(aa)) for aa in self.aa_to_codons)

    def get_aa(self, codon):
        return self.codon_to_aa[codon]

    def codon_adaption_weight(self, codon):
        return self.get_codon_freq(codon) / self.get_aa_max_freq(self.get_aa(codon))

    def codon_adaptation_index(self, cds) -> float:
        cai = 1
        for codon in cds:
            cai *= self.codon_adaption_weight(codon)
        return cai**(1/len(cds))

    def log_codon_adaptation_index(self, cds):
        cai = 0
        for codon in cds:
            cai += math.log(self.codon_adaption_weight(codon))
        return cai / len(cds)


@dataclasses.dataclass
class UniProtAASeq:
    seq: str
    uniprot_name: str
    protein_name: str


def read_cds(tsv_path):
    cds = []
    first = True
    for line in open(tsv_path, 'r'):
        if first:
            first = False
            continue
        tokens = line.strip(" \n").split("\t")
        if len(tokens) != 7:
            assert False, "Invalid CDS file format"
        if any(not is_valid_aa_letter(aa) for aa in tokens[6]):
            continue
        cds.append(UniProtAASeq(tokens[6], tokens[1], tokens[2]))
    return cds


def random_cds(aa_seq, freq_table) -> list[str]:
    import random
    cds = []
    for aa in aa_seq:
        cds.append(random.choice(list(freq_table.aa_to_codons[aa])))
    return cds


def encode_cds_one_hot(cds: list[str], freq_table: CodonFrequencyTable) -> list[list[float]]:
    width = freq_table.max_codons()
    res = []
    for codon in cds:
        enc = [0.0 for _ in range(width)]
        codons = sorted(freq_table.get_codons(freq_table.get_aa(codon)))
        enc[codons.index(codon)] = 1.0
        res.append(enc)
    return res


def main():
    """Test the codon frequency table class."""
    path = "../data/codon_tables/homosapiens.txt"
    ct = CodonFrequencyTable(path)
    print(ct.get_codon_freq("AUG"))

    cds = ["AUG", "GGC", "UUG", "UCC", "CGG", "AGC", "GAG"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "CCC", "CCC", "GGG", "GGC", "UUA", "AAA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "GUU", "AAA", "GUA", "GUA", "GGG", "GUA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)


if __name__ == "__main__":
    main()
