""" This module contains code for dealing with amino acid sequence and coding sequence (CDS) data. """
import math
import dataclasses
import python_codon_tables

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
    
    # static method to get list of all available tables
    @staticmethod
    def get_available_tables():
        return python_codon_tables.get_available_tables()
    
    def __init_from_file(self, file_path):
        if file_path.endswith("homosapiens.txt"): return 9606
        if file_path.endswith("mouse.txt"): return 10090
    
    def __init__(self, taxid):
        # if taxid is a string, check if it is a file path
        if isinstance(taxid, str):
            if taxid.endswith(".txt"):
                taxid = self.__init_from_file(taxid)
        self.taxid = taxid
        self.table = python_codon_tables.get_codons_table(taxid, replace_U_by_T=False)
        self.codon_to_aa = {}
        self.aa_to_codons = {}
        for aa in self.table.keys():
            self.aa_to_codons[aa] = set()
            for codon in self.table[aa].keys():
                self.codon_to_aa[codon] = aa
                self.aa_to_codons[aa].add(codon)
    
    def get_codon_freq(self, codon):
        total = 0
        for aa in self.table.keys():
            # if codon is in the table under the current amino acid
            if codon in self.table[aa]:
                total += self.table[aa][codon]
        return total
    
    def get_aa_max_freq(self, aa):
        return max(self.table[aa][codon] for codon in self.table[aa].keys())
    
    def get_codons(self, aa) -> set[str]:
        return self.table[aa].keys()
    
    def get_aa(self, codon):
        return self.codon_to_aa[codon]

    # Maximum number of codons for a single amino acid
    def max_codons(self) -> int:
        return max(len(self.get_codons(aa)) for aa in self.amino_acids)
    
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
    
    @property
    def amino_acids(self):
        return self.table.keys()
    
    @property
    def name(self):
        return self.table.name


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


def random_cds(aa_seq, freq_table, seed=None) -> list[str]:
    import random
    generator = random.Random(seed)
    cds = []
    for aa in aa_seq:
        cds.append(generator.choice(list(freq_table.aa_to_codons[aa])))
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
