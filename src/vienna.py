"""Interface to ViennaRNA package"""
from typing import Tuple
import RNA


class ViennaContext:
    def __init__(self, rna, temp=None, dangles=2, noLPs=False) -> None:
        md = RNA.md()
        md.uniq_ML = 1
        md.dangles = dangles
        md.noLP = noLPs
        if temp is not None:
            md.temperature = temp
        self.fc = RNA.fold_compound(rna, md)
        self.pf_computed = False

    def __ensure_pf(self):
        if self.pf_computed:
            return
        _, mfe_energy = self.fc.mfe()
        self.fc.exp_params_rescale(mfe_energy)
        _, efe = self.fc.pf()
        self.efe = efe
        self.pf_computed = True

    def free_energy(self, ss):
        return self.fc.eval_structure(ss)

    def prob(self, ss):
        self.__ensure_pf()
        return self.fc.pr_structure(ss)

    def make_bppt(self) -> list[list[float]]:
        self.__ensure_pf()
        bpp = self.fc.bpp()
        sz = self.fc.length
        res = [[0.0 for _ in range(sz)] for _ in range(sz)]
        for i in range(sz):
            for j in range(sz):
                if j < i:
                    res[i][j] = bpp[j+1][i+1]
                elif i < j:
                    res[i][j] = bpp[i+1][j+1]
        for i in range(sz):
            res[i][i] = 1-sum(res[i])
        return res

    def ensemble_free_energy(self):
        self.__ensure_pf()
        return self.efe

    def subopt(self, energy_delta):
        sub = self.fc.subopt(int(energy_delta*100), sorted=0)
        return [s.structure for s in sub]

    def ensemble_defect(self, ss):
        self.__ensure_pf()
        return self.fc.ensemble_defect(ss)

    def mfe(self):
        return self.fc.mfe()[0]

    def psample(self, samples=1, redundant=True):
        self.__ensure_pf()
        return self.fc.pbacktrack(samples, RNA.PBACKTRACK_DEFAULT if redundant else RNA.PBACKTRACK_NON_REDUNDANT)


def ensemble_unpaired(bppt: list[list[float]]):
    return sum([bppt[i][i] for i in range(len(bppt))])


def average_unpaired(rna_seq: str) -> float:
    ctx = ViennaContext(rna_seq)
    return ensemble_unpaired(ctx.make_bppt())/len(rna_seq)


def ensemble_free_energy(rna_seq: str) -> float:
    ctx = ViennaContext(rna_seq)
    return ctx.ensemble_free_energy()

def aup_and_efe(rna_seq: str) -> Tuple[float, float]:
    ctx = ViennaContext(rna_seq)
    return ensemble_unpaired(ctx.make_bppt())/len(rna_seq), ctx.ensemble_free_energy()


def cds_to_rna(cds):
    return ''.join(cds)
