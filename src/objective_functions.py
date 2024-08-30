"""Implementations of objective functions for use in ARWs."""
from typing import Callable, List
import math
import dataclasses
from .vienna import cds_to_rna, average_unpaired, ensemble_free_energy
from .protein import CodonFrequencyTable


@dataclasses.dataclass
class CAIThresholdObjectiveConfig:
    """Configuration for CAI threshold objective function."""
    freq_table: CodonFrequencyTable
    cai_threshold: float = 0.8
    cai_exp_scale: float = 1.0
    verbose: bool = False


def make_cai_threshold_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI up to the threshold"""
    def obj(cds: List[str]) -> float:
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)-1
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}")
        fitness = -cai_penalty
        measures = {"CAI": cai}
        return fitness, measures
    return obj


def make_cai_and_aup_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI and AUP: (1-aup)-(e^(max(0,threshold-cai)*scale)-1)"""
    def obj(cds: List[str]) -> float:
        rna_seq = cds_to_rna(cds)
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)-1
        aup = average_unpaired(rna_seq)
        ensemble_paired_prob = 1-aup
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}, AUP: {aup}")
        fitness = ensemble_paired_prob - cai_penalty
        measures = {"CAI": cai, "AUP": aup}
        return fitness, measures
    return obj


def make_cai_and_efe_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI and EFE: -efe*(1/e^(max(0,threshold-cai)*scale))"""
    def obj(cds: List[str]) -> float:
        rna_seq = cds_to_rna(cds)
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)
        efe = ensemble_free_energy(rna_seq)
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}, EFE: {efe}")
        fitness = -efe*(1/cai_penalty)
        measures = {"CAI": cai, "EFE": efe}
        return fitness, measures
    return obj
