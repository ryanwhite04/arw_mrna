""" Implements adaptive random walks on CDSs. """
from typing import Callable, List, Optional, Generator
from .protein import CodonFrequencyTable, random_cds
import random
from dataclasses import dataclass, asdict

@dataclass
class WalkConfig:
    aa_seq: str
    freq_table: CodonFrequencyTable
    objective: Callable[[List[str]], float]
    steps: int
    init_cds: Optional[List[str]] = None,
    verbose: bool = False
    seed: Optional[int] = None

@dataclass
class Measures:
    CAI: float
    AUP: Optional[float] = None
    EFE: Optional[float] = None

@dataclass
class WalkResult:
    cds: List[str]
    step: int
    measures: Measures
    fitness: float
    best_fitness: float
    
    def to_dict(self):
        return asdict(self)

def adaptive_random_walk(config: WalkConfig) -> Generator[WalkResult, None, None]:
    """Runs ARW to maximize objective function on CDSs."""
    if config.seed is not None:
        random.seed(config.seed)
    cds = config.init_cds
    if cds is None:
        cds = random_cds(config.aa_seq, config.freq_table, config.seed)
    # cds = random_cds(config.aa_seq, config.freq_table, config.seed) if config.init_cds is None else config.init_cds,
    best_fitness, measures = config.objective(cds)
    if config.verbose:
        print(f"Initial CDS: {cds}")
    yield WalkResult(cds, 0, measures, best_fitness, best_fitness)
    for step in range(config.steps):
        mutidx = random.randint(0, len(cds) - 1)
        while len(config.freq_table.get_codons(config.aa_seq[mutidx])) == 1:
            mutidx = random.randint(0, len(cds) - 1)
        mutcodons = list(config.freq_table.get_codons(config.aa_seq[mutidx]) - set([cds[mutidx]]))
        mutcodon = random.choice(mutcodons)
        new_cds = cds[:mutidx] + [mutcodon] + cds[mutidx+1:]
        fitness, measures = config.objective(new_cds)
        if fitness > best_fitness:
            cds = new_cds
            best_fitness = fitness
            if config.verbose:
                print(f"New CDS: {new_cds}")
        if config.verbose:
            print(f"Step: {step}, Fitness: {fitness}, Best Fitness: {best_fitness}")
        yield WalkResult(cds, step+1, measures, fitness, best_fitness)
