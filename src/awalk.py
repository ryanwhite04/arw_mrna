""" Implements adaptive random walks on CDSs. """
from typing import Callable, List, Optional
import protein
import random
import dataclasses


@dataclasses.dataclass
class WalkConfig:
    aa_seq: str
    freq_table: protein.CodonFrequencyTable
    objective: Callable[[List[str]], float]
    steps: int
    init_cds: Optional[List[str]] = None
    verbose: bool = False


@dataclasses.dataclass
class WalkResult:
    cds: List[str]
    fitness: float


def adaptive_random_walk(config: WalkConfig):
    """Runs ARW to maximize objective function on CDSs."""
    cds = protein.random_cds(
        config.aa_seq, config.freq_table) if config.init_cds is None else config.init_cds
    prev_fitness = config.objective(cds)
    if config.verbose:
        print(f"Initial CDS: {cds}")
    for step in range(config.steps):
        mutidx = random.randint(0, len(cds) - 1)
        while len(config.freq_table.get_codons(config.aa_seq[mutidx])) == 1:
            mutidx = random.randint(0, len(cds) - 1)
        mutcodons = list(config.freq_table.get_codons(
            config.aa_seq[mutidx])-set([cds[mutidx]]))
        mutcodon = random.choice(mutcodons)
        new_cds = cds[:mutidx] + [mutcodon] + cds[mutidx+1:]
        new_fitness = config.objective(new_cds)
        if new_fitness > prev_fitness:
            cds = new_cds
            prev_fitness = new_fitness
            if config.verbose:
                print(f"New CDS: {new_cds}")
        if config.verbose:
            print(
                f"Step: {step}, Fitness: {new_fitness}, Best Fitness: {prev_fitness}")
    return WalkResult(cds, prev_fitness)


def adaptive_random_walk_generator(config: WalkConfig):
    """Runs ARW to maximize objective function on CDSs."""
    cds = protein.random_cds(
        config.aa_seq, config.freq_table) if config.init_cds is None else config.init_cds
    prev_fitness = config.objective(cds)
    if config.verbose:
        yield {"type": "initial", "cds": cds, "fitness": prev_fitness}
    for step in range(config.steps):
        mutidx = random.randint(0, len(cds) - 1)
        while len(config.freq_table.get_codons(config.aa_seq[mutidx])) == 1:
            mutidx = random.randint(0, len(cds) - 1)
        mutcodons = list(config.freq_table.get_codons(
            config.aa_seq[mutidx]) - set([cds[mutidx]]))
        mutcodon = random.choice(mutcodons)
        new_cds = cds[:mutidx] + [mutcodon] + cds[mutidx+1:]
        new_fitness = config.objective(new_cds)
        if new_fitness > prev_fitness:
            cds = new_cds
            prev_fitness = new_fitness
            if config.verbose:
                yield {"type": "new_cds", "cds": new_cds, "fitness": new_fitness, "step": step, "best_fitness": prev_fitness}
        if config.verbose:
            yield {"type": "progress", "step": step, "fitness": new_fitness, "best_fitness": prev_fitness}
    yield {"type": "final", "cds": cds, "fitness": prev_fitness}