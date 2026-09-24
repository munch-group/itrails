"""
Simulate a dummy four-species MAF alignment for the iTRAILS workflow.

One haploid genome is sampled per species on the tree
(((hg38, panTro5), gorGor5), ponAbe2) with split times and ancestral
population sizes mirroring human/chimp/gorilla/orangutan (generation
time 25 years):

    human-chimp split      240,000 generations (~6 Mya),   N_AB   =  65,000
    +gorilla               350,000 generations (~8.75 Mya), N_ABC = 100,000
    +orangutan             700,000 generations (~17.5 Mya), N_root = 100,000

Extant population sizes do not affect the data (a single lineage per
species cannot coalesce before the splits) but are set to plausible
values anyway. The mutation rate matches the fixed `mu` in
config/example_config.yaml rather than the literature human rate, since
all iTRAILS estimates are scaled relative to that fixed value.

Regenerate with:  pixi run simulate-example
Output:           data/simulated_alignment.maf (1 Mb in 50 kb MAF blocks)
"""

from pathlib import Path

import msprime
import numpy as np
import argparse

def simulate():
    demography = msprime.Demography()
    demography.add_population(name=args.SPECIES_A, initial_size=args.N_A)
    demography.add_population(name=args.SPECIES_B, initial_size=args.N_B)
    demography.add_population(name=args.SPECIES_C, initial_size=args.N_C)
    demography.add_population(name=args.SPECIES_D, initial_size=args.N_D)
    demography.add_population(name='AB', initial_size=args.N_AB)
    demography.add_population(name='ABC', initial_size=args.N_ABC)
    demography.add_population(name='ABCD', initial_size=args.N_ABCD)
    demography.add_population_split(time=args.T_AB, derived=[args.SPECIES_A, args.SPECIES_B], ancestral='AB')
    demography.add_population_split(time=args.T_ABC, derived=['AB', args.SPECIES_C], ancestral='ABC')
    demography.add_population_split(time=args.T_ABCD, derived=['ABC', args.SPECIES_D], ancestral='ABCD')

    ts = msprime.sim_ancestry(
        samples={species: 1 for species in SPECIES},
        demography=demography,
        ploidy=1,
        sequence_length=args.SEQ_LEN,
        recombination_rate=args.REC,
        random_seed=args.SEED,
    )
    return msprime.sim_mutations(ts, rate=args.MU, random_seed=args.SEED + 1)


def alignment(ts):
    """Return one sequence per species: a shared random ancestral background
    with the simulated alleles written in at variant sites."""
    rng = np.random.default_rng(args.SEED + 2)
    ancestral = rng.choice(np.array(list('ACGT')), size=args.SEQ_LEN)
    seqs = {species: ancestral.copy() for species in SPECIES}

    sample_order = [ts.population(ts.node(u).population).metadata['name']
                    for u in ts.samples()]
    n_variants = 0
    for variant in ts.variants():
        pos = int(variant.site.position)
        for species, genotype in zip(sample_order, variant.genotypes):
            seqs[species][pos] = variant.alleles[genotype]
        n_variants += 1
    print(f'{n_variants} variant sites in {args.SEQ_LEN} bp')
    return seqs


def write_maf(seqs, path):
    with open(path, 'w') as f:
        f.write('##maf version=1 scoring=none\n')
        f.write('# simulated with msprime by scripts/simulate_data.py\n')
        for start in range(0, args.SEQ_LEN, args.BLOCK_LEN):
            end = min(start + args.BLOCK_LEN, args.SEQ_LEN)
            f.write('\na score=0\n')
            for species in SPECIES:
                block = ''.join(seqs[species][start:end])
                f.write(f's {species}.chr1 {start} {end - start} + {args.SEQ_LEN} {block}\n')
        f.write('\n')
    print(f'wrote {path}')


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Process a user name.")
    parser.add_argument("output", help="Path of output maf file")
    parser.add_argument("--seed", dest="SEED", type=int)
    parser.add_argument("--seq-len", dest="SEQ_LEN", type=int)
    parser.add_argument("--block-len", dest="BLOCK_LEN", type=int)
    parser.add_argument("--mu", dest="MU", type=float)
    parser.add_argument("--rec", dest="REC", type=float)
    parser.add_argument("--species-a", dest="SPECIES_A")
    parser.add_argument("--species-b", dest="SPECIES_B")
    parser.add_argument("--species-c", dest="SPECIES_C")
    parser.add_argument("--species-d", dest="SPECIES_D")
    parser.add_argument("--t-ab", dest="T_AB", type=int)
    parser.add_argument("--t-abc", dest="T_ABC", type=int)
    parser.add_argument("--t-abcd", dest="T_ABCD", type=int)
    parser.add_argument("--n-a", dest="N_A", type=int)
    parser.add_argument("--n-b", dest="N_B", type=int)
    parser.add_argument("--n-c", dest="N_C", type=int)
    parser.add_argument("--n-d", dest="N_D", type=int)
    parser.add_argument("--n-ab", dest="N_AB", type=int)
    parser.add_argument("--n-abc", dest="N_ABC", type=int)
    parser.add_argument("--n-abcd", dest="N_ABCD", type=int)
    args = parser.parse_args()

    SPECIES = [args.SPECIES_A, args.SPECIES_B, args.SPECIES_C, args.SPECIES_D]

    assert args.output.endswith('maf')

    OUT_PATH = Path(args.output)
    assert OUT_PATH.exists

    write_maf(alignment(simulate()), OUT_PATH)


    
"""
pixi run python scripts/simulate_data.py --seed 42\
    --seq-len 10000000 --block-len 250000 \
    --mu 1.2e-8 --rec 1e-9 \
    --species-a hg38 --n-a 15_000 \
    --species-b panTro5 --n-b 35_000 \
    --species-c gorGor5 --n-c 25_000 \
    --species-d ponAbe2 --n-d 25_000 \
    --t-ab 220000 --t-abc 340000 --t-abcd 1000000 \
    --n-ab 86000 --n-abc 67000 --n-abcd 167000 \
    data/hcg.maf
"""



    # SEED = args.seed
    # SEQ_LEN = args.seq_len
    # BLOCK_LEN = args.block_len
    # MU = args.mu    # per site per generation; must match the config's fixed mu
    # RHO = args.rho   # per site per generation

    # ##### RHO NOT R... OR WHAT? 


    # species order: A, B, C, outgroup

    # N_EXTANT = dict(zip(SPECIES, [args.N_A, args.N_B, args.N_C, args.N_D]))
#               'hg38', 'panTro5', 'gorGor5', 'ponAbe2']

    # T_AB = 220000     # human-chimp split (generations)
    # T_ABC = 340_000    # gorilla split
    # T_ROOT = 1_000_000   # orangutan split

    # # N_EXTANT = {'hg38': 15_000, 'panTro5': 35_000, 'gorGor5': 25_000, 'ponAbe2': 25_000}

    # N_AB = 86_000
    # N_ABC = 67_000
    # N_ROOT = 167_000

    # OUT_PATH = Path(__file__).parent.parent / 'data' / 'simulated_alignment.maf'

    # write_maf(alignment(simulate()), OUT_PATH)



# SEED = 42
# SEQ_LEN = 10_000_000
# BLOCK_LEN = 250_000
# MU = 1.2e-8    # per site per generation; must match the config's fixed mu
# RHO = 1e-8   # per site per generation

# # species order: A, B, C, outgroup
# SPECIES = ['hg38', 'panTro5', 'gorGor5', 'ponAbe2']

# T_AB = 220_000     # human-chimp split (generations)
# T_ABC = 340_000    # gorilla split
# T_ROOT = 1_000_000   # orangutan split

# N_EXTANT = {'hg38': 15_000, 'panTro5': 35_000, 'gorGor5': 25_000, 'ponAbe2': 25_000}
# N_AB = 86_000
# N_ABC = 67_000
# N_ROOT = 167_000

# OUT_PATH = Path(__file__).parent.parent / 'data' / 'simulated_alignment.maf'
