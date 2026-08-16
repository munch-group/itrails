
# iTRAILS workflow component

A [GWF](https://gwf.app/) workflow around [iTRAILS](https://itrails.readthedocs.io/en/docs-stable/), which fits the TRAILS coalescent hidden Markov model ([Rivas-González et al., 2024](https://doi.org/10.1371/journal.pgen.1010836)) to a four-genome MAF alignment (species A, B, C and an outgroup) and decodes gene-tree topologies along the genome.

For each analysis the workflow chains three targets, each writing to its own subfolder of the analysis folder:

1. `itrails-optimize` — maximum-likelihood fit of ancestral population sizes and split times → `optimize/{name}.best_model.yaml`
2. `itrails-viterbi` — most likely hidden-state sequence along the alignment → `viterbi/{name}.viterbi.csv`
3. `itrails-posterior` — per-position posterior state probabilities → `posterior/{name}.posterior.csv`

The repository works both **standalone** and as a **git submodule** component of a larger gwf project. iTRAILS itself runs in an isolated pixi environment defined in `pixi.toml` (it pins numpy/scipy/numba), so nothing needs to be installed beyond pixi.

📖 **Documentation:** [munch-group.github.io/itrails](https://munch-group.github.io/itrails/) — rendered from the `docs/` pages with Quarto and published on every push to `main`. Build it locally with `pixi run quarto render`.

## Standalone use

The repo ships a 1 Mb simulated great-ape-like alignment (`data/simulated_alignment.maf`, made with msprime by `scripts/simulate_data.py`; regenerate with `pixi run simulate-example`), and `analyses.yml` points at it out of the box. To use the real great-ape alignment instead (183 MB, hg38/panTro5/gorGor5/ponAbe2), download it with `pixi run download-example` and set `maf: data/example_alignment.maf` in `analyses.yml`.

Edit `analyses.yml` (analyses, output dir, Slurm account) and the per-analysis
optimization config (see `config/example_config.yaml`), then:

```bash
pixi run gwf run          # on a Slurm cluster: gwf config set backend slurm
pixi run gwf status
```

To test locally, start a worker pool in a second terminal with `pixi run gwf -b local workers` and use `gwf -b local run`.

Outputs land in `steps/itrails/{name}/`, one subfolder per step (`split/`, `optimize/`, `viterbi/`, `posterior/`, `concat/`), and `notebooks/example.ipynb` (executed as the last workflow step) gives a first look at them.

## Use as a submodule

Add this repository to the parent project:

```bash
git submodule add git@github.com:munch-group/itrails.git itrails
```

In the parent `workflow.py`, import the factory and graft the targets onto the parent's own `gwf` object:

```python
import os, sys
from gwf import Workflow

gwf = Workflow(defaults={'account': 'my-project'})

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'itrails'))
from itrails_workflow import itrails_workflow

gwf, itrails_targets = itrails_workflow(
    gwf=gwf,
    analyses=[
        dict(name='primates',
             maf='data/alignment.maf',
             config='config/itrails_config.yaml'),
    ],
    output_dir='steps/itrails',
)

# downstream parent targets can depend on e.g.:
#   itrails_targets['viterbi'][0].outputs['viterbi']
#   itrails_targets['optimize'][0].outputs['best_model']
```

Paths are resolved relative to where `gwf run` is invoked (the parent project root). Resource options can be overridden globally (`optimize_options={'cores': 64, 'memory': '128g'}`) or per analysis (an `optimize_options` key in the analysis dict); the same goes for `viterbi_options` and `posterior_options`.

By default each target runs its command inside this repo's pixi `itrails` environment (`pixi run -e itrails ...` with the submodule's manifest), which works on both the local and slurm backends. Pass `run_prefix=...` to use another environment instead.

## Windowed analyses

For a genome-scale alignment, set `window_size` (per analysis in `analyses.yml`, or as an argument to `itrails_workflow()`), e.g.:

```yaml
analyses:
  - name: primates
    maf: data/alignment.maf
    config: config/itrails_config.yaml
    window_size: 2000000
    species: [hg38, panTro5, gorGor5, ponAbe2]
```

The alignment is then split into 2 Mb windows at MAF-block boundaries (lossless — iTRAILS treats each block as an independent sequence), every window is decoded as its own job, and the per-window CSVs are concatenated into `concat/{name}.viterbi.csv` / `concat/{name}.posterior.csv`, each row prefixed with `window`, `chrom`, `block_start` and `src_block` columns that place it on the genome. The window MAFs and the `windows.tsv` / `blocks.tsv` manifests live in `split/`, and the per-window decode files in `viterbi/` and `posterior/`.

- **`fit: genome`** (default) fits the model once on the whole alignment and decodes every window with it — the standard coalescent-HMM setup.
- **`fit: window`** fits a separate model per window (one heavy `itrails-optimize` job each) for studying parameter variation along the genome; expect noisy estimates from small windows.

Notes:

- Windows are binned by the start coordinate of the `reference` species' line in each block (first line of the block if `reference` is unset). Blocks are never split, so window edges are approximate at block granularity.
- Set `species` to the iTRAILS `species_list` so the splitter drops the same incomplete blocks the iTRAILS parser would skip; otherwise the block-index annotation in the concatenated CSVs assumes every block parses.
- Window boundaries are derived from the MAF **when the workflow is defined** (the scan is cached by file size/mtime), so the MAF must already exist — if an upstream target produces it, run that part first.
- `viterbi_options`/`posterior_options` (and `optimize_options` with `fit: window`) apply to *each* per-window job, so size them for a single window. The concatenated posterior CSV has one row per alignment position and can get very large for a whole genome.

## Starting from phased VCFs instead of an alignment

For closely related species that were all mapped and variant-called against the same reference genome — e.g. baboon species against the macaque genome — the shared coordinate system makes the genomes implicitly aligned, and the workflow can build the four-species alignment itself from phased (g)VCFs. Give `vcf`, `fasta` and `samples` instead of `maf`:

```yaml
analyses:
  - name: baboons
    vcf: [data/papio.g.vcf.gz, data/theropithecus.g.vcf.gz]  # or a single path
    fasta: data/rheMac10.fa           # reference the VCFs were called against
    samples: [pAnubis1, pHamadryas1, pCynocephalus1, REF]    # A, B, C, outgroup
    config: config/baboons_config.yaml
    window_size: 2000000
```

An `alignment/` step (`itrails_vcf2maf_{name}`, running `scripts/vcf_to_maf.py`) reconstructs the **first haplotype** of each sample (`haplotype: 2` selects the second) by writing its phased variant alleles onto the reference sequence, and writes `alignment/{name}.maf` in 1 Mb blocks (`block_size`). Everything downstream — optimize, windowing, decode, concat — is unchanged.

**The `fasta` key is optional.** Without it, sequences are reconstructed from the VCF records alone — meant for all-sites *genome* VCFs, where every callable position has a record carrying its REF base. Any position without a record for a sample is written as `N`, and the `REF` pseudo-sample is assembled from the REF fields of the records. Contig names, lengths and order then come from the `##contig` header lines (lengths required). One caveat: banded gVCF non-variant blocks (`END=` records) only carry the base at their anchor position, so with such files the block interiors come out as `N` — the converter warns about this, and supplying `fasta` fills them from the reference.

Notes:

- **Samples double as species names**: use the sample IDs in the config's `species_list` (IDs must not contain `.`). The reserved name `REF` stands for the reference sequence itself — handy when the reference genome *is* the outgroup, as with macaque for baboons.
- Samples may be spread over multiple VCFs (per-sample or per-region files); each file contributes the samples in its header. Files must be coordinate-sorted in the same contig order as the FASTA (or the first file's header).
- **gVCF semantics by default**: positions not covered by any record (and any non-SNP or missing call, which reference coordinates cannot represent) become `N`. For plain VCFs where "no record" means reference-equal, set `mask_uncovered: false` (needs `fasta`). `FILTER` is ignored — pre-filter if needed.
- `regions: [chr1, chr2]` restricts the analysis to selected contigs (default: every contig in the FASTA or VCF headers).
- Windowed VCF analyses compute window boundaries from contig lengths, so the alignment does not have to exist when the workflow is defined: with `fasta` from the index `{fasta}.fai` (`samtools faidx`), without it from the VCF `##contig` lines (the VCFs must then exist when the workflow is defined).

## Reference notes

- `claude-itrails-ref.md` — iTRAILS CLI, config format, model parameters, and output files
- `claude-gwf-ref.md` / `gwf.md` — GWF API, CLI, and patterns

## Initial set up

```bash
pixi run init
```

## Incorporate updates to upstream fork

If you forked `itrails` rather than using it as template, you can incorporate changes/fixes made to `itrails`.

Add upstream if not already added

```bash
git remote add upstream https://github.com/munch-group/itrails.git
```

Fetch upstream changes

```bash
git fetch upstream
```

Either rebase your changes on top of upstream (cleaner history)

```bash
git rebase upstream/main
```

Or, merge upstream into your fork (preserves history)

```bash
git merge upstream/main
```

If you want to see what's changed upstream before applying:

```bash
git log HEAD..upstream/main
```

See the actual diff

```bash
git diff HEAD...upstream/main
```

Then push your updated fork:

```bash
git push origin main
```

If you rebased and need to force push

```bash
git push origin main --force-with-lease
```
