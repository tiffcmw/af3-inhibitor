# Systemic Discovery of Phage-Encoded Inhibitors in Oral Pathogens

A computational pipeline for discovering inhibitors of **gingipains**, the cysteine protease
virulence factors (RgpA, RgpB, Kgp) produced by *Porphyromonas gingivalis*, a key pathogen in
periodontal (gum) disease. The pipeline builds a candidate library of propeptide-derived and
protease-inhibitor sequences, filters/clusters it down with GPU-accelerated tools, and uses
AlphaFold3 / ColabFold complex prediction to rank candidates by how confidently they're predicted
to bind and block the gingipain active site.

All of the code lives in a single notebook: [`code/code.ipynb`](code/code.ipynb).

## Background

Gingipains are synthesized with an attached **propeptide** region that folds back and blocks their
own active site until the enzyme is secreted and activated. Recombinant forms of this propeptide
have been shown to bind to and inhibit the active enzyme. This suggests a general strategy for this
class of target: pathogens whose virulence depends on a protease that is natively autoinhibited by
its own propeptide are candidates for inhibition by recombinant propeptide-derived peptides (or
other peptides that competitively bind the same active site).

## Workflow

```
0) Collect phage/candidate protein library
1) Generate + curate a candidate library (sequence clustering & filtering, GPU-accelerated)
2) AlphaFold3 complex prediction on a high-confidence shortlist (expensive, accurate)
3) Fast docking / learned scoring for broad screening (cheap, approximate)
4) Feed docking hits back into AlphaFold3 for refined complex prediction
```

Concretely, the notebook implements this as four stages:

1. **Stage 1: Target Justification.** Establishes which defense systems, in which oral pathogens,
   are worth targeting. Worked example: *P. gingivalis* gingipains and their autoinhibitory
   propeptides.
2. **Stage 2: Building the Phage Protein Library.** Two independent tracks that each produce a
   FASTA candidate pool:
   - *Variant track*: extract gingipain propeptide sequences from UniProt, generate all
     single-point-substitution variants.
   - *Database track*: pull candidate inhibitor-like sequences directly from UniProt (cysteine
     protease inhibitor keyword, mechanism-matched to gingipains) and from MEROPS.
3. **Stage 3: Data Filtering.** Cuts the pool down before the expensive folding step, using
   BLOSUM62 comparative scoring, GPU-vectorized (cuDF/cuPy) biophysical filtering, top-N
   diversity-aware clustering, and MMseqs2-GPU clustering. Also extracts the gingipain catalytic
   domains that serve as the fold targets.
4. **Stage 4: Complex Prediction.** Predicts whether each surviving candidate binds the gingipain
   catalytic domain, via **AlphaFold3** (job generation + multi-GPU launch script for a remote
   server) or **ColabFold** (`colabfold_batch`, cheaper/faster, with built-in residue-numbering
   verification, active-site-overlap checking, results parsing, and an inline results viewer).
   Candidates are ranked by ipTM/pTM; top hits, ideally scoring better than the native propeptide
   itself, are the pipeline's discovered inhibitor candidates.

The two Stage 2 tracks and two Stage 4 backends are independent alternatives, not a strict chain.
You can, for example, run the database track straight through ColabFold without touching
AlphaFold3.

## Repository layout

```
af3-inhibitor/
├── code/
│   ├── code.ipynb        # the entire pipeline, organized by stage (see pipeline.md for details)
│   ├── images/           # figures referenced from the notebook
│   └── output/           # generated intermediates & results (gitignored, created by running cells)
└── datasets/
    ├── UniProt/          # canonical RgpA/RgpB/Kgp FASTA sequences (source of truth for
    │                     # propeptide/chain boundary numbering)
    └── MEROPS/           # MEROPS family C25 (clan CD) peptidase/inhibitor library
```

## Setup

Recommended: a Python virtual environment on a RAPIDS-compatible Python version (3.11–3.14).

```bash
# create & activate a venv
python3.11 -m venv myvenv
source myvenv/bin/activate   # macOS/Linux

# check your CUDA version (needed to pick a matching RAPIDS build)
nvcc --version
```

Then install, matched to your CUDA/Python version:

- **[RAPIDS](https://docs.rapids.ai/install/)** (cuDF/cuPy): used for GPU-accelerated filtering in
  Stage 3. Falls back to CPU/pandas automatically if unavailable.
- **[Docker + NVIDIA Docker support](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#installing-docker)**
  and **[AlphaFold3](https://github.com/google-deepmind/alphafold3)** (model weights + databases): required only if using the AlphaFold3 backend in Stage 4a.
- **[ColabFold](https://github.com/YoshitakaMo/localcolabfold)** (`colabfold_batch` on `PATH`):
  required only if using the ColabFold backend in Stage 4b.
- **MMseqs2** with GPU support, on `PATH`: required for Stage 3bii clustering:
  ```bash
  conda install -c bioconda mmseqs2
  ```

## Usage

Open `code/code.ipynb` and run the cells stage by stage. Each stage's CONFIG cell documents the
paths/parameters you need to set before running it (e.g. UniProt boundary maps, AlphaFold3
server paths, ColabFold extra args, MMseqs2 clustering thresholds). Outputs from each stage are
written under `code/output/` and consumed by later stages.

AlphaFold3 and ColabFold jobs are generated locally but intended to run on a remote GPU server.
the corresponding notebook cells write job files/launch scripts and print the command to run
there rather than executing the fold jobs in-notebook.

## Data sources

- **UniProt**: canonical (reviewed) RgpA, RgpB, and Kgp sequences for *P. gingivalis*.
- **[MEROPS](https://www.ebi.ac.uk/merops/download_list.shtml)**: family C25 (clan CD,
  gingipain/caspase-like cysteine proteases) peptidase and inhibitor unit library.

## References

- Interaction of gingipain catalytic domains with their propeptides:
  [PLOS ONE, 10.1371/journal.pone.0065447](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0065447)
- Propeptide-mediated gingipain inhibition:
  [PMC3677877](https://pmc.ncbi.nlm.nih.gov/articles/PMC3677877/),
  [PMC10403534](https://pmc.ncbi.nlm.nih.gov/articles/PMC10403534/)

## For contributors

See [`pipeline.md`](pipeline.md) for a deeper dive into the pipeline's internal architecture,
hardcoded boundary maps, and known gotchas when modifying the notebook.
