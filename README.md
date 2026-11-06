# Analysis of Nicking Sites in Circular Plasmids

Map sequencing fragments to a circular reference plasmid in both orientations, characterize terminal nucleotide composition, and report fragment mapping coordinates. Outputs include CSV tables, PNG figures, and editable SVG figures.

Only two FASTA inputs are required: a fragment file and a reference file containing one forward-strand plasmid sequence. The reference is supplied as an ordinary FASTA; circular wraparound and reverse-complement matching are handled automatically.

## Features

- Full-length circular fragment mapping, including matches spanning the reference origin.
- Forward and reverse-complement searches, retaining every position that meets the identity threshold.
- Mapping coordinates, strand, fragment length, identity, and candidate cleavage boundaries.
- A combined FASTA of matched fragments and separate FASTA files for their first and last six bases.
- Separate first-six and last-six frequency logos, information-content logos, and nucleotide-frequency bar charts, with plotting data.
- Explicit reporting of unmatched fragments and unsupported sequence characters.

## Repository layout

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── circular_fragments.py
├── plot_figures.py
└── examples/
    ├── fragments.fasta
    └── plasmid.fasta
```


`circular_fragments.py` runs the complete workflow and requires `plot_figures.py` in the same directory. The plotting script can also regenerate figures from existing tables without repeating mapping.

## Installation

Python 3.10 or later is required.

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Dependencies: Biopython, NumPy, pandas, and Matplotlib. Plotting uses a noninteractive backend and does not require a graphical desktop.

## Quick start

Run the synthetic example from the repository root:

```bash
python circular_fragments.py \
  --fragments examples/fragments.fasta \
  --reference examples/plasmid.fasta \
  --output results/example
```

Analyze your own data:

```bash
python circular_fragments.py \
  --fragments /path/to/fragments.fasta \
  --reference /path/to/plasmid.fasta \
  --output results/sample \
  --identity 0.95
```

| Argument | Description | Default |
|---|---|---|
| `--fragments` | Required fragment FASTA; record IDs must be unique | Required |
| `--reference` | Required FASTA containing exactly one forward reference sequence | Required |
| `--output` | Output directory | `circular_results` |
| `--identity` | Full-length nucleotide identity threshold, from 0 to 1 | `0.95` |

Progress is printed during processing. Existing output files with the same names are overwritten; use a separate output directory for each sample.


Regenerate figures from existing output tables:

```bash
python plot_figures.py \
  --output results/sample
```

This command requires `end6_counts.csv`, `end6_frequency.csv`, `end6_information_bits.csv` in the output directory. No reference FASTA is required for replotting.

Python or Jupyter usage:

```python
from circular_fragments import analyze

summary = analyze(
    'fragments.fasta', 'plasmid.fasta', 'results/sample', threshold=0.95
)
```

## Mapping method

The workflow performs full-query, ungapped circular-window matching. Identity is the number of matching bases divided by fragment length, with a default threshold of 95%. Substitutions are allowed; insertions and deletions are not. This is neither local alignment nor gapped PairwiseAligner alignment.

To accelerate candidate detection, a query of length m is divided into k+1 nonoverlapping seeds, where k is the maximum number of mismatches permitted by the threshold. Any qualifying window must contain at least one exact seed. Candidate windows are then checked for full-length identity. Both the original query and its reverse complement are searched against the circular reference, with all coordinates reported on the supplied forward reference.

All qualifying sites and orientations are retained. Duplicating the reference internally enables origin-spanning windows without reporting the same reference start twice.

## Main outputs

| File | Contents |
|---|---|
| `fragment_matches.csv` | Fragment ID, coordinates, strand, length, identity, origin-spanning flag, and boundary positions for each qualifying mapping |
| `matched_fragments_both_strands.fasta` | Matched original fragments, each input record written once in its original orientation |
| `first6.fasta`, `last6.fasta` | First and last six bases of matched fragments |
| `ends6_combined.fasta` | Combined terminal sequences with end-specific ID suffixes |
| `first6_frequency_logo.*`, `last6_frequency_logo.*` | Separate terminal nucleotide-frequency logos |
| `first6_information_bits_logo.*`, `last6_information_bits_logo.*` | Separate information-content logos; each y-axis extends 10% above its maximum stack height |
| `first6_frequency_bars.*`, `last6_frequency_bars.*` | Separate A/C/G/T frequency bar charts |
| `first6_counts.csv`, `last6_counts.csv` | Nucleotide counts at each terminal position |
| `first6_frequency.csv`, `last6_frequency.csv` | Nucleotide frequencies at each terminal position |
| `first6_information_bits.csv`, `last6_information_bits.csv` | Per-base information-content heights |
| `unmatched_fragments.csv` | Unmatched fragment IDs, lengths, and reasons |
| `summary.json` | Input counts, mapping counts, parameters, and counting conventions |

`.*` denotes PNG and SVG versions. The complete workflow also retains combined `end6_*` figures and matrices as supplementary outputs; use the `first6_*` and `last6_*` figures for separate terminal views.

## Coordinate conventions

- All base coordinates refer to the supplied forward reference and use **1-based, inclusive** intervals.
- `Start_bp` and `End_bp` describe the covered interval along the forward reference. Start exceeds End when a mapping crosses the reference origin.
- `Strand` is `+` or `-`. For a reverse-strand fragment, its original 5′ end corresponds to End and its 3′ end corresponds to Start. `Fragment_5prime_bp` and `Fragment_3prime_bp` explicitly report those terminal coordinates.
- `Cut_after_bp=k` denotes the boundary between forward-reference bases k and k+1. **Zero denotes the boundary between the final plasmid base and base 1.**
- These are candidate cleavage positions inferred from fragment boundaries. No enzyme-specific offset is assumed, and the positions are not independently validated cleavage events.

## Counting conventions

Terminal nucleotide distributions count each matched input fragment once in its original 5′→3′ orientation. Multiple mappings do not increase that fragment's weight in terminal statistics. First-six positions are labeled 1 through 6; last-six positions are labeled -6 through -1. Information-content heights are calculated as `p × (2 − H)` using a uniform DNA background without small-sample correction.



## Input constraints

- The reference must contain exactly one nonempty A/C/G/T DNA sequence.
- Fragment IDs must be unique; duplicate IDs cause an error.
- Non-A/C/G/T or empty fragments are reported as unmatched with reason `non_ACGT_or_empty` and are excluded from mapping.
- Fragments longer than one plasmid circumference are reported as `longer_than_plasmid`; multicycle fragments are unsupported.
- Matched fragments shorter than six bases remain in mapping results but are excluded from terminal-six-base analysis. Terminal windows may overlap for fragments of 6–11 bases.
- If no fragments match, the workflow still produces empty tables, empty FASTA files, and zero-height plots.
- Very low identity thresholds, short queries, or highly repetitive references may generate many candidates and reduce performance.

## Validation

Synthetic checks covered forward/reverse origin-spanning coordinates, the 95% identity threshold, repeated sites, palindromic matches in both orientations, oversized queries, short fragments, unmatched fragments, and terminal count/frequency totals. Files in `examples/` are synthetic demonstration data, not experimental sequencing data.

An internal run on 50,348 fragments and a 2,961 bp plasmid produced 41,293 matched fragments: 21,024 forward and 20,269 reverse, including 2,841 origin-spanning mappings. Another 8,641 fragments were below threshold and 414 contained non-ACGT characters. Coordinate, threshold, matched-FASTA count, checks passed. These experimental inputs are not distributed with this repository, so this result is not an independently reproducible public benchmark.

## License

MIT License
