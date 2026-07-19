# pyISFinder 🧬

**pyISFinder** is a fast, standalone, pure-Python pipeline for Insertion Sequence (IS) element identification, annotation, and visualization in bacterial genomes.

Unlike legacy C-based or Python 2 IS prediction tools, `pyISFinder` is written natively for modern Python 3.10+ environments using high-performance C-extensions (`pyrodigal` and `pyhmmer`), avoiding binary compilation errors and dependency conflicts.

---

## 🌟 Key Features

* **In-Memory Bacterial Gene Prediction**: Powered by `pyrodigal` (Prodigal Cython bindings).
* **Transposase Profile HMM Search**: Fast domain searching via `pyhmmer`.
* **Repeat Detection Engine**: Searches for flanking Terminal Inverted Repeats (TIRs) and Direct Repeats / Target Site Duplications (DRs/TSDs).
* **IS Family Classification**: Classifies IS elements into canonical families (`IS1`, `IS3`, `IS4`, `IS5`, `IS6`, `IS21`, `IS30`, `IS110`, `IS256`, `IS607`, `IS630`, `IS66`, `IS982`, `IS1380`, `ISL3`, etc.).
* **Multi-Format Export**: Generates `.tsv` prediction summaries, standard `.gff3` tracks for IGV/JBrowse/Artemis, and manuscript summary statistics.
* **Publication Figures**: Automatically outputs high-resolution linear genomic maps (`.png` & `.svg`) using `matplotlib`.

---

## 🚀 Quick Start

### 1. Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/MdUmar-tech/pyISFinder.git
cd pyISFinder

# Option A: Install via pip
pip install -r requirements.txt

# Option B: Create virtual environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2. Usage

Run `pyISFinder` on any bacterial genome FASTA file:

```bash
python3 pyISFinder.py -i path/to/genome.fasta -o output_directory
```

#### Optional Arguments:
* `-i`, `--input`: Input bacterial genome FASTA file (required).
* `-o`, `--outdir`: Output directory for generated reports and plots (default: `results`).
* `--hmm`: Path to custom transposase profile HMM file (default: `db/clusters.faa.hmm`).

---

## 📊 Output Files

`pyISFinder` produces 4 main output files in the output directory:

1. **`*predictions.tsv`**: Detailed tabular report containing IS ID, Scaffold, IS Family, Start, End, Strand, Length, ORF IDs, TIR Left/Right sequences, TIR Identity %, and Score.
2. **`*.gff3`**: Standard GFF3 browser track for visualization in IGV, JBrowse, or Artemis.
3. **`*_summary.txt`**: Summary statistics table detailing counts and percentages per IS family.
4. **`*_genome_map.png` / `.svg`**: Vector graphic genomic track map showing IS element distributions across scaffolds.

---

## 🛠 Repository Structure

```text
pyISFinder/
├── pyISFinder.py          # Main CLI executable script
├── pyISFinder/            # Core Python package
│   ├── __init__.py
│   ├── gene_pred.py       # Bacterial ORF prediction (pyrodigal)
│   ├── annotate.py        # HMM search engine (pyhmmer)
│   ├── repeat_finder.py   # Flanking TIR and DR search engine
│   ├── classifier.py      # IS family clustering and classification
│   ├── reporter.py        # TSV, GFF3, and summary table exporter
│   └── visualizer.py      # Matplotlib genomic map visualizer
├── db/                    # Transposase Profile HMMs and Reference Database
│   ├── clusters.faa.hmm
│   └── clusters.single.faa
├── requirements.txt       # Dependencies
├── LICENSE                # MIT License
└── README.md              # Documentation
```

---

## 📜 Acknowledgements & Attribution

`pyISFinder` utilizes transposase profile HMMs and protein reference models derived from:
- **ISEScan** (Xie & Tang, 2017), licensed under the **Apache License 2.0**.
- **ISfinder Database** (Siguier et al., 2006) and **ACLAME Database** (Leplae et al., 2004).

---

## 📄 Citation

If you use `pyISFinder` in your research, please cite this repository or include it in your manuscript's Software Availability section:

```bibtex
@software{pyISFinder2026,
  title = {pyISFinder: A Python pipeline for insertion sequence identification and annotation in bacterial genomes},
  author = {MD Umar},
  year = {2026},
  url = {https://github.com/MdUmar-tech/pyISFinder}
}
```
