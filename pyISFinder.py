#!/usr/bin/env python3
"""
pyISFinder: Pure-Python Insertion Sequence Identification & Annotation Pipeline
"""

import sys
import os
import argparse
from Bio import SeqIO

# Add local directory to python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pyISFinder.gene_pred import predict_genes
from pyISFinder.annotate import annotate_transposases
from pyISFinder.repeat_finder import find_tir_and_dr
from pyISFinder.classifier import classify_and_cluster_is_elements
from pyISFinder.reporter import export_tsv, export_gff3, export_summary
from pyISFinder.visualizer import plot_genome_is_map

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_hmm = os.path.join(script_dir, "db", "clusters.faa.hmm")
    default_ref = os.path.join(script_dir, "db", "clusters.single.faa")

    parser = argparse.ArgumentParser(description="pyISFinder: Insertion Sequence prediction pipeline")
    parser.add_argument("-i", "--input", required=True, help="Path to input bacterial genome FASTA file")
    parser.add_argument("-o", "--outdir", default="results", help="Output directory for reports & figures")
    parser.add_argument("--hmm", default=default_hmm, help="Path to transposase profile HMM database")
    parser.add_argument("--ref", default=default_ref, help="Path to reference transposase FASTA database")
    
    args = parser.parse_args()
    
    fasta_path = args.input
    outdir = args.outdir
    hmm_path = args.hmm
    ref_path = args.ref
    
    if not os.path.exists(fasta_path):
        print(f"Error: Input file '{fasta_path}' does not exist.")
        sys.exit(1)
        
    os.makedirs(outdir, exist_ok=True)
    
    print("=" * 60)
    print("      pyISFinder - Insertion Sequence Annotation Pipeline")
    print("=" * 60)
    print(f"Input Genome:    {fasta_path}")
    print(f"HMM Database:    {hmm_path}")
    print(f"Output Directory: {outdir}")
    print("=" * 60 + "\n")
    
    # 1. Load FASTA scaffolds
    scaffolds_dict = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        scaffolds_dict[rec.id] = str(rec.seq)
    print(f"[Step 1/5] Loaded {len(scaffolds_dict)} scaffolds/contigs from FASTA.")
    
    # 2. Predict ORFs / CDS
    print("\n[Step 2/5] Predicting genes & ORFs...")
    genes = predict_genes(fasta_path)
    genes_dict = {g["gene_id"]: g for g in genes}
    
    # 3. Transposase Homology & Domain Annotation
    print("\n[Step 3/5] Annotating transposases against reference HMM database...")
    annotated_hits = annotate_transposases(genes, hmm_path, ref_path)
    
    # 4. Classify IS elements & Find TIR / DR repeats
    print("\n[Step 4/5] Clustering transposases & searching for flanking TIR / DR repeats...")
    is_elements = classify_and_cluster_is_elements(annotated_hits, genes_dict, scaffolds_dict, find_tir_and_dr)
    print(f"Successfully identified {len(is_elements)} insertion sequence (IS) elements!")
    
    # 5. Export Reports & Genome Map
    print("\n[Step 5/5] Exporting deliverables & genome map...")
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    
    tsv_out = os.path.join(outdir, f"{base_name}_IS_predictions.tsv")
    gff_out = os.path.join(outdir, f"{base_name}_IS.gff3")
    summary_out = os.path.join(outdir, f"{base_name}_IS_summary.txt")
    map_png = os.path.join(outdir, f"{base_name}_IS_genome_map.png")
    map_svg = os.path.join(outdir, f"{base_name}_IS_genome_map.svg")
    
    export_tsv(is_elements, tsv_out)
    export_gff3(is_elements, gff_out)
    export_summary(is_elements, summary_out)
    plot_genome_is_map(is_elements, scaffolds_dict, map_png, map_svg)
    
    print("\n" + "=" * 60)
    print("Pipeline Execution Complete!")
    print(f"Predictions TSV: {tsv_out}")
    print(f"GFF3 Track:     {gff_out}")
    print(f"Summary Table:  {summary_out}")
    print(f"Genome Map:     {map_png}")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()
