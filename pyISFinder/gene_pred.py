import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

def predict_orfs_pyrodigal(fasta_path):
    """
    Predict ORFs using pyrodigal (Cython wrapper for Prodigal).
    """
    import pyrodigal
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    predicted_genes = []
    
    # pyrodigal v3.x uses GeneFinder or OrfFinder
    if hasattr(pyrodigal, "GeneFinder"):
        orf_finder = pyrodigal.GeneFinder(meta=True)
    elif hasattr(pyrodigal, "OrfFinder"):
        orf_finder = pyrodigal.OrfFinder(meta=True)
    else:
        orf_finder = pyrodigal.OrfFinder()
        
    gene_idx = 1
    for rec in records:
        scaffold_id = rec.id
        sequence = str(rec.seq)
        
        preds = orf_finder.find_genes(sequence)
        for p in preds:
            start = p.begin # 1-indexed / 0-indexed adjustment
            end = p.end
            strand = "+" if p.strand == 1 else "-"
            protein_seq = p.translate()
            dna_seq = p.sequence()
            
            gene_id = f"gene_{gene_idx:05d}"
            gene_idx += 1
            
            predicted_genes.append({
                "gene_id": gene_id,
                "scaffold": scaffold_id,
                "start": start,
                "end": end,
                "strand": strand,
                "protein": protein_seq,
                "dna": dna_seq,
                "length_bp": end - start + 1,
                "length_aa": len(protein_seq)
            })
            
    return predicted_genes

def predict_orfs_fallback(fasta_path, min_len_aa=80):
    """
    Fallback pure-Python ORF finder implementing Translation Table 11 (Bacterial).
    """
    start_codons = {"ATG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    predicted_genes = []
    gene_idx = 1
    
    for rec in records:
        scaffold_id = rec.id
        seq_str = str(rec.seq).upper()
        seq_len = len(seq_str)
        
        for strand, seq in [("+", seq_str), ("-", str(Seq(seq_str).reverse_complement()))]:
            for frame in range(3):
                i = frame
                while i < seq_len - 2:
                    codon = seq[i:i+3]
                    if codon in start_codons:
                        for j in range(i + 3, seq_len - 2, 3):
                            stop_candidate = seq[j:j+3]
                            if stop_candidate in stop_codons:
                                dna_len = (j + 3) - i
                                aa_len = dna_len // 3
                                if aa_len >= min_len_aa:
                                    orf_dna = seq[i:j+3]
                                    protein_seq = str(Seq(orf_dna).translate(table=11, to_stop=True))
                                    
                                    if strand == "+":
                                        genomic_start = i + 1
                                        genomic_end = j + 3
                                    else:
                                        genomic_start = seq_len - (j + 2)
                                        genomic_end = seq_len - i
                                        
                                    gene_id = f"gene_{gene_idx:05d}"
                                    gene_idx += 1
                                    
                                    predicted_genes.append({
                                        "gene_id": gene_id,
                                        "scaffold": scaffold_id,
                                        "start": genomic_start,
                                        "end": genomic_end,
                                        "strand": strand,
                                        "protein": protein_seq,
                                        "dna": orf_dna,
                                        "length_bp": dna_len,
                                        "length_aa": aa_len
                                    })
                                i = j
                                break
                    i += 3
                    
    return predicted_genes

def predict_genes(fasta_path):
    """
    Predict bacterial genes/ORFs from FASTA file using pyrodigal or fallback.
    """
    try:
        print("Running pyrodigal for gene prediction...")
        genes = predict_orfs_pyrodigal(fasta_path)
        print(f"Predicted {len(genes)} genes using pyrodigal.")
        return genes
    except Exception as e:
        print(f"pyrodigal error ({e}). Falling back to pure-Python ORF finder...")
        genes = predict_orfs_fallback(fasta_path)
        print(f"Predicted {len(genes)} genes using fallback ORF finder.")
        return genes
