import os
import sys
from Bio.Align import PairwiseAligner

def _to_str(val):
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="ignore")
    return str(val)

def search_hmmer_pyhmmer(predicted_genes, hmm_path, evalue_threshold=1e-5):
    """
    Search predicted gene protein sequences against profile HMM database using pyhmmer.
    """
    import pyhmmer
    
    # Prepare digital sequences for pyhmmer
    pyhmmer_seqs = []
    for g in predicted_genes:
        seq_bytes = g["protein"].encode("ascii")
        seq_id = g["gene_id"].encode("ascii")
        digital_seq = pyhmmer.easel.TextSequence(name=seq_id, sequence=seq_bytes).digitize(pyhmmer.easel.Alphabet.amino())
        pyhmmer_seqs.append(digital_seq)
        
    print(f"Loaded {len(pyhmmer_seqs)} protein sequences into pyhmmer digital alphabet.")
    
    hits_list = []
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmms = list(hmm_file)
        print(f"Loaded {len(hmms)} HMM profiles from {os.path.basename(hmm_path)}.")
        
        # Run hmmsearch
        for top_hits in pyhmmer.hmmsearch(hmms, pyhmmer_seqs):
            hmm_name = _to_str(top_hits.query.name)
            
            # Parse IS family from HMM name e.g. "IS110_1|IS110..." or "IS3_4|IS3..."
            is_family = "Unknown"
            if "|" in hmm_name:
                parts = hmm_name.split("|")
                if len(parts) >= 2 and parts[1]:
                    is_family = parts[1]
                else:
                    is_family = parts[0].split("_")[0]
            elif "_" in hmm_name:
                is_family = hmm_name.split("_")[0]
            else:
                is_family = hmm_name

            for hit in top_hits:
                if hit.evalue <= evalue_threshold:
                    gene_id = _to_str(hit.name)
                    score = hit.score
                    evalue = hit.evalue
                    
                    hits_list.append({
                        "gene_id": gene_id,
                        "hmm_name": hmm_name,
                        "is_family": is_family,
                        "evalue": evalue,
                        "score": score
                    })
                    
    return hits_list

def search_fasta_pairwise(predicted_genes, ref_faa_path, min_identity=0.30):
    """
    Fallback alignment engine searching against protein reference database (clusters.single.faa).
    """
    from Bio import SeqIO
    ref_records = list(SeqIO.parse(ref_faa_path, "fasta"))
    print(f"Loaded {len(ref_records)} reference IS proteins for pairwise similarity search.")
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    hits_list = []
    
    # Pre-index k-mers (length 5) for fast candidate filtering
    kmer_idx = {}
    for r in ref_records:
        r_seq = str(r.seq)
        family = "Unknown"
        header = r.id
        if "|" in header:
            parts = header.split("|")
            family = parts[1] if len(parts) >= 2 and parts[1] else parts[0].split("_")[0]
        elif "_" in header:
            family = header.split("_")[0]
            
        for i in range(len(r_seq) - 4):
            km = r_seq[i:i+5]
            if km not in kmer_idx:
                kmer_idx[km] = []
            kmer_idx[km].append((r.id, family, r_seq))
            
    for g in predicted_genes:
        g_seq = g["protein"]
        if len(g_seq) < 30:
            continue
            
        candidate_counts = {}
        for i in range(len(g_seq) - 4):
            km = g_seq[i:i+5]
            if km in kmer_idx:
                for ref_id, fam, r_seq in kmer_idx[km]:
                    candidate_counts[(ref_id, fam, r_seq)] = candidate_counts.get((ref_id, fam, r_seq), 0) + 1
                    
        top_candidates = sorted(candidate_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        
        best_hit = None
        best_score = -1
        
        for (ref_id, fam, r_seq), shared in top_candidates:
            if shared < 3:
                continue
            # Use aligner.score() directly instead of aligner.align()
            score = aligner.score(g_seq, r_seq)
            max_len = max(len(g_seq), len(r_seq))
            identity = (score / (2 * max_len))
            
            if identity >= min_identity and score > best_score:
                best_score = score
                best_hit = {
                    "gene_id": g["gene_id"],
                    "hmm_name": ref_id,
                    "is_family": fam,
                    "evalue": 1e-12,
                    "score": score,
                    "identity": identity
                }
                
        if best_hit:
            hits_list.append(best_hit)
            
    return hits_list

def annotate_transposases(predicted_genes, hmm_path, ref_faa_path=None):
    """
    Annotate transposase proteins among predicted genes.
    """
    hits = []
    try:
        print("Searching transposase HMM profiles using pyhmmer...")
        hits = search_hmmer_pyhmmer(predicted_genes, hmm_path)
        print(f"Identified {len(hits)} transposase hits via pyhmmer.")
        return hits
    except Exception as e:
        print(f"pyhmmer search error ({e}).")
        if ref_faa_path and os.path.exists(ref_faa_path):
            print("Falling back to pairwise reference similarity search...")
            hits = search_fasta_pairwise(predicted_genes, ref_faa_path)
            print(f"Identified {len(hits)} transposase hits via pairwise fallback search.")
            return hits
        else:
            print("ERROR: No valid HMM/reference database available for annotation.")
            return []
