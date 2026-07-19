from Bio.Seq import Seq

def reverse_complement(seq_str):
    return str(Seq(seq_str).reverse_complement())

def find_tir_and_dr(scaffold_seq, g_start, g_end, strand, max_flank=500):
    """
    Search for Terminal Inverted Repeats (TIR) and Direct Repeats (DR/TSD)
    flanking a predicted transposase gene.
    """
    seq_len = len(scaffold_seq)
    
    # 1-indexed to 0-indexed adjustment
    start_idx = max(0, g_start - 1)
    end_idx = min(seq_len, g_end)
    
    # Upstream and downstream flanking windows
    up_start = max(0, start_idx - max_flank)
    up_seq = scaffold_seq[up_start:start_idx].upper()
    
    down_end = min(seq_len, end_idx + max_flank)
    down_seq = scaffold_seq[end_idx:down_end].upper()
    
    down_rc = reverse_complement(down_seq)
    
    best_tir = {
        "found": False,
        "tir_left": "",
        "tir_right": "",
        "tir_length": 0,
        "identity": 0.0,
        "is_start": start_idx + 1,
        "is_end": end_idx
    }
    
    # Search for TIR matches between up_seq (from end backwards) and down_rc (from start forwards)
    # Target length 10 to 40 bp
    min_len = 10
    max_len = min(40, len(up_seq), len(down_rc))
    
    best_match_score = 0
    
    for l_len in range(max_len, min_len - 1, -1):
        for i in range(len(up_seq) - l_len, -1, -1):
            left_candidate = up_seq[i:i+l_len]
            for j in range(0, min(100, len(down_rc) - l_len)):
                right_rc_candidate = down_rc[j:j+l_len]
                
                # Count matching bases
                matches = sum(1 for a, b in zip(left_candidate, right_rc_candidate) if a == b)
                identity = matches / l_len
                
                if identity >= 0.70:
                    score = matches * identity
                    if score > best_match_score:
                        best_match_score = score
                        
                        # Compute genomic boundaries of the IS element including TIR
                        is_start = up_start + i + 1
                        # Right TIR in original downstream sequence
                        right_orig_start = end_idx + (len(down_seq) - (j + l_len))
                        is_end = right_orig_start + l_len
                        
                        best_tir = {
                            "found": True,
                            "tir_left": left_candidate,
                            "tir_right": reverse_complement(right_rc_candidate),
                            "tir_length": l_len,
                            "identity": round(identity * 100, 1),
                            "is_start": is_start,
                            "is_end": is_end
                        }
                        break
            if best_tir["found"]:
                break
        if best_tir["found"]:
            break
            
    # Search for Direct Repeats (DR / TSD) immediately flanking predicted IS start and end
    dr_info = {
        "found": False,
        "dr_seq": "",
        "dr_length": 0
    }
    
    is_s = best_tir["is_start"] - 1
    is_e = best_tir["is_end"]
    
    # DR length 3 to 12 bp
    for dr_len in range(12, 2, -1):
        if is_s - dr_len >= 0 and is_e + dr_len <= seq_len:
            left_dr = scaffold_seq[is_s - dr_len : is_s].upper()
            right_dr = scaffold_seq[is_e : is_e + dr_len].upper()
            
            if left_dr == right_dr and len(left_dr) >= 3:
                dr_info = {
                    "found": True,
                    "dr_seq": left_dr,
                    "dr_length": dr_len
                }
                break
                
    return best_tir, dr_info
