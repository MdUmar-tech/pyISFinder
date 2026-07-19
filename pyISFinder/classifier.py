def classify_and_cluster_is_elements(annotated_hits, genes_dict, scaffolds_dict, repeat_finder_func):
    """
    Cluster adjacent transposases into full IS elements and assign family taxonomy.
    """
    if not annotated_hits:
        return []
        
    # Group hits by scaffold
    by_scaffold = {}
    for hit in annotated_hits:
        g = genes_dict[hit["gene_id"]]
        scaf = g["scaffold"]
        if scaf not in by_scaffold:
            by_scaffold[scaf] = []
        by_scaffold[scaf].append({**hit, **g})
        
    final_is_list = []
    element_id_counter = 1
    
    for scaf, hits in by_scaffold.items():
        # Sort hits by genomic start coordinate
        hits.sort(key=lambda x: x["start"])
        
        clusters = []
        curr_cluster = [hits[0]]
        
        for h in hits[1:]:
            prev_h = curr_cluster[-1]
            # Merge if distance between ORFs <= 500 bp
            if h["start"] - prev_h["end"] <= 500:
                curr_cluster.append(h)
            else:
                clusters.append(curr_cluster)
                curr_cluster = [h]
        clusters.append(curr_cluster)
        
        scaffold_seq = scaffolds_dict.get(scaf, "")
        
        for cl in clusters:
            elem_id = f"IS_elem_{element_id_counter:04d}"
            element_id_counter += 1
            
            # Element boundaries from ORFs
            min_start = min(x["start"] for x in cl)
            max_end = max(x["end"] for x in cl)
            strand = cl[0]["strand"]
            families = [x["is_family"] for x in cl]
            primary_family = max(set(families), key=families.count)
            
            gene_ids = [x["gene_id"] for x in cl]
            
            # TIR and DR detection
            tir_info, dr_info = repeat_finder_func(scaffold_seq, min_start, max_end, strand)
            
            full_start = tir_info["is_start"] if tir_info["found"] else min_start
            full_end = tir_info["is_end"] if tir_info["found"] else max_end
            full_length = full_end - full_start + 1
            
            final_is_list.append({
                "is_id": elem_id,
                "scaffold": scaf,
                "is_family": primary_family,
                "start": full_start,
                "end": full_end,
                "strand": strand,
                "length_bp": full_length,
                "orfs": ",".join(gene_ids),
                "num_orfs": len(cl),
                "tir_found": tir_info["found"],
                "tir_left": tir_info["tir_left"],
                "tir_right": tir_info["tir_right"],
                "tir_length": tir_info["tir_length"],
                "tir_identity": tir_info.get("identity", 0.0),
                "dr_found": dr_info["found"],
                "dr_seq": dr_info["dr_seq"],
                "dr_length": dr_info["dr_length"],
                "best_evalue": min(x["evalue"] for x in cl),
                "best_score": max(x["score"] for x in cl)
            })
            
    return final_is_list
