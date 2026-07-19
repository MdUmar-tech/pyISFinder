import os

def export_tsv(is_elements, output_path):
    """
    Write detailed IS element predictions to TSV file.
    """
    headers = [
        "IS_ID", "Scaffold", "IS_Family", "Start", "End", "Strand",
        "Length_bp", "Num_ORFs", "ORFs", "TIR_Found", "TIR_Left",
        "TIR_Right", "TIR_Identity_%", "DR_Found", "DR_Seq", "Score"
    ]
    
    with open(output_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for elem in is_elements:
            row = [
                elem["is_id"],
                elem["scaffold"],
                elem["is_family"],
                str(elem["start"]),
                str(elem["end"]),
                elem["strand"],
                str(elem["length_bp"]),
                str(elem["num_orfs"]),
                elem["orfs"],
                "Yes" if elem["tir_found"] else "No",
                elem["tir_left"] if elem["tir_found"] else "-",
                elem["tir_right"] if elem["tir_found"] else "-",
                f"{elem['tir_identity']:.1f}" if elem["tir_found"] else "-",
                "Yes" if elem["dr_found"] else "No",
                elem["dr_seq"] if elem["dr_found"] else "-",
                f"{elem['best_score']:.1f}"
            ]
            f.write("\t".join(row) + "\n")
    print(f"Exported TSV report: {output_path}")

def export_gff3(is_elements, output_path):
    """
    Write GFF3 browser track for predicted IS elements.
    """
    with open(output_path, "w") as f:
        f.write("##gff-version 3\n")
        for elem in is_elements:
            seqid = elem["scaffold"]
            source = "pyISFinder"
            type_feat = "insertion_sequence"
            start = elem["start"]
            end = elem["end"]
            score = f"{elem['best_score']:.1f}"
            strand = elem["strand"]
            phase = "."
            attributes = f"ID={elem['is_id']};Name={elem['is_family']};family={elem['is_family']};orfs={elem['orfs']}"
            if elem["tir_found"]:
                attributes += f";TIR_left={elem['tir_left']};TIR_right={elem['tir_right']}"
            if elem["dr_found"]:
                attributes += f";DR={elem['dr_seq']}"
                
            f.write(f"{seqid}\t{source}\t{type_feat}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n")
    print(f"Exported GFF3 track: {output_path}")

def export_summary(is_elements, output_path):
    """
    Write manuscript summary table of IS families.
    """
    family_counts = {}
    for elem in is_elements:
        fam = elem["is_family"]
        family_counts[fam] = family_counts.get(fam, 0) + 1
        
    sorted_families = sorted(family_counts.items(), key=lambda x: x[1], reverse=True)
    total_is = len(is_elements)
    
    with open(output_path, "w") as f:
        f.write("=====================================================\n")
        f.write("         pyISFinder Summary Report\n")
        f.write("=====================================================\n\n")
        f.write(f"Total Predicted IS Elements: {total_is}\n\n")
        f.write(f"{'IS Family':<18} {'Count':<10} {'Percentage (%)':<15}\n")
        f.write("-" * 45 + "\n")
        for fam, count in sorted_families:
            pct = (count / total_is * 100) if total_is > 0 else 0
            f.write(f"{fam:<18} {count:<10} {pct:.1f}%\n")
        f.write("-" * 45 + "\n")
    print(f"Exported summary report: {output_path}")
