import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_genome_is_map(is_elements, scaffolds_dict, output_png, output_svg=None):
    """
    Generate publication-grade genome map of IS elements using matplotlib.
    """
    if not is_elements or not scaffolds_dict:
        print("No IS elements to plot.")
        return
        
    scaffold_lengths = {k: len(v) for k, v in scaffolds_dict.items()}
    # Top 5 longest scaffolds to display
    sorted_scaffolds = sorted(scaffold_lengths.items(), key=lambda x: x[1], reverse=True)[:5]
    
    # Unique colors for IS families
    families = sorted(list(set(e["is_family"] for e in is_elements)))
    cmap = plt.get_cmap("tab10" if len(families) <= 10 else "tab20")
    family_colors = {fam: cmap(i % 20) for i, fam in enumerate(families)}
    
    fig, ax = plt.subplots(figsize=(12, 2.5 * len(sorted_scaffolds)))
    
    y_pos = len(sorted_scaffolds) * 3
    y_ticks = []
    y_labels = []
    
    for scaf, length in sorted_scaffolds:
        # Draw main chromosome/scaffold backbone line
        ax.plot([0, length], [y_pos, y_pos], color="#4A5568", lw=4, zorder=1)
        ax.text(-length * 0.02, y_pos, scaf, va="center", ha="right", fontsize=11, fontweight="bold")
        
        # Filter elements on this scaffold
        scaf_elems = [e for e in is_elements if e["scaffold"] == scaf]
        
        for e in scaf_elems:
            start = e["start"]
            width = max(e["length_bp"], length * 0.005) # Minimum width for visibility
            fam = e["is_family"]
            color = family_colors[fam]
            
            # Position above (+) or below (-) strand
            offset = 0.4 if e["strand"] == "+" else -0.4
            rect_y = y_pos + offset - 0.15
            
            rect = patches.Rectangle(
                (start, rect_y), width, 0.3,
                facecolor=color, edgecolor="black", linewidth=0.8, alpha=0.9, zorder=2
            )
            ax.add_patch(rect)
            
            # Add label for IS family
            ax.text(start + width / 2, rect_y + 0.35 if e["strand"] == "+" else rect_y - 0.25,
                    fam, fontsize=8, ha="center", va="center", fontweight="semibold")
                    
        y_ticks.append(y_pos)
        y_labels.append(scaf)
        y_pos -= 3
        
    ax.set_ylim(-1, len(sorted_scaffolds) * 3 + 1)
    ax.set_xlabel("Genomic Position (bp)", fontsize=12, labelpad=10)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Legend for IS families
    legend_patches = [patches.Patch(color=color, label=fam) for fam, color in family_colors.items()]
    ax.legend(handles=legend_patches, loc="upper right", title="IS Family", bbox_to_anchor=(1.15, 1.0))
    
    plt.title("Genomic Map of Predicted Insertion Sequences (pyISFinder)", fontsize=14, fontweight="bold", pad=20)
    plt.tight_layout()
    
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Saved genome map figure (PNG): {output_png}")
    
    if output_svg:
        plt.savefig(output_svg, format="svg", bbox_inches="tight")
        print(f"Saved genome map figure (SVG): {output_svg}")
        
    plt.close()
