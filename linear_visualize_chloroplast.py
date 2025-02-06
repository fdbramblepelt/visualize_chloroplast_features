#!/usr/bin/env python3

# Required packages:
# - matplotlib
# - bcbio-gff
# Install using: pip install matplotlib bcbio-gff

# Import necessary libraries

import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import pandas as pd
import csv
import argparse

def parse_gff3_file(file_path):
    """
    Parse the GFF3 file to extract genome features (genes, regions, etc.).
    """
    features = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments and empty lines
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue  # Skip malformed lines
            
            start = parts[3]
            end = parts[4]

            # Skip entries with 'NA' in start or end
            if start == 'NA' or end == 'NA':
                continue

            feature = {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(start),
                "end": int(end),
                "score": parts[5],
                "strand": parts[6],
                "phase": parts[7],
                "attributes": {}
            }
            
            # Parse attributes
            attributes_str = parts[8]
            for attr in attributes_str.split(";"):
                if attr:
                    key_value = attr.split("=")
                    if len(key_value) == 2:
                        key = key_value[0].strip()
                        value = key_value[1].strip()
                        feature["attributes"][key] = value
            
            features.append(feature)
    
    return features

def parse_regions_file(region_file_path):
    """
    Parse the regions file to extract region sizes (LSC, SSC, IRs).
    """
    regions = {"LSC": None, "SSC": None, "IR": []}
    with open(region_file_path, 'r') as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments and empty lines
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue  # Skip malformed lines
            
            region_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            
            if region_type == "lsc":
                regions["LSC"] = (start, end)
            elif region_type == "ssc":
                regions["SSC"] = (start, end)
            elif region_type in ["ira", "irb"]:
                regions["IR"].append((start, end))
    
    return regions

def organize_features(features, regions):
    """
    Organize features into genes and tRNAs.
    """
    genes = []
    tRNAs = []
    
    # Process genes and tRNAs as before
    for feature in features:
        if feature['type'] == "gene":
            genes.append({
                "name": feature['attributes'].get("Name", "Unknown"),
                "start": feature['start'],
                "id": feature['attributes'].get("ID", "Unknown"),
                "end": feature['end'],
                "strand": feature['strand']
            })
        elif feature['type'] == "tRNA":
            tRNAs.append({
                "name": feature['attributes'].get("Name", "Unknown"),
                "id": feature['attributes'].get("ID", "Unknown"),
                "start": feature['start'],
                "end": feature['end'],
                "strand": feature['strand']
            })

    return regions, genes, tRNAs

def load_gene_colors():
    df = pd.read_csv("/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/chloroplast_gene.csv")
    return {row['gene']: row['hex_code'].strip() for _, row in df.iterrows()}

def get_gene_color(gene_name, color_map):
    """
    Find the color for a gene based on partial name matching.
    """
    # Remove 'trn' prefix for tRNA matching to avoid double matches
    for ref_gene, color in color_map.items():
        if gene_name.startswith('trn') and ref_gene.startswith('trn'):
            if gene_name.split('-')[0] == ref_gene.split('-')[0]:
                return color
        # For other genes, check if the gene name is contained in the reference name
        elif gene_name.lower() in ref_gene.lower() or ref_gene.lower() in gene_name.lower():
            return color
    return "#808080"  # Default gray color for unmatched genes

def parse_hairpin_csv(file_path, regions):
    """
    Parse the hairpin CSV file to extract hairpin features.
    """
    hairpins = []
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Skip rows with 'NA' in start_pos or end_pos
            if row['start_pos'] == 'NA' or row['end_pos'] == 'NA':
                continue

            region = row['region'].lower()
            start_pos = int(row['start_pos'])
            end_pos = int(row['end_pos'])
            
            # Adjust start and end positions based on the region
            if region == 'lsc':
                start_pos += regions['LSC'][0] - 1
                end_pos += regions['LSC'][0] - 1
            elif region == 'ssc':
                start_pos += regions['SSC'][0] - 1
                end_pos += regions['SSC'][0] - 1
            elif region == 'irb':
                start_pos += regions['IR'][1][0] - 1
                end_pos += regions['IR'][1][0] - 1
            
            hairpins.append({
                "id": row['ID'],
                "start": start_pos,
                "end": end_pos,
                "sequence": row['sequence'],
                "structure": row['structure'],
                "free_energy": float(row['free_energy']),
                "region": region,
                "genus": row['genus'],
                "species": row['species'],
                "sequence_id": row['sequence_id'],
                "major_clade": row['major_clade'],
                "subfamily": row['subfamily'],
                "tribe": row['tribe'],
                "subtribe": row['subtribe'],
                "photosynthetic_pathway": row['photosynthetic_pathway']
            })
    return hairpins

# Step 3: Plot the genome
def plot_genome(regions, genes, tRNAs, hairpins, genome_length):
    """
    Create a linear plot of the chloroplast genome including hairpins.
    """
    fig, ax = plt.subplots(figsize=(15, 5))
    print(f"Genome length: {genome_length}. Regions: {regions}")
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-1, 1)
    
    # Load gene colors from CSV
    gene_colors = load_gene_colors()
    
    # Plot regions (LSC, SSC, IRs)
    ax.add_patch(Rectangle((regions["LSC"][0], -0.5), regions["LSC"][1] - regions["LSC"][0], 1, color="lightgrey", alpha=0.5, label="LSC"))
    ax.add_patch(Rectangle((regions["SSC"][0], -0.5), regions["SSC"][1] - regions["SSC"][0], 1, color="lightgrey", alpha=0.5, label="SSC"))
    for i, ir in enumerate(regions["IR"]):
        label = "IR" if i == 0 else ""
        ax.add_patch(Rectangle((ir[0], -0.5), ir[1] - ir[0], 1, color="lightgrey", alpha=0.5, label=label))

    count_pos = 0
    count_neg = 0
    # Plot genes with simplified labeling
    for gene in genes:
        if gene['start'] < gene['end'] and gene['start'] >= 0 and gene['end'] <= genome_length:
            y_pos = 0.6 if gene['strand'] == '+' else 0
            height = 0.4
            
            # Get color for this gene
            color = get_gene_color(gene['name'], gene_colors)

            ax.add_patch(Rectangle((gene["start"], y_pos), gene["end"] - gene["start"], height, 
                                 color=color, alpha=1))
            
            if (gene["end"] - gene["start"]) > 50:  # Only label genes above minimum size
                display_name = gene["id"] if gene["name"] == "Unknown" else gene["name"]
                text_y = 1.2 if gene['strand'] == '+' else -0.2
                y_change = 0.12 if gene['strand'] == '+' else -0.12

                if gene['strand'] == '+':
                    count_pos += 1
                else:
                    count_neg += 1

                if count_pos == 2:
                    text_y = text_y + y_change
                elif count_neg == 2:
                    text_y = text_y - y_change
                
                if count_pos == 3:
                    text_y = text_y + (y_change * 2)
                    count_pos = 0
                elif count_neg == 3:
                    text_y = text_y - (y_change * 2)
                    count_neg = 0


                # Draw a line from gene to label
                mid_x = (gene["start"] + gene["end"]) / 2
                ax.plot([mid_x, mid_x], [y_pos + height/2, text_y], linewidth=0.5, color=color)
                
                ax.text(mid_x, text_y, display_name, 
                       ha="center", va="center", fontsize=8, color=color,
                       rotation=30 if gene['strand'] == '+' else -20)

    # Plot tRNAs with simplified labeling
    for tRNA in tRNAs:
        # Check if a gene exists with the same start and end position
        gene_exists = any(gene['start'] == tRNA['start'] and gene['end'] == tRNA['end'] for gene in genes)
        
        if not gene_exists and tRNA['start'] < tRNA['end'] and tRNA['start'] >= 0 and tRNA['end'] <= genome_length:
            y_pos = 0.8 if tRNA['strand'] == '+' else -0.2
            height = 0.4
            
            ax.add_patch(Rectangle((tRNA["start"], y_pos), tRNA["end"] - tRNA["start"], height, 
                                 color="red", alpha=0.5))
            
            if (tRNA["end"] - tRNA["start"]) > 50:  # Only label tRNAs above minimum size
                display_name = tRNA["id"] if tRNA["name"] == "Unknown" else tRNA["name"]
                text_y = 1.8 if tRNA['strand'] == '+' else -0.8
                y_change = 0.16 if tRNA['strand'] == '+' else -0.16
            
                if tRNA['strand'] == '+':
                    count_pos += 1
                else:
                    count_neg += 1

                if count_pos == 2:
                    text_y = text_y + y_change

                if count_pos == 3:
                    text_y = text_y + (y_change * 2)
                    count_pos = 0
                # Draw a line from tRNA to label
                mid_x = (tRNA["start"] + tRNA["end"]) / 2
                ax.plot([mid_x, mid_x], [y_pos + height/2, text_y], linewidth=0.5)
                
                ax.text(mid_x, text_y, display_name, 
                       ha="center", va="center", fontsize=8,
                       rotation=30 if tRNA['strand'] == '+' else -30)

    # Plot hairpins
    for index, hairpin in enumerate(hairpins, start=1):
        y_center = 0.45
        height = 0.65
        ax.add_patch(Rectangle((hairpin["start"], y_center - height / 2), hairpin["end"] - hairpin["start"], height + 0.2, 
                             color="black", alpha=0.8, edgecolor="hotpink", linewidth=1))
        ax.add_patch(Rectangle((hairpin["start"], y_center - height / 2), hairpin["end"] - hairpin["start"], height, 
                             color="hotpink", alpha=1, edgecolor="hotpink", linewidth=1))
        mid_x = (hairpin["start"] + hairpin["end"]) / 2
        y_center = 0.2
        if index % 3 == 0:
            y_center = y_center - 0.3

        elif index % 2 == 0:
            y_center = y_center - 0.15

        ax.text(mid_x, y_center - height, hairpin["id"], 
                ha="center", va="center", fontsize=8, color="hotpink")

    # Add genome summary text
    summary_text = [
        f"LSC: {regions['LSC'][0]} - {regions['LSC'][1]} ({regions['LSC'][1] - regions['LSC'][0]} bp)",
        f"SSC: {regions['SSC'][0]} - {regions['SSC'][1]} ({regions['SSC'][1] - regions['SSC'][0]} bp)"
    ]
    for i, ir in enumerate(regions["IR"]):
        summary_text.append(f"IR{i+1}: {ir[0]} - {ir[1]} ({ir[1] - ir[0]} bp)")
    num_hairpins = len(hairpins)  # Define num_hairpins as the length of the hairpins list
    summary_text.extend([
        f"Total Genes: {len(genes)}",
        f"Total tRNAs: {len(tRNAs)}",
        f"Total Hairpins: {num_hairpins}",
        f"Genome Length: {genome_length}"
    ])

    # Position the summary text in the upper left corner
    summary_y = 3.5
    summary_x = 0
    for line in summary_text:
        ax.text(summary_x, summary_y, line, 
                fontsize=8, va='top', ha='left')
        summary_y -= 0.15
        if summary_y <= 3.1:
            summary_y = 3.5
            summary_x += 25000

    gene_colors = {
        "ATP-synthase": "#4b726e",
        "Calvin-Cycle": "#4d4539",
        "Cytochrome-complex": "#8caba1",
        "Initiation-factor": "#79444a",
        "Large-subunit": "#ae5d40",
        "NADH-dehydrogenase": "#927441",
        "Other": "#847875",
        "Photosystem-I": "#77743b",
        "Photosystem-II": "#b3a555",
        "Small-subunit": "#c77b58",
        "Transcription-and-splicing": "#d1b187",
        "ribosomal-RNA": "#927441",
        "transfer-RNAs": "#ba9158"
    }
    
    # Print legend for each gene color
    for function, color in gene_colors.items():
        print(f"Function: {function}, Color: {color}")

    # Create legend handles for function-color pairs
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=function,
                                  markerfacecolor=color, markersize=10)
                       for function, color in gene_colors.items()]
    # Set plot limits and labels
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-1, 3)  # Adjusted y-limits to accommodate labels
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1), fontsize=8, ncol=5)
    # Define colors for genome regions
    region_colors = {
        "LSC": "#f5d893",
        "SSC": "#e8b26f",
        "IRa": "#b6834c",
        "IRb": "#b6834c"
    }

    # Create legend handles for region-color pairs
    region_legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=region,
                                         markerfacecolor=color, markersize=10)
                              for region, color in region_colors.items()]

    # Add legend for gene colors to the plot
    ax.legend(handles=region_legend_elements + legend_elements, loc='upper right', bbox_to_anchor=(1, 1), fontsize=8, ncol=5)


    # Set the title and subtitle using the first hairpin's metadata
    if hairpins:
        first_hairpin = hairpins[0]
        title = f"{first_hairpin['genus']} {first_hairpin['species']} ({first_hairpin['sequence_id']})"
        subtitle = f"{first_hairpin['major_clade']} | {first_hairpin['subfamily']} | {first_hairpin['tribe']} | {first_hairpin['subtribe']}. " \
                   f"Pathway = {first_hairpin['photosynthetic_pathway']}"
        ax.text(90000, 3.5, title, fontsize=14, va='top', ha='left', style='italic')
        ax.text(0.8, 1.05, subtitle, transform=ax.transAxes, fontsize=10, va='top', ha='center')

    plt.tight_layout()
    plt.show()

def filter_by_sequence_id(data, sequence_id):
    """
    Filter data by sequence ID, allowing partial and case-insensitive matches.
    """
    sequence_id_lower = sequence_id.lower()
    return [item for item in data if sequence_id_lower in item['sequence_id'].lower()]

def main(sequence_id, region_file=None, gff3_file=None, hairpin_file=None):
    """
    Main function to load, parse, and visualize the genome.
    """
    # Use default files if not provided
    if not region_file:
        region_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/all.filtered_chloroplast_genome_region_lengths.gff3"
    if not gff3_file:
        gff3_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/Sorghum_bicolor/Sorghum.bicolor.all_hairpins.csv"
    if not hairpin_file:
        hairpin_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/Sorghum_bicolor/Sorghum.bicolor.all_hairpins.csv"

    # Parse the regions file
    regions = parse_regions_file(region_file)

    # Parse the GFF3 file and filter by sequence ID
    features = parse_gff3_file(gff3_file)
    features = filter_by_sequence_id(features, sequence_id)

    # Organize features into regions, genes, and tRNAs
    regions, genes, tRNAs = organize_features(features, regions)

    # Parse the hairpin file and filter by sequence ID
    hairpins = parse_hairpin_csv(hairpin_file, regions)
    hairpins = filter_by_sequence_id(hairpins, sequence_id)

    # Calculate genome length
    lsc_length = regions["LSC"][1] - regions["LSC"][0]
    ssc_length = regions["SSC"][1] - regions["SSC"][0]
    ir_length = sum(ir[1] - ir[0] for ir in regions["IR"])
    genome_length = lsc_length + ssc_length + ir_length

    # Plot the genome
    plot_genome(regions, genes, tRNAs, hairpins, genome_length)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize chloroplast genome.")
    parser.add_argument("sequence_id", help="Sequence ID to search for in the input files.")
    parser.add_argument("--region_file", help="Path to the region file.", default=None)
    parser.add_argument("--gff3_file", help="Path to the GFF3 file.", default=None)
    parser.add_argument("--hairpin_file", help="Path to the hairpin CSV file.", default=None)

    args = parser.parse_args()
    main(args.sequence_id, args.region_file, args.gff3_file, args.hairpin_file)